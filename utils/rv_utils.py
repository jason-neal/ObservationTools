import numpy as np

# #######################################################
# Functions for RV calculations
# #######################################################
def true_anomaly(ma, ecc, niterationmax=10000):
    """Compute the true anomaly using the Newton-Raphson method.

    Parameters
    ----------
    ma: array-like
        Mean anomaly.
    ecc: float
        Orbital eccentricity.
    niterationmax: int
        Maximum number of iterations for N-R method.

    Returns
    -------
    ta: array-like
        True anomaly

    Notes
    -----
    Adapted from Rodrigo Diaz.

    """
    if not isinstance(ma, (int, float)):
        ea = ma
    else:
        ea = np.array([ma, ])

    if isinstance(ea, list):
        raise TypeError("Unsupported type 'list', input a numpy array or an int/float.")
    if len(ea) == 0:
        raise ValueError("A empty array was given.")

    # Initialise at ea0 = ma
    niteration = 0
    ea0 = ma

    while np.linalg.norm(ea - ea0, ord=1) > 1e-5 or niteration == 0:
        ea0 = ea

        ff = ea - ecc * np.sin(ea) - ma   # Function
        dff = 1 - ecc * np.cos(ea)        # Derivative

        # Use Newton method
        ea = ea0 - ff / dff

        # Increase iteration number; if above limit, break with exception.
        niteration += 1
        if niteration >= niterationmax:
            raise RuntimeError('Eccentric anomaly computation'
                               'not converged.')

    # Compute true anomaly from eccentric anomaly
    return 2. * np.arctan2(np.sqrt(1. + ecc) * np.sin(ea / 2.),
                           np.sqrt(1. - ecc) * np.cos(ea / 2.))


def mean_anomaly(times, t0, period):
    """Calculate mean anomaly using period, tau and a time value.

    Parameters
    ----------
    times: array-like
        Times to compute mean anomaly.
    t0: float
        Time of periastron passage. (Julian days)
    period: float
        Period of orbit.

    Returns
    -------
    ma: array-like
        Mean anomaly.

    """
    if not isinstance(times, (int, float)):
        times = times
    else:
        times = np.array(times)
    return 2 * np.pi * (times - t0) / period


def radial_velocity(gamma, k, ta, omega, ecc):
    """Radial velocity equation.

    Parameters
    ----------
    gamma: float
        Mean RV motion of system.
    k: float
        RV amplitude.
    ta: array-like, float
        True anomaly.
    omega: float
        Argument of periastron. (radians)
    ecc: float
        Eccentricity of orbit.

    Returns
    -------
    RV: array-like, float
        Radial veloctity values

    Notes
    -----
    RV = gamma + k *(np.cos(ta + omega) + ecc * np.cos(omega)).

    """
    # Calculate radial velocity of star
    return gamma + k * (np.cos(ta + omega) + ecc * np.cos(omega))


# RV calculation done in python (for when ajplanet is not available)
def rv_curve_py(times, gamma, k, omega, ecc, t0, period):
    """Generate values for Radial Velocity curve.

    Parameters
    times: array-like, float
        Times to calcualte radial velocity values(Julian days))
    gamma: float
        Mean RV offset value.
    k: float
        Radial velocity amplitude.
    omega: float
        Argument of periastron (radians).
    ecc: float
        Eccentricity.
    t0: float
        Time of periastron passage.
    period: float
        Period of orbit.

    Returns
    -------
    rv: array-like
        Radial velocity values evaulated at the given times.

    """
    ma = mean_anomaly(times, t0, period)
    ta = true_anomaly(ma, ecc)
    rv = radial_velocity(gamma, k, ta, omega, ecc)
    return rv


def RV_from_params(t, params, ignore_mean=False, companion=False):
    """Get radial velocity values at given times using the orbital parameters.

    Parameters
    ----------
    t: array-like
        The times at which to calculate the RV.
    params: dict
        Orbtial parameters required for RV.(mean_val, k1, omega, e, tau, period, optinal [k2]).
    ignore_mean: bool
        Ignore the average radial velocity motion of the system
    companion: bool
        Calculate RV for companion instead.

    Notes:
    1) Omega should be given in degrees. This function converts it to radians.
    2) The units of mean_val and k1,k2 should be the same e.g. both km/s

    Returns:
    rvs: array-like
        The radial velocity values evaluated at the provided times.

    """
    if not isinstance(t, np.ndarray):
        t = np.array(t)

    # Select needed entries from dict to calcualte rv
    if companion:
        param_list = [params["mean_val"], params["k2"], params["omega"],
                      params["eccentricity"], params["tau"], params["period"]]
    else:
        param_list = [params["mean_val"], params["k1"], params["omega"],
                      params["eccentricity"], params["tau"], params["period"]]

    param_list[2] = np.deg2rad(param_list[2])   # convert omega to radians

    if not ignore_mean:
        # Set the mean rv veloctiy to zero
        param_list[0] = 0

    # if use_ajplanet:   # Use ajplanet if available
    #     rvs = pl_rv_array(t, *param_list[:])  # *unpacks parameters from list
    # else:              # Use python version
    rvs = rv_curve_py(t, *param_list[:])  # *unpacks parameters from list

    return rvs