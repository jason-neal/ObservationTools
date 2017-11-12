import copy

import numpy as np

from observationtools.utils.parse import parse_paramfile


class RV(object):
    """

    Notes
    -----
    1) Omega should be given in degrees. This function converts it to radians.
    2) The units of mean_val and k1, k2 should be the same e.g.both km / s
    3) Tau should be the julian date, and the period given in days.
    """

    def __init__(self, semi_amp=0.0, period=0.0, ecc=0.0,
                 tau=0.0, gamma=0.0, omega=0.0, **other_params):
        self.semi_amp = semi_amp
        self.period = period
        self.ecc = ecc
        self.tau = tau
        self.gamma = gamma
        self.omega = omega
        self._params = self.orbit_dict()

        self.ignore_mean = other_params.get("ignore_mean", False)
        if other_params is not None:
            self._params.update(other_params)

    def __repr__(self):
        other_params = ""
        for key, value in self._params.items():
            if key not in ["k1", "eccentricity", "period", "mean_val", "tau", "omega"]:
                other_params += ", {}={}".format(key, value)

        return "RV(semi_amp={1}, period={2}, ecc={3}, tau={4}, omega={5}, gamma={6}{7})".format(
            self.__class__, self.semi_amp, self.period, self.ecc, self.tau, self.omega, self.gamma, other_params)

    def orbit_dict(self):
        return {"k1": self.semi_amp, "period": self.period, "eccentricity": self.ecc,
                "tau": self.tau, "mean_val": self.gamma, "omega": self.omega}

    def to_dict(self):
        self._params.update(self.orbit_dict())
        self._params.update({"ignore_mean": self.ignore_mean})
        return self._params

    @classmethod
    def from_dict(cls, params):
        other_params = params.copy()
        for par in ["k1", "eccentricity", "period", "mean_val", "tau", "omega"]:
            other_params.pop(par)
        return cls(semi_amp=params["k1"], period=params["period"], ecc=params["eccentricity"],
                   tau=params["tau"], gamma=params["mean_val"], omega=params["omega"], **other_params)

    @classmethod
    def from_file(cls, filename):
        """Parameters in key = val\n text file."""
        param_dict = parse_paramfile(filename)
        return cls.from_dict(param_dict)

    def create_companion(self, mass_ratio=None):
        """Create companion RV object.

        It has 3 ways to determine the amplitude of the companion in order of priority:
            1: using a mass_ratio m1/m2 passed into the method.
            2: If "k2" is already a parameter use that.
            3: Calculate the mass ratio from m1 and m2. (In same units)

        Switches k1 and k2 and m1 an m2 parameters. (m1 refers to self, while m2 the other body in orbit.)

        Inputs
        ------
        mass_ratio: float
            m_1 / m_2

        Returns
        ------
        companion: RV
            Rv object for the companion.
        """
        params = copy.copy(self._params)
        if mass_ratio is not None:
            k2 = -params["k1"] * mass_ratio
            params["ratio_flag"] = True
        elif params.get("k2") is not None:
            k2 = params.get("k2")
            params["k2_flag"] = True
        else:
            # Make from masses in parameters
            M1 = params.get("m1")
            M2 = params.get("m2")
            if (M1 is None) or (M2 is None):
                print("params =", params)
                raise ValueError("A mass parameter (m1 or m2) was not provided")
            else:
                # M1_jup = M1 * (M_sun / M_jup).value
                mass_ratio = M1 / M2  # assuming m1 and m2 have same units
            k2 = -params["k1"] * mass_ratio
            params["m1"], params["m2"] = params.get("m2"), params.get("m1")

        params["k2"], params["k2_old"] = k2, params.get("k2")

        params["k1"], params["k2"] = params["k2"], params["k1"]  # Switch to Companion
        return RV.from_dict(params)

    @property
    def ignore_mean(self):
        try:
            return self._ignore_mean
        except AttributeError:
            return False

    @ignore_mean.setter
    def ignore_mean(self, value=None):
        if value is None:
            val = self._params.get("ignore_mean", False)
        self._ignore_mean = value

    def rv_at_phase(self, phase):
        t = phase * self.period + self.tau
        return self.rv_at_times(t)

    def rv_at_times(self, t):
        """Evaluate RV at the provided times."""
        true_anomaly = self.true_anomaly(self.mean_anomaly(t, self.tau, self.period), self.ecc)
        return self.radial_velocity(self.gamma, self.semi_amp,
                                    true_anomaly, self.omega, self.ecc)

    def rv_full_phase(self, center=0, points=100):
        """Return RV curve evaluated one full phase."""
        phase = np.linspace(0, 1, points) + center
        return self.rv_at_phase(phase)

    def max_amp(self):
        amp_1 = self.semi_amp * (1 + self.ecc * np.cos(self.omega * np.pi / 180))
        amp_2 = self.semi_amp * (1 - self.ecc * np.cos(self.omega * np.pi / 180))
        return np.max([np.abs(amp_1), np.abs(amp_2)])

    @staticmethod
    def true_anomaly(ma, ecc, niterationmax=10000):
        # type: (Any, float, int) -> Any
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
            ea = np.array([ma, ], dtype=np.float)

        if isinstance(ea, list):
            raise TypeError("Unsupported type 'list', input a numpy array or an int/float.")
        if len(ea) == 0:
            raise ValueError("A empty array was given.")

        # Initialise at ea0 = ma
        niteration = 0
        ea0 = ma

        while np.linalg.norm(ea - ea0, ord=1) > 1e-5 or niteration == 0:
            ea0 = ea

            ff = ea - ecc * np.sin(ea) - ma  # Function
            dff = 1 - ecc * np.cos(ea)  # Derivative

            # Use Newton method
            ea = ea0 - ff / dff

            # Increase iteration number; if above limit, break with exception.
            niteration += 1
            if niteration >= niterationmax:
                raise RuntimeError('Eccentric anomaly computation '
                                   'not converged.')

        # Compute true anomaly from eccentric anomaly
        return 2. * np.arctan2(np.sqrt(1. + ecc) * np.sin(ea / 2.),
                               np.sqrt(1. - ecc) * np.cos(ea / 2.))

    @staticmethod
    def mean_anomaly(times, t0, period):
        # type: (Any, float, float) -> Any
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

    @staticmethod
    def radial_velocity(gamma, k, ta, omega, ecc):
        # type: (float, float, Any, float, float, float) -> Any
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
            Radial velocity values

        Notes
        -----
        RV = gamma + k *(np.cos(ta + omega) + ecc * np.cos(omega)).

        """
        # Calculate radial velocity of star
        return gamma + k * (np.cos(ta + omega) + ecc * np.cos(omega))

    def __eq__(self, other):
        # all properties the same
        checks = [self.semi_amp == other.semi_amp,
                  self.period == other.period,
                  self.omega == other.omega,
                  self.tau == other.tau,
                  self.ecc == other.ecc,
                  self.gamma == other.gamma]

        return all(checks)

    def __ne__(self, other):
        return not self.__eq__(other)


