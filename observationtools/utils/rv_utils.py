import copy
import datetime
import logging
from typing import Any

from astropy.constants import M_jup, M_sun

from observationtools.rv.rv import RV
from observationtools.utils.parse import parse_obslist


# TODO: Replace "Any" with numpy type hint when available


# RV calculation done in python (for when ajplanet is not available)
def rv_curve_py(times, gamma, k, omega, ecc, t0, period):
    # type: (Any, float, float, float, float, float, float) -> Any
    """Generate values for Radial Velocity curve.

    Parameters
    times: array-like, float
        Times to calculate radial velocity values(Julian days))
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
        Radial velocity values evaluated at the given times.

    """
    ma = RV.mean_anomaly(times, t0, period)
    ta = RV.true_anomaly(ma, ecc)
    rv = RV.radial_velocity(gamma, k, ta, omega, ecc)
    return rv


class JulianDate(object):
    """Handle julian dates."""
    julian_epoch_dt = datetime.datetime(2000, 1, 1, 12)  # noon
    julian_epoch_jd = datetime.timedelta(2451545)  # julian epoch in julian dates
    reduce_jd = 2400000
    strformat = "%Y-%m-%d %H:%M:%S"
    strformat2 = "%Y-%m-%d"

    def __init__(self, jd, reduced=False):
        self.jd = jd
        self.reduced = reduced

    @classmethod
    def from_datetime(cls, dt, reduced=False):
        # type: (Any, bool) -> JulianDate
        """Convert from a datetime to a jd object.

        Test against pyehem.julian_date()

        Parameters
        ----------
        dt: datetime object
            Datetime for date to calculate jd.
        reduced: bool
            Return reduced JD, (JD-2400000)

        Returns
        -------
        jd: JulianDate
            JulianDate object.
        Inspiration from https://stackoverflow.com/questions/13943062/
        """
        jd = dt - cls.julian_epoch_dt + cls.julian_epoch_jd
        jd = float(jd.total_seconds()) / (24 * 60 * 60)  # Turn timedelta into a float
        if reduced:
            jd -= cls.reduce_jd
        return cls(jd, reduced=reduced)

    def to_datetime(self):
        # type: None -> datetime.datetime
        """ Return JulianDate as a datetime.datetime object.

        Returns
        -------
        dt: datetime object
            Datetime of julian date.
        Inspiration from https://stackoverflow.com/questions/13943062/
        """
        print("self.jd", self.jd)
        if self.reduced:
            _jd = self.jd + self.reduce_jd
        else:
            _jd = self.jd
        print("_jd", _jd)
        dt = datetime.timedelta(_jd) + self.julian_epoch_dt - self.julian_epoch_jd
        return dt

    @classmethod
    def now(cls):
        return cls.from_datetime(datetime.datetime.now())

    @classmethod
    def from_str(cls, time_str, format=strformat):
        """Return JulianDate from a time string.

        Inputs
        ------
        time_str: str
        format: str
            Format of time string.

        Returns
        -------
        dt: datetime object
            Datetime of julian date.
        """
        if format is None:
            format = cls.strformat
        try:
            dt = datetime.datetime.strptime(time_str, format)
        except ValueError:
            try:
                dt = datetime.datetime.strptime(time_str, cls.strformat)
            except ValueError:
                dt = datetime.datetime.strptime(time_str, cls.strformat2)
        return cls.from_datetime(dt)

    def to_str(self, format=None):
        """Return date string from a JulianDate.

        Input
        -----
        format: str
             String datetime format.

        Returns
        -------
        datesting: str
            String with date.
        """
        if format is None:
            format = self.strformat
        dt = self.to_datetime()
        return dt.strftime(format)

    def reduce(self):
        if not self.reduced:
            self.jd -= self.reduce_jd
            self.reduced = True


def strtimes2jd(obs_times, reduced=False, format=None):
    # type: (Union[str, List[str]], bool, Union[None, str]) -> List[float]
    """Convenience function for convert str times to reduced JD.
    If reduced=True returns JD-2400000
    """
    change_flag = False
    if obs_times is not None:
        print("obs times", obs_times)

        if isinstance(obs_times, str):
            obs_times = [obs_times]
            change_flag = True

        if isinstance(obs_times, (list, tuple)):
            jds = []
            for obs in obs_times:
                jd = JulianDate.from_str(obs, format)
                print(jd)
                if reduced:
                    jd.reduce()
                print(jd)
                jds.append(jd.jd)
        print("obs jd times", jds)
        if change_flag == 1:
            jds = jds[0]
        return jds
    else:
        return None


def join_times(obs_times=None, obs_list=None):
    # type: (List[str], str) -> List[str]
    """Combine observation dates and turn to jd.

    Parameters
    ----------
    obs_times: list of str or None
        List of dates entered at command line.
    obs_list: str or None
        Filename to observation list.

    Returns
    -------
    obs_times: list of str
        Combined list of date strings.

    """
    if obs_times is None:
        obs_times = []

    if obs_list is None:
        obs_list = []
    else:
        obs_list = parse_obslist(obs_list)

    logging.debug("obs list = {}", format(obs_list))
    obs_times = obs_times + obs_list

    logging.debug("obs_times = {}".format(obs_times))
    if not obs_times:  # An empty list is "Falsely"
        return None
    else:
        return obs_times


def prepare_mass_params(params, only_msini=True):
    """Update parameter dictionary to set m1 and m2 if not given."""
    if params.get("m1") is None:
        params["m1"] = params["m_sun"]  # solar mass

    # Convert m_sun to jupyter masses
    params["m1"] = params["m1"] * M_sun / M_jup  # jupyter mass

    if params.get("m2") is None:
        params["m2"] = params["m_sini"] if only_msini else params["m_true"]
        # should give a key error if the correct mass not given

    params["msini_flag"] = only_msini

    params["k2_flag"] = False if params.get("k2") is None else True

    return params


def check_core_parameters(params):
    """Test core rv parameters."""
    for key in ["name", "k1", "eccentricity", "omega", "tau", "period"]:
        if key not in params.keys():
            raise ValueError("A core parameter was not provided in the param file, '{}'".format(key))

    if "mean_val" not in params.keys():
        logging.info("mean_val parameter was not provided so set to 0 km/s")
        params["mean_val"] = 0.0
    elif params["mean_val"] == "":
        logging.info("mean_val parameter was blank so set to 0 km/s")
        params["mean_val"] = 0.0

    return params


def generate_companion_label(companion):
    msini_flag = companion._params.get("msini_flag", False)
    k2_flag = companion._params.get("k2_flag", False)
    ratio_flag = companion._params.get("ratio_flag", False)

    if not msini_flag and not k2_flag and not ratio_flag:
        label = "M2 Companion"
    elif msini_flag and not k2_flag and not ratio_flag:
        label = "M2sini Companion"
    elif k2_flag and not ratio_flag:
        label = "Given k2 Companion"
    else:
        label = "Mass ratio Companion"

    return label
