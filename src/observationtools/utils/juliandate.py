import datetime


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