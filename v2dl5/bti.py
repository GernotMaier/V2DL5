"""Bad time interval definition."""

import logging

import astropy.units as u
import numpy as np
from gammapy.data import GTI


class BTI:
    """
    Bad time intervals (BTI).

    Add bad time intervals to gammapy GTI objects.

    """

    def __init__(self, obs=None):
        self._logger = logging.getLogger(__name__)

        self.gti = obs.gti
        self.tstart = obs.tstart
        self.tstop = obs.tstop

    def max_time_interval(self):
        """Return time between beginning of first and end of last GTI."""
        duration = round((np.max(self.gti.met_stop) - np.min(self.gti.met_start)).value)
        self._logger.info("Max time interval [s]: %d", duration)

        return duration

    def _create_on_time_list(self, bti):
        """
        Create list of on and off times in 1 s bin.

        Parameters
        ----------
        bti : list
            List of bad time intervals
            (given as a list of start and stop times relative to the observation start time).

        Returns
        -------
        on_time_s : list
            List of on and off times in 1 s bins.

        """
        on_time_s = [0] * (self.max_time_interval() + 1)

        start_times = (self.gti.table["START"] - self.tstart).to_value("sec")
        end_times = (self.gti.table["STOP"] - self.tstart).to_value("sec")
        self._logger.info(f"start times: {start_times}")
        self._logger.info(f"stop times: {end_times}")

        for start, end in zip(start_times, end_times):
            i_1 = round(start)
            i_2 = round(end)
            on_time_s[i_1 : i_2 + 1] = [1] * (i_2 - i_1 + 1)

        for start, end in bti:
            start, end = max(start, 0), max(end, 0)
            if round(start) < len(on_time_s):
                i_1 = round(start)
                i_2 = round(end) if round(end) < len(on_time_s) else len(on_time_s) - 1
                on_time_s[i_1 : i_2 + 1] = [0] * (i_2 - i_1 + 1)

        return on_time_s

    def _extract_gti_start_stop(self, on_time_s):
        """
        Extract start and stop times from list of on times.

        Parameters
        ----------
        on_time_s : list
            List of on and off times in 1 s bins.

        Returns
        -------
        _gti_start : list
            List of start times.
        _gti_end : list
            List of stop times.

        """
        _gti_start = []
        _gti_end = []
        _duration = self.max_time_interval()

        if on_time_s[0] != 0:
            _gti_start.append(0 * u.s)

        for i in range(1, _duration - 1):
            if on_time_s[i] == 0 and on_time_s[i - 1] == 1:
                _gti_end.append((i - 1) * u.s)
            if on_time_s[i] == 0 and on_time_s[i + 1] == 1:
                _gti_start.append((i + 1) * u.s)

        if on_time_s[-1] != 0:
            _gti_end.append(_duration * u.s)

        self._logger.info(f"start times (updated): {_gti_start}")
        self._logger.info(f"stop times (updated): {_gti_end}")

        _met_tstart = self._met_tstart()
        _gti_start = [start + _met_tstart for start in _gti_start]
        _gti_end = [end + _met_tstart for end in _gti_end]

        return _gti_start, _gti_end

    def _met_tstart(self):
        """
        Return start time in seconds since MET.

        Returns
        -------
        tstart : astropy.units.Quantity
            Start time in seconds since MET.

        """
        self._logger.info(f'Seconds since MET {(self.tstart - self.gti.time_ref).to_value("sec")}')

        return (self.tstart - self.gti.time_ref).to_value("sec") * u.s

    def update_gti(self, bti):
        """
        Update good time intervals by removing bad time intervals.

        Parameters
        ----------
        bti : list
            List of bad time intervals
            (given as a list of start and stop times relative to the observation start time).

        Returns
        -------
        gti : gammapy.data.GTI
            Updated good time intervals.

        """
        self._logger.info(f"BTI: {bti}")

        if bti is None:
            return self.gti

        on_time_s = self._create_on_time_list(bti)
        gti_start, gti_end = self._extract_gti_start_stop(on_time_s)

        return GTI.create(gti_start, gti_end, self.gti.time_ref)
