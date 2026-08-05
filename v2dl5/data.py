"""Data class holding data store and observations."""

import logging

import numpy as np
from astropy import units as u
from gammapy.data import DataStore

import v2dl5.bti as BTI  # noqa: N812


class Data:
    """
    Data class holding data store and observations.

    Allows to select data from run list or based on
    target coordinates and observation cone.

    Parameters
    ----------
    run_list : str
        Path to run list.
    data_directory : str
        Path to data directory (holding hdu-index.fits.gz and obs-index.fits.gz).
    target : SkyCoord
        Target coordinates.
    obs_cone_radius : float
        observation cone radius (deg).

    """

    def __init__(self, args_dict, target=None):
        """
        Initialize Data object.

        Uses 'run_list' from args_dict if not set to None, otherwise selects data
        according to target coordinates and observation cone.

        """
        self._logger = logging.getLogger(__name__)
        self._observation_cache = {}
        self._bti = args_dict.get("bti")

        self._logger.info(
            "Initializing data object from %s", args_dict["observations"]["datastore"]
        )
        self._data_store = DataStore.from_dir(args_dict["observations"]["datastore"])
        self.target = target
        if args_dict.get("run_list") is None:
            self.runs = self._from_target(
                args_dict["observations"].get("obs_cone_radius", 5.0 * u.deg)
            )
        else:
            self.runs = self._from_run_list(args_dict.get("run_list"))
        self._update_gti(args_dict.get("bti", None))

    def get_data_store(self):
        """Return data store."""
        return self._data_store

    def get_observations(self, reflected_region=True, skip_missing=False):
        """
        Return list of observations.

        Parameters
        ----------
        reflected_region : bool
            Reflected region analysis.
        skip_missing : bool
            Skip missing observations.

        Returns
        -------
        observations : list of `~gammapy.data.Observation`
            List of observations.

        """
        cache_key = (reflected_region, skip_missing)
        if cache_key not in self._observation_cache:
            required_irf = "point-like" if reflected_region else "full-enclosure"
            self._observation_cache[cache_key] = self._data_store.get_observations(
                self.runs,
                required_irf=required_irf,
                skip_missing=skip_missing,
            )
            self._apply_bti(self._observation_cache[cache_key])
        return self._observation_cache[cache_key]

    def _from_run_list(self, run_list):
        """
        Read run list from file and select data.

        Parameters
        ----------
        run_list : str
            Path to run list.

        """
        if run_list is None:
            return None

        try:
            _runs = np.loadtxt(run_list, dtype=int, usecols=0)
        except OSError:
            self._logger.error("Run list %s not found.", run_list)
            raise

        _runs = [_runs] if _runs.ndim == 0 else np.ndarray.tolist(_runs)

        self._logger.info("Reading run list with %d observations from %s", len(_runs), run_list)
        if len(_runs) == 0:
            self._logger.error("Run list is empty.")
            raise ValueError
        return _runs

    def _from_target(self, obs_cone_radius):
        """
        Select data based on target coordinates and observation cone.

        Parameters
        ----------
        obs_cone_radius : float
            observation cone radius (deg).

        """
        if self.target is None:
            raise ValueError("A target is required for target-based data selection")
        observations = self._data_store.obs_table
        obs_cone_radius = u.Quantity(obs_cone_radius).to(u.deg)
        mask = self.target.separation(observations.pointing_radec) < obs_cone_radius
        _runs = observations[mask]["OBS_ID"].data

        self._logger.info(
            "Selecting %d runs from observation cone around %s", len(_runs), self.target
        )
        if len(_runs) == 0:
            raise ValueError(f"No observations found within {obs_cone_radius} of {self.target}")
        return _runs

    def get_on_region_radius(self):
        """
        Return on region radius.

        Simplest case. Ignores possible energy and offset dependence.

        """
        observations = self.get_observations()
        try:
            rad_max = {obs.rad_max.data[0][0] for obs in observations}
        except IndexError:
            self._logger.error("On region radius not found in observations.")
            raise

        if len(rad_max) > 1:
            self._logger.error("On region radius not the same for all observations.")
            raise ValueError

        on_region = rad_max.pop() * u.deg
        self._logger.info(f"On region radius: {on_region}")

        return on_region

    def get_max_wobble_distance(self, fov=3.5 * u.deg):
        """
        Return maximum distance from target position.

        Add if necessary the telescope field of view.

        Parameters
        ----------
        fov : astropy.units.Quantity
            Telescope field of view.

        Returns
        -------
        max_offset : astropy.units.Quantity
            Maximum offset (radius of FoV).

        """
        woff = np.array(
            [
                self.target.separation(obs.pointing.get_icrs()).degree
                for obs in self.get_observations()
            ]
        )
        return np.max(woff) * u.deg + fov / 2.0

    def _update_gti(self, bti):
        """
        Update good time intervals by removing bad time intervals.

        Parameters
        ----------
        bti : list of dict
            List of bad time intervals
            Given us {"run": run, "bti_start": start, "bti_length": length}

        """
        if bti is not None:
            self._bti = bti

        for observations in self._observation_cache.values():
            self._apply_bti(observations)

    def _apply_bti(self, observations):
        """Apply configured bad-time intervals to an observation collection."""
        if self._bti is None:
            return

        for obs in observations:
            bti_pairs = [
                (item["bti_start"], item["bti_start"] + item["bti_length"])
                for item in self._bti
                if item["run"] == obs.obs_id
            ]
            if len(bti_pairs) == 0:
                self._logger.debug("No BTI found for %s", obs.obs_id)
                continue
            self._logger.debug("Updating GTI for %s with %s", obs.obs_id, bti_pairs)
            obs.gti = BTI.BTI(obs).update_gti(bti_pairs)
