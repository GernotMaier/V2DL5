"""Sky regions definition."""

import logging
from importlib.resources import files
from pathlib import Path

from astropy import units as u
from astropy.coordinates import Angle, SkyCoord, name_resolve
from astropy.io import fits
from astropy.table import Table
from gammapy.maps import WcsGeom
from regions import CircleSkyRegion


class SkyRegions:
    """
    Define sky regions required for the analysis.

    Sky regions include:

    - target coordinates
    - on region
    - exclusion regions

    Parameters
    ----------
    args_dict : dict
        Dictionary of configuration arguments.

    """

    def __init__(self, args_dict=None):
        self._logger = logging.getLogger(__name__)

        self.target = self.get_target(sky_coord=args_dict["on_region"])
        self.on_region = self.define_on_region(on_region_dict=args_dict["on_region"])
        self.exclusion_mask = None

    def get_target(self, sky_coord=None, print_target_info=True):
        """
        Define SkyCoord object for the target.

        Reads target coordinates from Simbad if target name is given
        (use coordinates, if given)

        Parameters
        ----------
        sky_coord : dict
            sky_coord dictionary.

        print_target_info: bool
            Print info about target.

        """
        target = None
        lon, lat = sky_coord.get("lon", None), sky_coord.get("lat", None)
        if lon is not None and lat is not None:
            if sky_coord.get("frame", None) == "icrs":
                target = SkyCoord(
                    ra=lon,
                    dec=lat,
                    frame=sky_coord.get("frame", None),
                )
            elif sky_coord.get("frame", None) == "galactic":
                target = SkyCoord(
                    l=lon,
                    b=lat,
                    frame=sky_coord.get("frame", None),
                )
            else:
                raise ValueError("Unsupported coordinate frame")
            self._logger.debug(f"Target coordinates from configuration: {target}")
        else:
            try:
                target = SkyCoord.from_name(sky_coord["target"])
                self._logger.debug("Target %s found in Simbad.", sky_coord["target"])
            except name_resolve.NameResolveError:
                self._logger.error('Target "%s" not found in Simbad.', sky_coord["target"])
                raise

        if print_target_info:
            self._logger.info("Target name: %s", sky_coord.get("target", "target name not given"))
            self._logger.info("Target coordinates: %s", target)

        return target

    def define_on_region(self, on_region_dict=None):
        """
        Define CircleSkyRegion object for the on region.

        Parameters
        ----------
        on_region : dict
            on_region dictionary.

        Returns
        -------
        on_region : CircleSkyRegion
            on region.

        """
        self.on_region = CircleSkyRegion(
            center=self.get_target(sky_coord=on_region_dict, print_target_info=False),
            radius=Angle(on_region_dict.get("radius", 0.5 * u.deg)),
        )

        self._logger.info(f"On region: {self.on_region}")

        return self.on_region

    def get_exclusion_mask(self, args_dict, max_wobble_distance=None):
        """
        Define mask for the exclusion regions.

        Parameters
        ----------
        args_dict: dict
            Dictionary of configuration arguments.
        max_wobble_distance: Angle
            maximum wobble distance (with FOV radius added).

        """
        exclusion_regions = []

        # on region
        on_region_dict = args_dict["on_region"]
        on_region_exclusion_radius = args_dict["datasets"]["exclusion_region"]["on_radius"]
        if on_region_dict is not None and on_region_exclusion_radius is not None:
            exclusion_regions.append(
                CircleSkyRegion(
                    center=self.get_target(sky_coord=on_region_dict, print_target_info=False),
                    radius=Angle(on_region_exclusion_radius),
                )
            )
            self._logger.info(f"On region exclusion: {exclusion_regions[-1]}")

        # bright star exclusion
        exclusion_regions.extend(
            self._read_bright_star_catalogue(
                exclusion_region_dict=args_dict["datasets"]["exclusion_region"],
                max_wobble_distance=max_wobble_distance,
            )
        )

        # exclusion mask
        geom = WcsGeom.create(
            npix=(150, 150), binsz=0.05, skydir=self.target.icrs, proj="TAN", frame="icrs"
        )

        self._logger.info("Number of exclusion regions: %d", len(exclusion_regions))
        return ~geom.region_mask(exclusion_regions)

    def _read_bright_star_catalogue(self, exclusion_region_dict=None, max_wobble_distance=None):
        """
        Read bright star catalogue from file.

        Catalogue files are expected to be in the v2dl5/data directory and are distributed as part
        of the v2dl5 package.

        Parameters
        ----------
        exclusion_region_dict: dict
            Dictionary of configuration arguments.
        max_wobble_distance: Angle
            maximum wobble distance (with FOV radius added).

        Returns
        -------
        List of CircleSkyRegion
            List exclusion regions due to bright stars.

        """
        _exclusion_regions = []
        if exclusion_region_dict["star_file"] is None:
            return _exclusion_regions

        self._logger.info(
            f"Reading bright star catalogue from {exclusion_region_dict['star_file']}"
            f" (max wobble distance: {max_wobble_distance})"
        )
        configured_path = Path(exclusion_region_dict["star_file"])
        star_file = (
            configured_path
            if configured_path.is_absolute() or configured_path.exists()
            else files("v2dl5").joinpath("data", configured_path.name)
        )
        with fits.open(star_file) as hip:
            catalogue = Table(hip[1].data)
        catalogue = catalogue[
            catalogue["Vmag"] + catalogue["B-V"] < exclusion_region_dict["magnitude_B"]
        ]

        catalogue["angular_separation"] = self.target.separation(
            SkyCoord(catalogue["_RA_icrs"], catalogue["_DE_icrs"], unit=u.deg)
        )

        catalogue = catalogue[catalogue["angular_separation"] < u.Quantity(max_wobble_distance)]
        self._logger.info(
            "Number of stars in the catalogue passing cuts on magnitude and FOV: %d", len(catalogue)
        )

        for row in catalogue:
            _exclusion_regions.append(
                CircleSkyRegion(
                    center=SkyCoord(row["_RA_icrs"], row["_DE_icrs"], unit=u.deg),
                    radius=Angle(exclusion_region_dict["star_exclusion_radius"]),
                )
            )

        return _exclusion_regions

    def update_regions(self, args_dict, on_region_radius=None, max_wobble_distance=None):
        """
        Update on and exclusion region definitions.

        Parameters
        ----------
        args_dict: dict
            Dictionary of configuration arguments.
        on_region_radius: Angle
            on region radius.
        max_wobble_distance: Angle
            maximum wobble distance (with FOV radius added).

        """
        args_dict["on_region"]["radius"] = on_region_radius
        self.define_on_region(on_region_dict=args_dict["on_region"])

        self.exclusion_mask = self.get_exclusion_mask(args_dict, max_wobble_distance)
