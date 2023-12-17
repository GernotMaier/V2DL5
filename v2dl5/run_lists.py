"""
Run list selection from observation table.

"""

import logging

import astropy.table
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import ascii

_logger = logging.getLogger(__name__)


def generate_run_list(args_dict, target):
    """
    Read observation index table, apply selection cuts and write run list.

    """

    _logger.info("Generate run list. from %s", args_dict["obs_table"])

    obs_table = _read_observation_table(args_dict["obs_table"])

    obs_table = _apply_selection_cuts(obs_table, args_dict, target)

    _logger.info("Selected %d runs.", len(obs_table))

    _dqm_report(obs_table)

    _write_run_list(obs_table, args_dict["output_dir"])


def _read_observation_table(obs_table_file_name):
    """
    Read observation table from obs_index file.

    Fill masked values for the following fields with default values:

    - DQMSTAT: "unknown"

    """

    obs_table = astropy.table.Table.read(obs_table_file_name)
    obs_table["DQMSTAT"].fill_value = "unknown"

    return obs_table.filled()


def _apply_selection_cuts(obs_table, args_dict, target):
    """
    Apply selection cuts to observation table.

    Parameters
    ----------
    obs_table : `~astropy.table.Table`
        Observation table.

    """

    _logger.info("Apply selection cuts.")

    _obs_table = _apply_cut_target(obs_table, args_dict, target)
    _obs_table = _apply_cut_atmosphere(_obs_table, args_dict)
    _obs_table = _apply_cut_dqm(_obs_table, args_dict)
    _obs_table = _apply_cut_ontime_min(_obs_table, args_dict)
    _obs_table = _apply_cut_ntel_min(_obs_table, args_dict)

    return _obs_table


def _apply_cut_ntel_min(obs_table, args_dict):
    """
    Apply mininimum telescope cut cut.

    """

    try:
        ntel_min = args_dict["dqm"]["ntel_min"]
    except KeyError:
        _logger.error("KeyError: dqm.ntel_min")
        raise

    mask = [row["N_TELS"] >= ntel_min for row in obs_table]
    _logger.info(f"Remove {mask.count(False)} runs with ntel < {ntel_min}")
    obs_table = obs_table[mask]
    _logger.info(f"Minimum number of telescopes: {np.min(obs_table['N_TELS'])}")

    return obs_table


def _apply_cut_ontime_min(obs_table, args_dict):
    """
    Apply ontime min cut.

    """

    try:
        ontime_min = u.Quantity(args_dict["dqm"]["ontime_min"]).to(u.s)
    except KeyError:
        _logger.error("KeyError: dqm.ontime_min")
        raise

    mask = [row["ONTIME"] > ontime_min.value for row in obs_table]
    _logger.info(f"Remove {mask.count(False)} runs with ontime < {ontime_min}")
    obs_table = obs_table[mask]
    _logger.info(f"Minimum run time: {np.min(obs_table['ONTIME'])} s")

    return obs_table


def _apply_cut_dqm(obs_table, args_dict):
    """
    Apply dqm cuts

    """

    try:
        dqm_stat = args_dict["dqm"]["dqmstat"]
    except KeyError:
        _logger.error("KeyError: dqm.dqmstat")
        raise

    mask = [row["DQMSTAT"] in dqm_stat for row in obs_table]
    _logger.info(f"Remove {mask.count(False)} runs with dqm status not not in {dqm_stat}")
    obs_table = obs_table[mask]
    _logger.info(f"Selected dqm status {np.unique(obs_table['DQMSTAT'])}")

    return obs_table


def _apply_cut_atmosphere(obs_table, args_dict):
    """
    Remove all fields in column "WEATHER" which are not in the list of args_dict.atmosphere.weather

    """

    try:
        weather = args_dict["atmosphere"]["weather"]
    except KeyError:
        _logger.error("KeyError: atmosphere.weather")
        raise

    try:
        mask = [row["WEATHER"][0] in weather for row in obs_table]
    except IndexError:
        _logger.error("IndexError: weather")
        raise
    _logger.info(f"Remove {mask.count(False)} runs with weather not in {weather}")
    obs_table = obs_table[mask]
    _logger.info(f"Selected weather conditions {np.unique(obs_table['WEATHER'])}")

    return obs_table


def _apply_cut_target(obs_table, args_dict, target):
    """
    Apply target cut.

    """

    try:
        obs_cone_radius = u.Quantity(args_dict["observations"]["obs_cone_radius"])
    except KeyError:
        _logger.error("KeyError: observations.obs_cone_radius")
        raise

    angular_separation = target.separation(
        SkyCoord(obs_table["RA_PNT"], obs_table["DEC_PNT"], unit=(u.deg, u.deg))
    )
    obs_table = obs_table[angular_separation < obs_cone_radius]

    _logger.info("Selecting %d runs from observation cone around %s", len(obs_table), target)

    return obs_table


def _write_run_list(obs_table, output_dir):
    """
    Write run list.

    Parameters
    ----------
    obs_table : `~astropy.table.Table`
        Observation table.
    output_dir : str
        Output directory.

    """

    _logger.info(f"Write run list to {output_dir}/run_list.txt")

    column_data = obs_table[np.argsort(obs_table["OBS_ID"])]["OBS_ID"]
    ascii.write(
        column_data,
        f"{output_dir}/run_list.txt",
        overwrite=True,
        format="no_header",
        delimiter="\n",
    )

    _logger.info(f"Write run table with selected runs to {output_dir}/run_list.fits.gz")

    obs_table.write(f"{output_dir}/run_list.fits.gz", overwrite=True)


def _dqm_report(obs_table):
    """
    Print list of selected runs to screen

    """

    obs_table.sort("OBS_ID")

    obs_table[
        "OBS_ID",
        "RUNTYPE",
        "DATACAT",
        "N_TELS",
        "TELLIST",
        "DQMSTAT",
        "WEATHER",
        "L3RATE",
        "L3RATESD",
        "FIRMEAN1",
        "FIRCORM1",
        "FIRSTD1",
    ].pprint_all()

    mask_V4, mask_V5, mask_V6, mask_V6_redHV = _epoch_masks(obs_table)

    for epoch_mask, epoch in zip(
        [mask_V4, mask_V5, mask_V6, mask_V6_redHV], ["V4", "V5", "V6", "V6_redHV"]
    ):
        _print_min_max(obs_table, epoch_mask, "L3RATE", f"{epoch} (Hz)")
        _print_min_max(obs_table, epoch_mask, "L3RATESD", f"{epoch} (Hz)")
        _print_min_max(obs_table, epoch_mask, "FIRMEAN1", f"{epoch} (deg)")
        _print_min_max(obs_table, epoch_mask, "FIRSTD1", f"{epoch} (deg)")
        _print_min_max(obs_table, epoch_mask, "FIRCORM1", f"{epoch} (deg)")

        _print_outlier(obs_table, epoch_mask, "L3RATE", f"{epoch} (Hz)")
        _print_outlier(obs_table, epoch_mask, "FIRMEAN1", f"{epoch} (deg)")
        _print_outlier(obs_table, epoch_mask, "FIRCORM1", f"{epoch} (deg)")
        _print_outlier(obs_table, epoch_mask, "FIRSTD1", f"{epoch} (deg)")


def _reject_outliers(data, m=3.0):
    """
    from
    stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list

    """
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d / mdev if mdev else np.zeros(len(d))
    print("AAAA", mdev, data[s > m])
    return data[s < m]


def _print_outlier(obs_table, mask, column, string, sigma=2):
    """
    Print OBS_ID with more than sigma deviation from mean

    """

    _obs_table_cleaned = obs_table[(obs_table[column] > -9998.0) & mask]

    _mean = np.mean(_obs_table_cleaned[column])
    _std = np.std(_obs_table_cleaned[column])
    _ff = np.ndarray.flatten(_obs_table_cleaned[column])
    _rr = _reject_outliers(_ff)
    print("BBBB", len(_ff), len(_rr))
    print(f"Mean {column} for {string}: {_mean:.2f} +- {_std:.2f}")

    # get list of obs_ids with more than sigma deviation from mean
    _outlier_list = [
        row["OBS_ID"] for row in _obs_table_cleaned if abs(row[column] - _mean) > sigma * _std
    ]
    print(f"{column} for {string}:")
    print(f"    Outliers: {_outlier_list}")
    for _outlier in _outlier_list:
        _tmp_obs = _obs_table_cleaned[_obs_table_cleaned["OBS_ID"] == _outlier]
        _tmp_obs.pprint_all()


def _print_min_max(obs_table, mask, column, string):
    """
    Print min/max entry for a specific column

    """

    _obs_table_cleaned = obs_table[(obs_table[column] > -9998.0) & mask]
    try:
        min_index = np.argmin(_obs_table_cleaned[column])
        max_index = np.argmax(_obs_table_cleaned[column])

        print(f"{column} for {string}:")
        print(
            f"    Max for obs_id {_obs_table_cleaned['OBS_ID'][max_index]}: "
            f"{_obs_table_cleaned[column][max_index]:.2f}"
        )
        print(
            f"    Min for obs_id {_obs_table_cleaned['OBS_ID'][min_index]}: "
            f"{_obs_table_cleaned[column][min_index]:.2f}"
        )
    except ValueError:
        _logger.warning(f"Empty list for min/max determination of {column}")


def _epoch_masks(obs_table):
    """
    Return VERITAS Epochs as table mask

    """

    mask_V4 = [row["OBS_ID"] < 46642 for row in obs_table]
    mask_V5 = [row["OBS_ID"] < 63372 and row["OBS_ID"] > 46642 for row in obs_table]
    mask_V6 = [row["OBS_ID"] > 63372 and row["RUNTYPE"] == "observing" for row in obs_table]
    mask_V6_redHV = [row["OBS_ID"] > 63372 and row["RUNTYPE"] == "obsLowHV" for row in obs_table]

    return mask_V4, mask_V5, mask_V6, mask_V6_redHV
