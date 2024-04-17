"""
Run list selection from observation table.

"""

import logging

import astropy.io
import astropy.table
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord

_logger = logging.getLogger(__name__)


def generate_run_list(args_dict, target):
    """
    Read observation index table, apply selection cuts and write run list.

    """

    _logger.info("Generate run list. from %s", args_dict["obs_table"])
    obs_table = _read_observation_table(args_dict["obs_table"])
    obs_table = _apply_selection_cuts(obs_table, args_dict, target)
    _logger.info("Selected %d runs.", len(obs_table))
    _dqm_report(obs_table, args_dict["output_dir"])
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
    _obs_table = _apply_cut_l3rate(_obs_table, args_dict)

    return _obs_table


def _apply_cut_l3rate(obs_table, args_dict):
    """
    Apply epoch and observation mode dependent minimum L3 Rate cut

    """

    mask_V4, mask_V5, mask_V6, mask_V6_redHV = _epoch_masks(obs_table)
    mask = np.ones(len(obs_table), dtype=bool)
    for epoch_mask, epoch in zip(
        [mask_V4, mask_V5, mask_V6, mask_V6_redHV], ["V4", "V5", "V6", "V6_redHV"]
    ):
        try:
            l3rate_min = u.Quantity(args_dict["dqm"]["l3_rate_min"][epoch]).to(u.Hz)
        except KeyError:
            l3rate_min = 0.0 * u.Hz
        mask = mask & [((obs_table["L3RATE"] > l3rate_min.value) & epoch_mask) | ~epoch_mask]
    return obs_table[mask[0]]


def _apply_cut_ntel_min(obs_table, args_dict):
    """
    Apply minimum telescope cut cut.

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
        obs_cone_radius_min = u.Quantity(args_dict["observations"]["obs_cone_radius_min"])
    except KeyError:
        obs_cone_radius_min = 0.0 * u.deg
    try:
        obs_cone_radius_max = u.Quantity(args_dict["observations"]["obs_cone_radius_max"])
    except KeyError:
        _logger.error("KeyError: observations.obs_cone_radius_mix missing")
        raise

    angular_separation = target.separation(
        SkyCoord(obs_table["RA_PNT"], obs_table["DEC_PNT"], unit=(u.deg, u.deg))
    )
    obs_table = obs_table[
        (angular_separation > obs_cone_radius_min) & (angular_separation < obs_cone_radius_max)
    ]

    _logger.info(
        f"Selecting {len(obs_table)} runs from observation cone around {target}"
        f"(min {obs_cone_radius_min}, max {obs_cone_radius_max})"
    )

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
    astropy.io.ascii.write(
        column_data,
        f"{output_dir}/run_list.txt",
        overwrite=True,
        format="no_header",
        delimiter="\n",
    )

    _logger.info(f"Write run table with selected runs to {output_dir}/run_list.fits.gz")

    obs_table.write(f"{output_dir}/run_list.fits.gz", overwrite=True)


def _dqm_report(obs_table, output_dir):
    """
    Print list of selected runs to screen

    """

    obs_table.sort("OBS_ID")
    # print general info
    obs_table[
        "OBS_ID",
        "RUNTYPE",
        "ONTIME",
        "RA_PNT",
        "DEC_PNT",
        "RA_OBJ",
        "DEC_OBJ",
    ].pprint_all()

    # print dqm
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
        print("Outliers for epoch", epoch)
        _print_min_max(obs_table, epoch_mask, "L3RATE", f"{epoch} (Hz)")
        _print_min_max(obs_table, epoch_mask, "L3RATESD", f"{epoch} (Hz)")
        _print_outlier(obs_table, epoch_mask, "L3RATE", f"{epoch} (Hz)", output_dir, "min")
        _print_outlier(obs_table, epoch_mask, "L3RATESD", f"{epoch} (Hz)", output_dir, "max", True)

    _print_min_max(obs_table, [True], "FIRMEAN1", "all (deg)")
    _print_min_max(obs_table, [True], "FIRSTD1", "all (deg)")
    _print_min_max(obs_table, [True], "FIRCORM1", "all (deg)")

    _print_outlier(obs_table, [True], "FIRMEAN1", "all (deg)", output_dir, "max")
    _print_outlier(obs_table, [True], "FIRCORM1", "all (deg)", output_dir, "max")
    _print_outlier(obs_table, [True], "FIRSTD1", "all (deg)", output_dir, "max", True)


def _reject_outliers(data, m=3.0):
    """
    from
    stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list

    """
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d / mdev if mdev else np.zeros(len(d))
    return data[s < m], data[s > m], np.median(data), mdev


def _outlier_lists(obs_table, string, column, center_measure, outlier_type, sigma=3.0):
    """
    Return list of outliers for mean and median cuts.

    """

    if center_measure == "Mean":
        m = np.mean(obs_table[column])
        s = np.std(obs_table[column])
    elif center_measure == "Median":
        _data = np.ndarray.flatten(obs_table[column])
        m = np.median(_data)
        s = np.median(np.abs(_data - m))
    print(f"{center_measure} {column} for {string}: {m:.2f} +- {s:.2f}")
    if outlier_type == "bounds":
        outlier_list = [row["OBS_ID"] for row in obs_table if abs(row[column] - m) > sigma * s]
    elif outlier_type == "max":
        outlier_list = [row["OBS_ID"] for row in obs_table if row[column] - m > sigma * s]
    elif outlier_type == "min":
        outlier_list = [row["OBS_ID"] for row in obs_table if m - row[column] > sigma * s]

    return m, s, outlier_list


def _print_outlier(
    obs_table, mask, column, string, output_dir, outlier_type="bounds", log_axis=False
):
    """
    Print OBS_ID with more than sigma deviation from mean

    """
    _obs_table_cleaned = obs_table[(obs_table[column] > -9998.0) & mask]
    if log_axis:
        _obs_table_cleaned[column] = np.log10(_obs_table_cleaned[column])

    _mean, _std, _outlier_list_mean = _outlier_lists(
        _obs_table_cleaned, string, column, "Mean", outlier_type, sigma=2.0
    )
    _median, _abs_deviation, _outlier_list_median = _outlier_lists(
        _obs_table_cleaned, string, column, "Median", outlier_type, sigma=3.0
    )

    print(f"{column} for {string}:")
    print(f"    Outliers (mean,std): {_outlier_list_mean}")
    _outlier_mask_mean = np.in1d(_obs_table_cleaned["OBS_ID"], _outlier_list_mean)
    _obs_table_cleaned[_outlier_mask_mean].pprint_all()
    print(f"    Outliers (median,abs): {_outlier_list_median}")
    _outlier_mask_median = np.in1d(_obs_table_cleaned["OBS_ID"], _outlier_list_median)
    _obs_table_cleaned[_outlier_mask_median].pprint_all()

    _plot_outliers(
        _obs_table_cleaned,
        column,
        string,
        _outlier_mask_mean,
        _outlier_mask_median,
        _mean,
        _std,
        _median,
        _abs_deviation,
        output_dir,
    )


def _plot_outliers(
    obs_table, column, string, mask_std, mask_med, mean, std, median, abs_deviation, output_dir
):
    """
    Plot distribution of column values include outliers.
    Indicate mean and std as background color (gray)

    """

    plt.axvspan(mean - std, mean + std, color="red", alpha=0.15)
    plt.axvline(mean, color="red", linestyle="--", linewidth=2)
    plt.axvspan(median - abs_deviation, median + abs_deviation, color="green", alpha=0.15)
    plt.axvline(median, color="green", linestyle="--", linewidth=2)
    hist_bins = plt.hist(obs_table[column], bins=50)[1]
    plt.hist(obs_table[column][mask_med], bins=hist_bins, color="green")
    plt.hist(obs_table[column][mask_std], bins=hist_bins, color="red")
    plt.title(f"{column} for {string}")
    plt.xlabel(column)
    plt.ylabel("Number of runs")
    plt.savefig(f"{output_dir}/outliers_{column}_{string}.png")
    plt.close()


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

    return (
        np.array(mask_V4),
        np.array(mask_V5),
        np.array(mask_V6),
        np.array(mask_V6_redHV),
    )
