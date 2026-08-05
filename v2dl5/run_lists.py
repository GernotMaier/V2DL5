"""Run list selection from observation table."""

import logging
from contextlib import ExitStack

import astropy.io
import astropy.table
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time

from v2dl5 import binaries, orbital_phase

_logger = logging.getLogger(__name__)


class ZeroRunLengthError(Exception):
    """Exception raised for errors in the run length being zero."""

    def __init__(self, message="Run length is zero."):
        self.message = message
        super().__init__(self.message)


def generate_run_list(args_dict, target):
    """Read observation index table, apply selection cuts and write run list."""
    calculate_averages(args_dict)

    _logger.info("Generate run list. from %s", args_dict["obs_table"])
    obs_table = _read_observation_table(args_dict["obs_table"])
    obs_table = _apply_selection_cuts(obs_table, args_dict, target)
    try:
        _logger.info("Selected %d runs.", len(obs_table))
    except TypeError:
        _logger.warning("No runs selected (zero length run list).")
        return
    _dqm_report(obs_table, args_dict["output_dir"])
    _write_run_list(obs_table, args_dict["output_dir"])
    _print_summary(obs_table, args_dict["on_region"]["target"], target)


def calculate_averages(args_dict):
    """
    Determine mean / medium for columns relevant for RMS-dependent selection cuts.

    All other DQM cuts are applied (except target cut).

    """
    _tmp_table = _read_observation_table(args_dict["obs_table"])
    _tmp_table = _apply_selection_cuts(_tmp_table, args_dict, None)
    _logger.info(f"Runs selected to calculate averages: {len(_tmp_table)}")

    # TODO
    # add a class handling averages / outliers plus application of log / lin
    # selection cuts?


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
    try:
        if target is not None:
            obs_table = _apply_cut_target(obs_table, args_dict, target)
        obs_table = _apply_cut_mjd(obs_table, args_dict)
        obs_table = _apply_cut_atmosphere(obs_table, args_dict, target)
        obs_table = _apply_cut_dqm(obs_table, args_dict, target)
        obs_table = _apply_cut_ontime_min(obs_table, args_dict, target)
        obs_table = _apply_cut_ntel_min(obs_table, args_dict, target)
        obs_table = _apply_cut_l3rate(obs_table, args_dict, target)
    except ZeroRunLengthError as e:
        _logger.warning(e.message)
        return None
    return obs_table


def _apply_cut_l3rate(obs_table, args_dict, target):
    """Apply epoch and observation mode dependent minimum L3 Rate cut."""
    mask_V4, mask_V5, mask_V6, mask_V6_redHV = _epoch_masks(obs_table)
    mask = np.ones(len(obs_table), dtype=bool)
    for epoch_mask, epoch in zip(
        [mask_V4, mask_V5, mask_V6, mask_V6_redHV], ["V4", "V5", "V6", "V6_redHV"]
    ):
        try:
            l3rate_min = u.Quantity(args_dict["dqm"]["l3_rate_min"][epoch]).to(u.Hz)
        except KeyError:
            l3rate_min = 0.0 * u.Hz
        mask = np.array(
            mask & [((obs_table["L3RATE"] > l3rate_min.value) & epoch_mask) | ~epoch_mask]
        )
        _print_removed_runs(
            obs_table, mask[0], "L3RATE", f"{epoch} L3Rate > {l3rate_min}", target is not None
        )
    return obs_table[mask[0]]


def _apply_cut_ntel_min(obs_table, args_dict, target):
    """Apply minimum telescope cut cut."""
    try:
        ntel_min = args_dict["dqm"]["ntel_min"]
    except KeyError:
        _logger.error("KeyError: dqm.ntel_min")
        raise

    mask = np.array([row["N_TELS"] >= ntel_min for row in obs_table])
    _print_removed_runs(obs_table, mask, "N_TELS", f"NTel >= {ntel_min}", target is not None)
    obs_table = obs_table[mask]
    if target is not None:
        _logger.info(f"Minimum number of telescopes: {np.min(obs_table['N_TELS'])}")

    if len(obs_table) == 0:
        raise ZeroRunLengthError()
    return obs_table


def _apply_cut_mjd(obs_table, args_dict):
    """Apply cut on MJD (min and max)."""
    if "mjd_min" in args_dict["observations"]:
        _logger.info(f"Selecting runs after MJD {args_dict['observations']['mjd_min']}")
        mjd_min = args_dict["observations"]["mjd_min"]
        mask = np.array([Time(row["DATE-OBS"], scale="utc").mjd > mjd_min for row in obs_table])
        obs_table = obs_table[mask]
    if "mjd_max" in args_dict["observations"]:
        _logger.info(f"Selecting runs after MJD {args_dict['observations']['mjd_max']}")
        mjd_max = args_dict["observations"]["mjd_max"]
        mask = np.array([Time(row["DATE-END"], scale="utc").mjd < mjd_max for row in obs_table])
        obs_table = obs_table[mask]

    if len(obs_table) == 0:
        raise ZeroRunLengthError()
    return obs_table


def _apply_cut_ontime_min(obs_table, args_dict, target):
    """Apply ontime min cut."""
    try:
        ontime_min = u.Quantity(args_dict["dqm"]["ontime_min"]).to(u.s)
    except KeyError:
        _logger.error("KeyError: dqm.ontime_min")
        raise

    mask = np.array([row["ONTIME"] > ontime_min.value for row in obs_table])
    _print_removed_runs(obs_table, mask, "ONTIME", f"ontime < {ontime_min} m", target is not None)
    obs_table = obs_table[mask]
    if target is not None and len(obs_table["ONTIME"]) > 0:
        _logger.info(f"Minimum run time: {np.min(obs_table['ONTIME'])} s")

    if len(obs_table) == 0:
        raise ZeroRunLengthError()
    return obs_table


def _apply_cut_dqm(obs_table, args_dict, target):
    """Apply dqm cuts."""
    try:
        dqm_stat = args_dict["dqm"]["dqmstat"]
    except KeyError:
        _logger.error("KeyError: dqm.dqmstat")
        raise

    mask = np.array([row["DQMSTAT"] in dqm_stat for row in obs_table])
    _print_removed_runs(obs_table, mask, "DQMSTAT", "dqm status", target is not None)
    obs_table = obs_table[mask]
    if target is not None:
        _logger.info(f"Selected dqm status {np.unique(obs_table['DQMSTAT'])}")

    if len(obs_table) == 0:
        raise ZeroRunLengthError()
    return obs_table


def _apply_cut_atmosphere(obs_table, args_dict, target):
    """Remove fields in column "WEATHER" which are not in args_dict.atmosphere.weather."""
    try:
        weather = args_dict["atmosphere"]["weather"]
    except KeyError:
        _logger.error("KeyError: atmosphere.weather")
        raise

    try:
        mask = np.array([row["WEATHER"][0] in weather for row in obs_table])
    except IndexError:
        _logger.error("IndexError: weather")
        raise
    _print_removed_runs(obs_table, mask, "WEATHER", "weather", target is not None)
    obs_table = obs_table[mask]
    if target is not None:
        _logger.info(f"Selected weather conditions {np.unique(obs_table['WEATHER'])}")

    if len(obs_table) == 0:
        raise ZeroRunLengthError()
    return obs_table


def _apply_cut_target(obs_table, args_dict, target):
    """Apply target cut."""
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

    if len(obs_table) == 0:
        raise ZeroRunLengthError()
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
    _logger.info(f"Write run list with {len(obs_table)} runs to {output_dir}/run_list.txt")

    column_data = obs_table[np.argsort(obs_table["OBS_ID"])]["OBS_ID"]
    table_with_single_row = Table([column_data], names=["Item"])
    astropy.io.ascii.write(
        table_with_single_row,
        f"{output_dir}/run_list.txt",
        overwrite=True,
        format="no_header",
    )

    _logger.info(f"Write run table with selected runs to {output_dir}/run_list.fits.gz")

    obs_table.write(f"{output_dir}/run_list.fits.gz", overwrite=True)


def _dqm_report(obs_table, output_dir):
    """Print list of selected runs to screen."""
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
    Reject outliers from data.

    from
    stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list

    """
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d / mdev if mdev else np.zeros(len(d))
    return data[s < m], data[s > m], np.median(data), mdev


def _outlier_lists(obs_table, string, column, center_measure, outlier_type, sigma=3.0):
    """Return list of outliers for mean and median cuts."""
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
    """Print OBS_ID with more than sigma deviation from mean."""
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
    _outlier_mask_mean = np.isin(_obs_table_cleaned["OBS_ID"], _outlier_list_mean)
    _obs_table_cleaned[_outlier_mask_mean].pprint_all()
    print(f"    Outliers (median,abs): {_outlier_list_median}")
    _outlier_mask_median = np.isin(_obs_table_cleaned["OBS_ID"], _outlier_list_median)
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
    """Print min/max entry for a specific column."""
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
    """Return VERITAS Epochs as table mask."""
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


def _print_removed_runs(obs_table, mask, column_name, cut_type, print_runs=True):
    """Print removed runs."""
    if not print_runs:
        return

    if len(mask) == 0:
        _logger.info(f"No runs removed by {cut_type} cut")
        return
    _logger.info(f"Remove {len(mask)-mask.sum(~False)} runs failing {cut_type} cut")
    _removed_runs = [f"{row['OBS_ID']} ({row[column_name]})" for row in obs_table[~mask]]
    _logger.info(f"Removed following run: {_removed_runs}")
    _logger.info(f"Keep {mask.sum(~False)} runs after application of {cut_type} cut")


def _print_summary(obs_table, target_name, target):
    """Print a summary for the given run list and target."""
    _logger.info(f"Summary for target {target_name}")
    _logger.info(f"Target coordinates ra={target.ra}, dec={target.dec}")
    _logger.info(f"Selected {len(obs_table)} runs.")
    total_on_time_min = np.sum(obs_table["ONTIME"]) / 60.0
    total_on_time_hour = np.sum(obs_table["ONTIME"]) / 3600.0
    _logger.info(f"Total on time: {total_on_time_min:.2f} min" f" ({total_on_time_hour:.2f} hours)")


def split_binary_run_list(run_list_file, obs_table, binary_name, orbital_bins):
    """
    Split a run list into run list per orbital phase bins.

    Parameters
    ----------
    run_list_file : str
        Path to the run list.
    obs_table : str
        Path to observation table.
    binary_name : str
        Binary name (e.g., LS I +61 303; see v2dl5.binaries for definition).
    orbital_bins : float
        Number of bins in orbital period for averaging.

    """
    with open(run_list_file, encoding="utf-8") as f:
        run_list = [line.strip() for line in f]
    _logger.info(f"Splitting run list of length {len(run_list)} into {orbital_bins} bins")
    obs_table = astropy.table.Table.read(obs_table)
    obs_table = obs_table[np.isin(obs_table["OBS_ID"], run_list)]

    _logger.info(
        f"Using {binary_name} for orbital phase calculation (period: "
        f"{binaries.binary_properties()[binary_name]['orbital_period']} days)"
    )

    live_times = [0] * orbital_bins
    with ExitStack() as stack:
        output_files = [
            stack.enter_context(
                open(
                    f"{run_list_file.removesuffix('.txt')}_orbital_bin_{i:02d}.txt",
                    "w",
                    encoding="utf-8",
                )
            )
            for i in range(orbital_bins)
        ]

        for row in obs_table:
            phase = orbital_phase.get_orbital_phase_from_iso_time(
                iso_time=row["DATE-OBS"],
                orbital_period=binaries.binary_properties()[binary_name]["orbital_period"],
                mjd_0=binaries.binary_properties()[binary_name]["mjd_0"],
            )
            bin_index = int(phase * orbital_bins) % orbital_bins
            live_times[bin_index] += row["LIVETIME"]
            output_files[bin_index].write(f"{row['OBS_ID']}\n")

    plt.figure()
    # Create histogram with bin edges from 0 to 1
    plt.hist(
        np.linspace(0, 1, orbital_bins + 1)[:-1],
        bins=np.linspace(0, 1, orbital_bins + 1),
        weights=np.array(live_times) / 3600.0,
        histtype="step",
        linewidth=2,
    )
    plt.xlim(0, 1)
    plt.xlabel("orbital phase")
    plt.ylabel("live time (h)")
    plot_file = f"{run_list_file.removesuffix('.txt')}_live_time_per_orbital_phase_bin.png"
    plt.savefig(plot_file)
    plt.close()
    _logger.info(f"Live time plot saved to {plot_file}")

    _logger.info(f"Run lists written to {run_list_file.removesuffix('.txt')}_orbital_bin_*.txt")
    _logger.info(f"Live times per orbital phase bin: {live_times}")
