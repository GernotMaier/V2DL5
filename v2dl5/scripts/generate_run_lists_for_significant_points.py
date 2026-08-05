#!/usr/bin/python

"""
Generate run lists for time periods with significant detections.

Read a light curve ECSV file and an anasum ROOT file, identify time windows
with significant detections, and generate run lists for those time windows. The run lists
are saved in a specified output directory. Used as input for e.g., time-dependent spectral
analysis.

"""
import argparse
import logging
from pathlib import Path

import numpy as np
import uproot
from astropy.table import Table

logger = logging.getLogger(__name__)


def _parse():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate run lists for time periods with significant detections."
    )

    parser.add_argument(
        "--min_significance",
        type=float,
        default=5.0,
        help="Minimum significance for light curve points to be considered",
    )
    parser.add_argument(
        "--light_curve_file",
        type=str,
        required=True,
        help="Path to light curve ECSV file",
    )
    parser.add_argument(
        "--anasum_file",
        type=str,
        required=True,
        help="Path to anasum ROOT file",
    )
    parser.add_argument(
        "--output_directory",
        type=str,
        default="./run_lists",
        help="Output directory for run lists (default: ./run_lists)",
    )

    return parser.parse_args()


def main():
    """Generate run lists for significant light-curve points."""
    logging.root.setLevel(logging.INFO)
    args = _parse()

    table = Table.read(args.light_curve_file, format="ascii.ecsv")
    time_windows = [
        (row["time_min"], row["time_max"])
        for row in table
        if row["significance"] > args.min_significance
    ]
    logger.info("Found %d significant (>%s) time windows.", len(time_windows), args.min_significance)

    tree = uproot.open(args.anasum_file)["total_1/stereo/tRunSummary"]
    mjd_on = tree["MJDOn"].array()
    run_on = tree["runOn"].array()

    output_dir = Path(args.output_directory)
    output_dir.mkdir(parents=True, exist_ok=True)

    for tmin, tmax in time_windows:
        mask = (mjd_on >= tmin) & (mjd_on <= tmax)
        matched = run_on[mask]

        output_file = output_dir / f"runs_MJD_{int(tmin)}_{int(tmax)}.txt"
        np.savetxt(output_file, matched, fmt="%d")
        print(f"Runs for MJD {tmin}-{tmax} with {len(matched)} runs written to {output_file}")


if __name__ == "__main__":
    main()
