#!/usr/bin/python

"""
Generate list of observation nights from anasum ROOT file.

Read an anasum ROOT file and extract the list of unique observation nights
using the start and stop MJD times from the run summary tree.

Example:
    python get_observation_nights.py --anasum_file analysis.anasum.root --output_file nights.txt
"""

import argparse
import logging
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


def _parse():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract list of observation nights from anasum ROOT file."
    )
    parser.add_argument(
        "--anasum_file",
        type=str,
        required=True,
        help="Path to anasum ROOT file",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Output file for list of nights",
    )
    return parser.parse_args()


def get_unique_nights(start_mjd, stop_mjd):
    """
    Get list of unique observation nights.

    Parameters
    ----------
    start_mjd : numpy.ndarray
        Array of run start times in MJD
    stop_mjd : numpy.ndarray
        Array of run stop times in MJD

    Returns
    -------
    numpy.ndarray
        Array of unique observation nights (integer MJD)
    """
    # Get integer part of MJD for start and stop times
    return np.unique(np.concatenate([
        np.floor(ak.to_numpy(start_mjd)).astype(int),
        np.floor(ak.to_numpy(stop_mjd)).astype(int)
    ]))


def main():
    """Extract list of observation nights."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    logger = logging.getLogger(__name__)

    args = _parse()

    logger.info(f"Reading anasum file: {args.anasum_file}")
    tree = uproot.open(args.anasum_file)["total_1/stereo/tRunSummary"]

    run_on = tree["runOn"].array()
    mjd_start = tree["MJDrunstart"].array()
    mjd_stop = tree["MJDrunstop"].array()

    mask = run_on >= 0
    mjd_start = mjd_start[mask]
    mjd_stop = mjd_stop[mask]

    logger.info(f"Found {len(run_on)} total runs")
    logger.info(f"Found {ak.sum(mask)} valid runs (runOn >= 0)")

    nights = get_unique_nights(mjd_start, mjd_stop)
    logger.info(f"Found {len(nights)} unique observation nights")

    night_bins = np.column_stack((nights, nights + 1))

    output_path = Path(args.output_file)
    np.savetxt(output_path, night_bins, fmt='%d %d')
    logger.info(f"Written night list to {output_path}")

    logger.info(f"First night: MJD {nights[0]}-{nights[0]+1}")
    logger.info(f"Last night:  MJD {nights[-1]}-{nights[-1]+1}")
    logger.info(f"Time span:   {nights[-1] - nights[0]} days")


if __name__ == "__main__":
    main()
