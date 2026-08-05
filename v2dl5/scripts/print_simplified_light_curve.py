#!/usr/bin/python
"""Read gammapy-generated light curve files and print a simplified/readable version of the data."""

import argparse
import logging

import numpy as np
from astropy.table import Table


def extract_single_value(column):
    """Extract a single value from a column."""
    return [
        float(value[0]) if isinstance(value, (list, np.ndarray)) else float(value)
        for value in column
    ]


def main():
    """Read and print light curve data."""
    parser = argparse.ArgumentParser(description="Read and print light curve data.")

    parser.add_argument(
        "--file",
        type=str,
        required=True,
        help="File (ecsv format) containing light curve data.",
    )

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    complicated_table = Table.read(args.file, format="ascii.ecsv")
    data_short = {
        "e_min": extract_single_value(complicated_table["e_min"]),
        "time_min": complicated_table["time_min"],
        "time_max": complicated_table["time_max"],
        "flux": extract_single_value(complicated_table["flux"]),
        "flux_err": extract_single_value(complicated_table["flux_err"]),
        "flux_ul": extract_single_value(complicated_table["flux_ul"]),
        "significance": extract_single_value(complicated_table["sqrt_ts"]),
    }

    short_table = Table(data_short)
    print(short_table)


if __name__ == "__main__":
    main()
