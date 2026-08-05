#!/usr/bin/python3
"""
Add spectral fit parameters from YAML files to light curve table.

Example:
    python add_spectral_fits_to_light_curve_table.py \
        --light-curve-file nightly.ecsv \
        --nightly-spectra-directory spectral_results/nightly \
        --output-file output.ecsv
"""

import argparse
import logging
from pathlib import Path

import numpy as np
import yaml
from astropy.table import Table

_logger = logging.getLogger(__name__)


def _parse():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Add spectral fit parameters to light curve table")
    parser.add_argument(
        "--light_curve_file",
        type=str,
        required=True,
        help="Input light curve ECSV file",
    )
    parser.add_argument(
        "--nightly_spectra_directory",
        type=str,
        required=True,
        help="Directory containing nightly spectral YAML files",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Output ECSV file",
    )
    return parser.parse_args()


def get_mjd_range(filename):
    """Extract MJD range from filename.

    Parameters
    ----------
    filename : str
        Filename in format *_MJD_<start>_<end>.yaml

    Returns
    -------
    tuple
        (mjd_min, mjd_max)
    """
    parts = Path(filename).stem.split("_")
    return int(parts[-2]), int(parts[-1])


def add_spectral_columns(table):
    """Add new columns for spectral parameters.

    Parameters
    ----------
    table : astropy.table.Table
        Input table

    Returns
    -------
    astropy.table.Table
        Table with new columns
    """
    n = len(table)
    new_columns = {
        "flux_norm": np.full(n, np.nan),
        "flux_norm_err": np.full(n, np.nan),
        "index": np.full(n, np.nan),
        "index_err": np.full(n, np.nan),
        "chi2": np.full(n, np.nan),
        "ndf": np.full(n, np.nan),
        "e0": np.full(n, np.nan),
    }
    for name, values in new_columns.items():
        table[name] = values
    return table


def process_yaml_file(yaml_path, table):
    """Process single YAML file and update table.

    Parameters
    ----------
    yaml_path : pathlib.Path
        Path to YAML file
    table : astropy.table.Table
        Table to update
    """
    mjd_min, mjd_max = get_mjd_range(yaml_path)
    _logger.info(f"Processing {yaml_path} for MJD range {mjd_min}-{mjd_max}")

    with open(yaml_path) as f:
        data = yaml.safe_load(f)

    mask = (table["time_min"] == mjd_min) & (table["time_max"] == mjd_max)
    idx = np.where(mask)[0]
    if len(idx) == 0:
        _logger.warning(f"No matching entry in table for MJD range {mjd_min}-{mjd_max}")
        return
    i = idx[0]

    fit = data["fit_parameters"][0]
    table["chi2"][i] = fit.get("chi2", np.nan)
    table["ndf"][i] = fit.get("ndf", np.nan)
    table["e0"][i] = fit.get("e0_TeV", np.nan)

    for p in fit.get("parameters", []):
        if p["name"] == "p0":
            table["flux_norm"][i] = p.get("value", np.nan)
            table["flux_norm_err"][i] = p.get("error", np.nan)
        elif p["name"] == "p1":
            table["index"][i] = p.get("value", np.nan)
            table["index_err"][i] = p.get("error", np.nan)


def main():
    """Add spectral fit parameters to light curve table."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    args = _parse()

    _logger.info(f"Reading light curve table from {args.light_curve_file}")
    table = Table.read(args.light_curve_file, format="ascii.ecsv")
    table = add_spectral_columns(table)

    yaml_dir = Path(args.nightly_spectra_directory)
    _logger.info(f"Processing YAML files from {yaml_dir}")
    for yaml_path in yaml_dir.glob("runs_MJD_*.yaml"):
        process_yaml_file(yaml_path, table)

    table.write(args.output_file, format="ascii.ecsv", overwrite=True)
    _logger.info(f"Written updated table to {args.output_file}")


if __name__ == "__main__":
    main()
