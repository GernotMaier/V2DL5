#!/usr/bin/python

"""
Plotting of binary light curves.

Includes the following plots:

    - Light curve vs time
    - Light curve vs orbital phase
    - Light curve vs orbital phase for each orbit
"""

import argparse
import logging

import v2dl5.light_curves.binary_plotting
import v2dl5.light_curves.data_reader
from v2dl5 import binaries


def _parse():
    """
    Parse command line arguments.

    Returns
    -------
    dict
        Command line parameters.

    """
    parser = argparse.ArgumentParser(description="Binary light-curve plotting.")

    parser.add_argument(
        "--instrument",
        type=str,
        required=True,
        help="Source of data (VERITAS, HESS, Gamma, Swift XRT).",
    )
    parser.add_argument(
        "--configuration",
        type=str,
        required=True,
        help="Configuration describing data files and plotting",
    )
    parser.add_argument(
        "--binary_name",
        type=str,
        required=True,
        default="LS I +61 303",
        help="Binary name (e.g., LS I +61 303; see v2dl5.binaries for definition).",
    )
    parser.add_argument(
        "--orbital_bins",
        type=int,
        required=False,
        default=16,
        help="Number of orbital bins",
    )
    parser.add_argument(
        "--plot_type",
        type=str,
        required=False,
        default=".pdf",
        help="File type for plots (e.g., '.pdf', '.png')",
    )
    parser.add_argument(
        "--figure_dir",
        type=str,
        required=False,
        default="./figures/",
        help="Directory for saving figures",
    )

    return parser.parse_args()


def main():
    """Binary light-curve plotting."""
    logging.root.setLevel(logging.INFO)

    args = _parse()
    if args.orbital_bins < 1:
        if args.instrument in ["VERITAS", "Gamma", "HESS"]:
            args.orbital_bins = 10
        else:
            args.orbital_bins = 20

    logging.info("Light Curve Analysis - run parameters")
    logging.info(f"instrument: {args.instrument}")
    logging.info(f"instrument list: {args.configuration}")
    logging.info(f"number of bins for averaging: {args.orbital_bins}")

    try:
        binary = binaries.binary_properties()[args.binary_name]
    except KeyError:
        raise KeyError(f"Binary {args.binary_name} not found in binaries.py")

    data_reader = v2dl5.light_curves.data_reader.LightCurveDataReader(
        args.configuration, binary=binary
    )
    data_reader.read_data()

    plotter = v2dl5.light_curves.binary_plotting.BinaryLightCurvePlotter(
        data=data_reader.data_dict, config=data_reader.config, binary=binary
    )
    for y_key in ["flux", "index", "live_time"]:
        plotter.plot_flux_vs_time(
            "MJD", y_key, None, None, file_type=args.plot_type, figure_dir=args.figure_dir
        )
        plotter.plot_flux_vs_time(
            "orbital phase", y_key, None, None, file_type=args.plot_type, figure_dir=args.figure_dir
        )
        plotter.plot_flux_vs_phase_for_individual_orbits(
            instrument=args.instrument,
            y_axis=y_key,
            file_type=args.plot_type,
            figure_dir=args.figure_dir,
        )
        plotter.plot_flux_vs_orbit_number(
            instrument=args.instrument,
            y_axis=y_key,
            phase_bins=args.orbital_bins,
            file_type=args.plot_type,
            figure_dir=args.figure_dir,
        )
        plotter.plot_distribution(
            instrument=args.instrument,
            y_axis=y_key,
            file_type=args.plot_type,
            figure_dir=args.figure_dir,
        )

    plotter.plot_index_vs_flux(
        instrument=args.instrument,
        file_type=args.plot_type,
        figure_dir=args.figure_dir,
    )

    plotter.plot_live_time_vs_phase_bin(
        instrument=args.instrument,
        phase_bins=args.orbital_bins,
        file_type=args.plot_type,
        figure_dir=args.figure_dir,
    )


if __name__ == "__main__":
    main()
