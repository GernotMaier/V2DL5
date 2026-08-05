"""Binary light curve plotting."""

import logging
import math

import matplotlib.pyplot as plt
import numpy as np

import v2dl5.plotting.utilities as plotting_utilities


class BinaryLightCurvePlotter:
    """
    Binary light curve plotter.

    Parameters
    ----------
    data : dict
        Data dictionary with light-curve data.
    config : dict
        Configuration dictionary.
    binary : dict
        Binary parameters.

    """

    def __init__(self, data, config, binary):
        self._logger = logging.getLogger(__name__)

        self.data = data
        self.config = config
        self.binary = binary

    def plot_flux_vs_time(
        self,
        time_axis="MJD",
        y_axis="flux",
        mjd_min=None,
        mjd_max=None,
        orbit_number=None,
        phase_min=None,
        phase_max=None,
        file_type=".pdf",
        figure_dir="./figures/",
        axes=None,
        fontsize=None,
        markersize=None,
    ):
        """
        Plot data (e.g., 'flux') vs time (MJD or orbital phase).

        Parameters
        ----------
        time_axis: str
            Time axis (MJD or orbital phase).
        y_axis: str
            Y-axis data (e.g. 'flux').
        mjd_min: float
            Minimum MJD value.
        mjd_max: float
            Maximum MJD value.
        phase_min: float
            Minimum phase value.
        phase_max: float
            Maximum phase value.
        orbit_number: int
            Orbit number to plot.
        file_type: str
            File type for plots (e.g., '.pdf', '.png').
        figure_dir: str
            Directory for saving figures.
        axes: matplotlib.axes.Axes
            Axes object to plot on.
        fontsize: int
            Font size for labels.
        markersize: int
            Marker size for points.
        """
        self._logger.info(
            f"Plotting {y_axis} vs {time_axis}(MJD {mjd_min}, {mjd_max}, orbit id {orbit_number})"
        )

        ax = axes if axes else plotting_utilities.paper_figures(None, None)

        for idx, (data_instrument, data) in enumerate(self.data.items()):
            if not self.plot_this_instrument(self.config[idx], y_axis):
                continue
            color, marker = self.get_marker_and_color(idx)

            x, y, e, x_ul, y_ul = self._get_light_curve_in_mjd_limits(
                data,
                y_axis,
                time_axis,
                mjd_min,
                mjd_max,
                orbit_number,
                phase_min,
                phase_max,
                self.config[idx].get("significance_min", None),
            )
            plt.errorbar(
                x,
                y,
                e,
                None,
                label=(
                    data_instrument
                    if self.config[idx].get("plot_label") is None
                    else self.config[idx]["plot_label"]
                ),
                color=color,
                marker=marker,
                linestyle="none",
                fillstyle="full",
                linewidth=plotting_utilities.get_line_width(),
                markersize=(
                    plotting_utilities.get_marker_size() if markersize is None else markersize
                ),
            )
            if len(y_ul) > 0:
                plt.errorbar(
                    x_ul,
                    y_ul,
                    yerr=0.1 * max(y_ul),
                    color=color,
                    fmt="_",
                    linestyle="none",
                    fillstyle="none",
                    uplims=True,
                    linewidth=plotting_utilities.get_line_width(),
                    markersize=plotting_utilities.get_marker_size() * 0.5,
                )

        ax.set_ylim([self.config[0].get("flux_axis_min"), self.config[0].get("flux_axis_max")])
        ax.axhline(0, color="lightgray", linestyle="--")
        plt.xlabel(self._get_time_axis_label(time_axis), fontsize=fontsize)
        if mjd_min is not None and mjd_max is not None:
            ax.set_xlim([mjd_min, mjd_max])
        plt.ylabel(self.config[0].get(y_axis + "_axis_label", ""), fontsize=fontsize)
        if axes is None:
            plt.legend()
            plotting_utilities.print_figure(
                f"Light-Curve-{self.binary['name']}-{y_axis}-vs-{time_axis.replace(' ', '-')}",
                file_type=file_type,
                figure_dir=figure_dir,
            )

    def _get_number_columns_and_rows(self, number):
        """Get number of columns and rows for plotting."""
        n_columns = 4
        fig_size = (10, 4)
        if number == 9:
            n_columns = 3
            fig_size = (10, 8)
        if number > 9:
            n_columns = 6
            fig_size = (12, 6)
        if number == 16:
            n_columns = 4
            fig_size = (10, 8)
        if number > 16:
            n_columns = 5
            fig_size = (16, 10)
        if number > 21:
            n_columns = 8
            fig_size = (16, 10)
        return n_columns, math.ceil(number / n_columns), fig_size

    def plot_flux_vs_phase_for_individual_orbits(
        self, instrument, y_axis="flux", file_type=".pdf", figure_dir="./figures/"
    ):
        """Plot flux vs phase with one plot per orbit."""
        self._logger.info(f"Plotting {y_axis} vs phase for individual orbits")
        time_axis = "orbital phase"

        data = self.data[instrument]
        if y_axis not in data:
            self._logger.warning(f"Y-axis {y_axis} not found in data for {instrument}")
            return
        orbits = sorted(set(data["orbit_number"]))
        self._logger.info(f"Orbits for {instrument} (total number {len(orbits)}): {orbits}")

        y_min, y_max = self._global_value_extrema(data[y_axis], data.get(y_axis + "_err", 0.0))
        fontsize = 6
        n_columns, n_rows, figsize = self._get_number_columns_and_rows(len(orbits))
        self._logger.info(f"Number of columns: {n_columns}, number of rows: {n_rows}")
        plt.figure(figsize=figsize)

        for i, orbit_id in enumerate(orbits):
            axes = plt.subplot(n_rows, n_columns, i + 1)
            axes.set_xlim([0.0, 1.0])
            axes.set_ylim([y_min, y_max])
            self.plot_flux_vs_time(
                time_axis=time_axis,
                y_axis=y_axis,
                mjd_min=None,
                mjd_max=None,
                orbit_number=orbit_id,
                axes=axes,
                fontsize=fontsize,
                markersize=plotting_utilities.get_marker_size() * 0.5,
            )
            plt.rc("xtick", labelsize=fontsize)
            plt.rc("ytick", labelsize=fontsize)
            plt.text(
                0.05,
                0.95,
                f"Orbit {orbit_id}",
                transform=plt.gca().transAxes,
                fontsize=fontsize,
                verticalalignment="top",
            )

        plt.tight_layout()
        plotting_utilities.print_figure(
            f"Light-Curve-Orbits-{self.binary['name']}-{y_axis}-vs-{time_axis.replace(' ', '-')}",
            file_type=file_type,
            figure_dir=figure_dir,
        )

    def plot_flux_vs_orbit_number(
        self, instrument, y_axis="flux", phase_bins=10, file_type=".pdf", figure_dir="./figures/"
    ):
        """Plot flux vs orbit number."""
        self._logger.info(f"Plotting {y_axis} vs orbit number")
        time_axis = "orbit number"

        data = self.data[instrument]
        if y_axis not in data:
            self._logger.warning(f"Y-axis {y_axis} not found in data for {instrument}")
            return
        orbits = sorted(set(data["orbit_number"]))
        self._logger.info(f"Orbits for {instrument} (total number {len(orbits)}): {orbits}")

        y_min, y_max = self._global_value_extrema(data[y_axis], data.get(y_axis + "_err", 0.0))
        fontsize = 6
        n_columns, n_rows, figsize = self._get_number_columns_and_rows(phase_bins)
        self._logger.info(f"Number of columns: {n_columns}, number of rows: {n_rows}")
        plt.figure(figsize=figsize)

        for i in range(phase_bins):
            axes = plt.subplot(n_rows, n_columns, i + 1)
            axes.set_xlim([min(orbits)-1., max(orbits)+1.])
            axes.set_ylim([y_min, y_max])

            phase_min = i * (1.0 / phase_bins)
            phase_max = (i + 1) * (1.0 / phase_bins)

            self.plot_flux_vs_time(
                time_axis=time_axis,
                y_axis=y_axis,
                mjd_min=None,
                mjd_max=None,
                orbit_number=None,
                phase_min=phase_min,
                phase_max=phase_max,
                axes=axes,
                fontsize=fontsize,
                markersize=plotting_utilities.get_marker_size() * 0.5,
            )
            plt.rc("xtick", labelsize=fontsize)
            plt.rc("ytick", labelsize=fontsize)
            plt.text(
                0.05,
                0.95,
                f"Orbit phase {phase_min:.2f} - {phase_max:.2f}",
                transform=plt.gca().transAxes,
                fontsize=fontsize,
                verticalalignment="top",
            )

        plt.tight_layout()
        plotting_utilities.print_figure(
            (
                f"Light-Curve-Orbit-Number-{self.binary['name']}-"
                f"{y_axis}-vs-{time_axis.replace(' ', '-')}"
            ),
            file_type=file_type,
            figure_dir=figure_dir,
        )

    def _global_value_extrema(self, y, y_err=0.0):
        """Absolute min and max of y values."""
        y = np.array(y)
        y_err = np.array(y_err)

        y = y[~np.isnan(y) & np.isfinite(y)]
        y_err = y_err[~np.isnan(y_err) & np.isfinite(y_err)]

        if len(y) == 0:
            return -1, 1

        return (
            min([min(y) - max(y_err), 0.0]),
            max(y + y_err) * 1.2,
        )

    def _get_time_axis_label(self, time_axis):
        """Return time axis label."""
        if time_axis == "orbital phase":
            return "phase"
        if time_axis == "orbit number":
            return "orbit number"
        return "Modified Julian Date (MJD)"

    def _get_light_curve_in_mjd_limits(
        self,
        data,
        y_key,
        time_axis,
        mjd_min=None,
        mjd_max=None,
        orbit_number=None,
        phase_min=None,
        phase_max=None,
        min_significance=None,
    ):
        """
        Get light curve restricted in MJD or for a certain orbit number.

        Parameters
        ----------
        data : dict
            Light-curve data.
        y_key : str
            Key for y-axis data (e.g. 'flux').
        time_axis : str
            Time axis (MJD or orbital phase)
        mjd_min : float
            Minimum MJD value.
        mjd_max : float
            Maximum MJD value.
        orbit_number : int
            Select orbit number.
        phase_min : float
            Minimum phase value.
        phase_max : float
            Maximum phase value.
        min_significance : float
            Minimum significance for values.

        Returns
        -------
        list
            Light-curve data as lists.
        """
        x, y, e, x_ul, y_ul = [], [], [], [], []
        if y_key not in data:
            self._logger.warning(f"Y-axis {y_key} not found in data")
            return x, y, e, x_ul, y_ul
        mjd = data["MJD"]
        for i, t in enumerate(mjd):
            if (
                min_significance is not None
                and data.get("significance", [None] * len(mjd))[i] < min_significance
            ):
                continue
            if (mjd_min is not None and t < mjd_min) or (mjd_max is not None and t > mjd_max):
                continue
            if orbit_number is not None and data["orbit_number"][i] != orbit_number:
                continue
            if (phase_min is not None and data["phase"][i] < phase_min) or (
                phase_max is not None and data["phase"][i] > phase_max
            ):
                continue
            if time_axis == "MJD":
                w_x = data["MJD"][i]
            elif time_axis == "orbital phase":
                w_x = data["phase"][i]
            elif time_axis == "orbit number":
                w_x = data["orbit_number"][i]
            else:
                raise ValueError(f"Unknown time axis: {time_axis}")
            y_key_ul = y_key + "_ul"
            y_key_err = y_key + "_err"
            _ul = data.get(y_key_ul, [None] * len(mjd))[i]
            if _ul is not None and _ul > 0:
                x_ul.append(w_x)
                y_ul.append(data[y_key_ul][i])
            elif y_key_ul not in data or data[y_key_ul][i] < 0:
                x.append(w_x)
                y.append(data[y_key][i])
                if y_key_err not in data:
                    e.append(0.0)
                else:
                    e.append(data[y_key_err][i])
        return x, y, e, x_ul, y_ul

    def plot_this_instrument(self, config, y_axis):
        """Return if this instrument/axis should be plotted."""
        plot_this = config.get("plot_instrument", True)
        if y_axis not in config.get("plot_axis", []):
            plot_this = False
        return plot_this

    def get_marker_and_color(self, idx):
        """Return marker and color."""
        colors = plotting_utilities.get_color_list(len(self.data))
        color = (
            colors[idx]
            if self.config[idx].get("marker_color") is None
            else self.config[idx]["marker_color"]
        )
        marker = (
            plotting_utilities.get_marker_list()[idx]
            if self.config[idx].get("marker_type") is None
            else self.config[idx]["marker_type"]
        )
        return color, marker

    def plot_live_time_vs_phase_bin(
        self,
        instrument,
        phase_bins=10,
        file_type=".pdf",
        figure_dir="./figures/",
    ):
        """Plot live time histogram vs phase bin."""
        self._logger.info("Plotting live time histogram vs phase bin")

        ax = plotting_utilities.paper_figures(None, None)
        ax.set_xlim([0, 1])

        for idx, (data_instrument, data) in enumerate(self.data.items()):
            if self.config[idx].get("plot_live_time_histogram", False) is False:
                continue
            color, _ = self.get_marker_and_color(idx)

            x, y, _, _, _ = self._get_light_curve_in_mjd_limits(
                data, "live_time", "orbital phase",
            )
            if len(x) == 0:
                continue
            x = np.array(x)
            y = np.array(y)
            bin_edges = np.linspace(0, 1, phase_bins + 1)
            bin_width = bin_edges[1] - bin_edges[0]
            bin_heights, _ = np.histogram(x, bins=bin_edges, weights=y)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            plt.bar(
                bin_centers,
                bin_heights,
                width=bin_width*0.9,
                label=(
                    data_instrument
                    if self.config[idx].get("plot_label") is None
                    else self.config[idx]["plot_label"]
                ),
                color=color,
                alpha=0.7,
            )

        plt.xlabel("Phase")
        plt.ylabel("Live Time (h)")

        plt.legend()
        plotting_utilities.print_figure(
            f"Light-Curve-{self.binary['name']}-Live-Time-vs-Phase-Bin",
            file_type=file_type,
            figure_dir=figure_dir,
        )

    def plot_index_vs_flux(
        self,
        instrument,
        file_type=".pdf",
        figure_dir="./figures/",
    ):
        """Plot spectral index vs flux."""
        self._logger.info("Plotting spectral index vs flux")

        _ = plotting_utilities.paper_figures(None, None)

        for idx, (data_instrument, data) in enumerate(self.data.items()):
            if self.config[idx].get("plot_flux_vs_index", False) is False:
                continue
            color, _ = self.get_marker_and_color(idx)

            _, x, x_e, _, _ = self._get_light_curve_in_mjd_limits(
                data, "flux", "orbital phase",
            )
            _, y, y_e, _, _ = self._get_light_curve_in_mjd_limits(
                data, "index", "orbital phase",
            )
            if len(x) == 0 or len(y) == 0:
                continue
            x = np.array(x)
            y = np.array(y)
            plt.errorbar(
                x,
                y,
                xerr=x_e,
                yerr=y_e,
                label=(
                    data_instrument
                    if self.config[idx].get("plot_label") is None
                    else self.config[idx]["plot_label"]
                ),
                color=color,
                linestyle="none",
                marker='o',
                alpha=0.7,
            )

        plt.xlabel(self.config[0].get("flux_axis_label", ""))
        plt.ylabel(self.config[0].get("index_axis_label", ""))

        plt.legend()
        plotting_utilities.print_figure(
            f"Light-Curve-{self.binary['name']}-Flux-vs-Index",
            file_type=file_type,
            figure_dir=figure_dir,
        )

    def plot_distribution(
        self,
        instrument,
        y_axis="flux",
        file_type=".pdf",
        figure_dir="./figures/",
    ):
        """Plot distribution of y-axis data."""
        self._logger.info(f"Plotting {y_axis} distribution")

        _ = plotting_utilities.paper_figures(None, None)

        for idx, (data_instrument, data) in enumerate(self.data.items()):
            if self.config[idx].get("plot_1d_distribution", False) is False:
                continue
            color, _ = self.get_marker_and_color(idx)

            _, y, _, _, _ = self._get_light_curve_in_mjd_limits(
                data, y_axis, "orbital phase",
            )
            if len(y) == 0:
                continue
            y = np.array(y)
            plt.hist(
                y,
                bins=25,
                label=(
                    data_instrument
                    if self.config[idx].get("plot_label") is None
                    else self.config[idx]["plot_label"]
                ),
                color=color,
                alpha=0.7,
            )

        plt.xlabel(self.config[0].get(y_axis + "_axis_label", ""))
        plt.ylabel("Counts")

        plt.legend()
        plotting_utilities.print_figure(
            f"Light-Curve-{self.binary['name']}-{y_axis}-Distribution",
            file_type=file_type,
            figure_dir=figure_dir,
        )
