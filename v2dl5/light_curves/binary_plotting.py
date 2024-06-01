"""
Binary light curve plotting.
"""

import logging

import matplotlib.pyplot as plt

import v2dl5.plotting.utilities as plotting_utilities


class BinaryLightCurvePlotter:
    """
    Binary light curve plotter

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

    def plot_flux_vs_time(self, time_axis="MJD", mjd_min=None, mjd_max=None, file_type=".pdf"):
        """
        Plot flux vs time (MJD or orbital phase)

        Parameters
        ----------
        mjd_min: float
            Minimum MJD value.
        mjd_max: float
            Maximum MJD value.

        """
        self._logger.info(f"Plotting flux vs {time_axis}")

        ax = plotting_utilities.paper_figures(None, None)
        colors = plotting_utilities.get_color_list(len(self.data))

        for idx, (instrument, data) in enumerate(self.data.items()):
            x, y, e, x_ul, y_ul = self._get_light_curve_in_MJD_limits(
                data, time_axis, mjd_min, mjd_max
            )

            if instrument in self.config and self.config["plot_instrument"] is False:
                continue
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
            plt.errorbar(
                x,
                y,
                e,
                None,
                label=(
                    instrument
                    if self.config[idx].get("plot_label") is None
                    else self.config[idx]["plot_label"]
                ),
                color=color,
                marker=marker,
                linestyle="none",
                fillstyle="none",
                linewidth=plotting_utilities.get_line_width(),
                markersize=plotting_utilities.get_marker_size(),
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
                    markersize=plotting_utilities.get_marker_size(),
                )

        plt.xlabel(self._get_time_axis_label(time_axis))
        if mjd_min is not None and mjd_max is not None:
            ax.set_xlim([mjd_min, mjd_max])
        # TMP - y-axis: use first entry in configuration dict
        plt.ylabel(self.config[0]["flux_axis_label"])
        plt.legend()

        plotting_utilities.print_figure(
            f"Light-Curve-{self.binary['name']}-Flux-vs-{time_axis.replace(' ', '-')}",
            file_type=file_type,
        )

    def _get_time_axis_label(self, time_axis):
        """
        Return time axis label

        """
        if time_axis == "orbital phase":
            return "Orbital phase"
        return "Modified Julian Date (MJD)"

    def _get_light_curve_in_MJD_limits(self, data, time_axis, mjd_min, mjd_max):
        """
        Get light curve restricted in MJD.

        Parameters
        ----------
        data : dict
            Light-curve data.
        time_axis : str
            Time axis (MJD or orbital phase)
        mjd_min : float
            Minimum MJD value.
        mjd_max : float
            Maximum MJD value.

        Returns
        -------
        list
            Light-curve data as three lists.

        """

        x = []
        y = []
        e = []
        x_ul = []
        y_ul = []
        for i in range(len(data["MJD"])):
            if mjd_min is not None and data["MJD"][i] < mjd_min:
                continue
            if mjd_max is not None and data["MJD"][i] > mjd_max:
                continue
            if time_axis == "MJD":
                w_x = data["MJD"][i]
            elif time_axis == "orbital phase":
                w_x = data["phase"][i]
            else:
                self._logger.error(f"Unknown time axis: {time_axis}")
                raise ValueError
            if "flux_ul" in data and data["flux_ul"][i] > 0:
                x_ul.append(w_x)
                y_ul.append(data["flux_ul"][i])
            elif "flux_ul" not in data or data["flux_ul"][i] < 0:
                x.append(w_x)
                y.append(data["flux"][i])
                e.append(data["flux_err"][i])
        return x, y, e, x_ul, y_ul


#
# def getNumberOfOrbits(fDataDict, PlotInstruments, OrbitalPeriod_BL, FixedMinMax = True ):
#     """
#     return number of orbits observed by the
#     data set.
#     """
#
#     # calculate number of orbits
#     min_MJD = 99999
#     max_MJD = 0
#     for I in PlotInstruments:
#         min_MJD = min(min_MJD, min(fDataDict[I]['MJD']))
#         max_MJD = max(max_MJD, max(fDataDict[I]['MJD']))
#     if FixedMinMax:
#         # hardwired min and max MJD for nicer plotting
#         # MJD min is fixed to first gamma-ray observations
#         min_MJD = 53080.
#         # 54540 is for VTS plotting
#         # min_MJD = 54540.
#         # value without 2019 data
#         # max_MJD = 58346.
#         # value to be used for 2019 data
#         max_MJD = 58485.
#         # adding NuSTAR data
#         # max_MJD = 58910.
#
#
#     MJD_firstOrbit = lightCurveAnalysisOrbitalPeriod.getMJDOrbitZeroPhase(
#         min_MJD, OrbitalPeriod_BL)
#     N_orbits = lightCurveAnalysisOrbitalPeriod.getNumberOfOrbits(
#         MJD_firstOrbit, max_MJD, OrbitalPeriod_BL)
#
#     print("Plot flux vs phase:")
#     print("\t orbital period (%.1f d)" % OrbitalPeriod_BL)
#     print("\t total number of orbits (%d)" % N_orbits)
#     if FixedMinMax:
#         print("\t MJD min (%.1f) max (%.1f) (hardwired!)" % (min_MJD, max_MJD))
#     else:
#         print("\t MJD min (%.1f) max (%.1f) " % (min_MJD, max_MJD))
#     print("\t MJD first orbit (%.1f)" % (MJD_firstOrbit))
#
#     return MJD_firstOrbit, N_orbits
#
#
# def plotLightCurve_fluxvsPhase_eachOrbit(
#         fDataDict,
#         PlotInstruments,
#         F,
#         lc_spline_bin_centers,
#         lc_spline_sv,
#         lc_spline_sv_err,
#         OrbitalPeriod_BL=315.,
#         OrbitalPeriodBins=20,
#         PlotRatio=False,
#         icrc2019Plots=False):
#     """plot flux vs orbital phase (separated in orbits, each orbit one plot)
#     """
#
#     if len(PlotInstruments) < 1:
#         return
#
#     # calculate min/max values in flux
#     f_ymin = 1.e9
#     f_ymax = 0.
#     for I in PlotInstruments:
#         f_ymin = min(f_ymin, min(
#             fDataDict[I]['flux']) - max(fDataDict[I]['flux_err']))
#         f_ymax = max(f_ymax, max(
#             fDataDict[I]['flux']) + max(fDataDict[I]['flux_err']))
#     if f_ymin > 0.:
#         f_ymin = 0.
#     f_ymax *= 1.2
#     f_ymin = max(-2.e-12, f_ymin)
#     f_ymin = max(-1.e-12, f_ymin)
#
#     MJD_firstOrbit, N_orbits = getNumberOfOrbits(
#         fDataDict, PlotInstruments,
#         OrbitalPeriod_BL, True)
#
#     if icrc2019Plots:
#         lightCurvePlottingUtilities.paper_figures(12, 6)
#     else:
#         lightCurvePlottingUtilities.paper_figures(8, 10)
#
#     if PlotInstruments[0].find('Swift XRT') >= 0:
#         colors = lightCurvePlottingUtilities.getColorList(len(PlotInstruments), "xray")
#     else:
#         colors = lightCurvePlottingUtilities.getColorList(len(PlotInstruments))
#
#     # sub plots: number of rows:
#     p_ncol = 3
#     if N_orbits > 12:
#         p_ncol = 4
#     if icrc2019Plots:
#         p_ncol = 6
#     p_nrow = math.ceil(N_orbits / p_ncol)
#
#     # plot one light curve per orbital period
#     for i in range(0, N_orbits):
#         # empty plot
#         axes = plt.subplot(p_nrow, p_ncol, i + 1)
#         # axis limits
#         axes.set_xlim([0., 1.])
#         if not PlotRatio:
#             axes.set_ylim([f_ymin, f_ymax])
#             if f_ymin < 0.:
#                 plt.axhline(y=0., linestyle=':', color='tab:gray')
#         else:
#             axes.set_ylim([0., 6.])
#
#         orbit_min = MJD_firstOrbit + i * OrbitalPeriod_BL
#         orbit_max = MJD_firstOrbit + (i + 1) * OrbitalPeriod_BL - 1.
#
#         OrbitPhStr = "Orbit %d:\nMJD %d" % (i + 1, orbit_min)
#         plt.text(
#             0.1,
#             0.83,
#             OrbitPhStr,
#             transform=axes.transAxes,
#             size='smaller')
#
#         r_av = 0.
#         r_n = 0.
#         c = 0
#         for I in PlotInstruments:
#             x = []
#             y = []
#             ex = []
#             ey = []
#             mj_min = []
#             mj_max = []
#             for p in range(len(fDataDict[I]['phaseN'])):
#                 if fDataDict[I]['MJD'][p] >= orbit_min \
#                         and fDataDict[I]['MJD'][p] < orbit_max:
#
#                     x.append(fDataDict[I]['phase'][p])
#                     # plot time periods as 'x-errors'
#                     ex.append(fDataDict[I]['MJD_err'][p]/OrbitalPeriod_BL)
#                     if not PlotRatio:
#                         y.append(fDataDict[I]['flux'][p])
#                         ey.append(fDataDict[I]['flux_err'][p])
#                     else:
#                         y.append(fDataDict[I]['flux'][p] /
#                                  F(1. + fDataDict[I]['phase'][p]))
#                         ey.append(fDataDict[I]['flux_err'][p] /
#                                   F(1. + fDataDict[I]['phase'][p]))
#                         if ey[-1] > 0.:
#                             r_av += y[-1] / (ey[-1] * ey[-1])
#                             r_n += 1. / (ey[-1] * ey[-1])
#                     mj_min.append(
#                         fDataDict[I]['MJD'][p] -
#                         fDataDict[I]['MJD_err'][p])
#                     mj_max.append(
#                         fDataDict[I]['MJD'][p] +
#                         fDataDict[I]['MJD_err'][p])
#
#             # average plot only for VERITAS and XRT
#             if I.find("VERITAS") >= 0 or I.find("XRT") >= 0:
#                 glabel = I + " (average)"
#                 plt.plot(
#                     lc_spline_bin_centers,
#                     lc_spline_sv,
#                     color='tab:gray',
#                     linestyle='--',
#                     linewidth=lightCurvePlottingUtilities.getLineWidth())
#                 plt.fill_between(
#                     lc_spline_bin_centers,
#                     lc_spline_sv - lc_spline_sv_err,
#                     lc_spline_sv + lc_spline_sv_err,
#                     color='tab:gray',
#                     linestyle='--',
#                     linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                     alpha=0.3)
#
#             # plot points
#             if len(x):
#                 pLabel=I
#                 if pLabel.find( 'XRT' )>=0:
#                     pLabel="$\it{Swift}$-XRT"
#                 plt.errorbar(
#                     x,
#                     y,
#                     ey,
#                     ex,
#                     color=colors[c],
#                     marker='o',
#                     linestyle='none',
#                     label=pLabel,
#                     linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                     markersize=lightCurvePlottingUtilities.getMarkerSize())
#             c = c + 1
#
#         plt.locator_params(axis='x', nbins=4)
#         plt.xlabel(
#             lightCurvePlottingUtilities.getOrbitalPhaseAxisString(OrbitalPeriod_BL),
#             fontsize=8)
#         plt.ylabel(lightCurvePlottingUtilities.getFluxAxisString(PlotInstruments[0]),
#         fontsize=8)
#         plt.legend(loc=2)
#         plt.legend(prop={'size': 5}, framealpha=0.2)
#
#         if PlotRatio:
#             plt.axhline(y=1, linestyle=':', color='tab:gray')
#             # plot average ratio + error
#             if r_n > 0.:
#                 r_av /= r_n
#                 r_er = math.sqrt(1. / r_n)
#                 r_x = [0, 1]
#                 r_y1 = [r_av - r_er, r_av - r_er]
#                 r_y2 = [r_av + r_er, r_av + r_er]
#                 axes.fill_between(
#                     r_x, r_y1, r_y2, facecolor='tab:gray', alpha=0.3)
#                 plt.axhline(y=r_av, linestyle=':', color='b')
#                 r_str = 'Ratio to average: %.2f+-%.2f' % (r_av, r_er)
#                 # plt.text(0.1, 0.1, r_str)
#                 print('Orbit %d: %s' % (i,r_str))
#             plt.ylabel('Ratio to average')
#
#     ratiostr = ""
#     if PlotRatio:
#         ratiostr = "Ratio"
#     lightCurvePlottingUtilities.printFigure(
#         getPrintInstrumentName(PlotInstruments) + "-HESSJ0632p057-LC-phaseFolded-%dd-perOrbit%s" %
#         (OrbitalPeriod_BL, ratiostr))
#
#
# def plotLightCurve_fluxvsPhase_inOrbits(
#         fDataDict,
#         PlotInstruments,
#         F,
#         lc_spline_bin_centers,
#         lc_spline_sv,
#         lc_spline_sv_err,
#         OrbitalPeriod_BL=315.,
#         OrbitalPeriodBins=20,
#         PlotVariable=None):
#     """plot flux vs orbital phase (separated in orbits, all in one plot)
#     """
#
#     MJD_firstOrbit, N_orbits = getNumberOfOrbits(
#         fDataDict, PlotInstruments,
#         OrbitalPeriod_BL, False)
#
#     lightCurvePlottingUtilities.paper_figures(4, 4)
#     colors = lightCurvePlottingUtilities.getColorList(N_orbits)
#     if N_orbits < 6:
#         markers = lightCurvePlottingUtilities.getMarkerList()
#     else:
#         markers = ['o']*N_orbits
#
#     if len(PlotInstruments) < 1:
#         return
#
#     # plot average and interpolated light curves
#     glabel = PlotInstruments[0] + " (average)"
#     if lc_spline_bin_centers and len(lc_spline_bin_centers) > 0:
#         plt.plot(
#             lc_spline_bin_centers,
#             lc_spline_sv,
#             color='tab:gray',
#             linestyle='--',
#             linewidth=lightCurvePlottingUtilities.getLineWidth())
#         plt.fill_between(lc_spline_bin_centers,
#                         lc_spline_sv - lc_spline_sv_err,
#                         lc_spline_sv + lc_spline_sv_err,
#                         color='tab:gray', linestyle='--',
#                         linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                         alpha=0.3)
#
#     # plot one light curve per orbital period
#     for i in range(0, N_orbits):
#         x = []
#         y = []
#         ex = []
#         mj = []
#         orbit_min = MJD_firstOrbit + i * OrbitalPeriod_BL
#         orbit_max = MJD_firstOrbit + (i + 1) * OrbitalPeriod_BL - 1.
#         for j in range(len(PlotInstruments)):
#             i_plotValue, i_plotError = lightCurvePlottingUtilities.getPlottingVariable(
#                 PlotVariable, j)
#
#             for p in range(len(fDataDict[PlotInstruments[j]]['phaseN'])):
#                 if fDataDict[PlotInstruments[j]]['MJD'][p] >= orbit_min \
#                         and fDataDict[PlotInstruments[j]]['MJD'][p] < orbit_max + 1.:
#
#                     x.append(fDataDict[PlotInstruments[j]]['phase'][p])
#                     y.append(fDataDict[PlotInstruments[j]][i_plotValue][p])
#                     ex.append(fDataDict[PlotInstruments[j]][i_plotError][p])
#                     mj.append(fDataDict[PlotInstruments[j]]['MJD'][p])
#
#         if len(mj):
#             OrbitPhStr = "MJD %d - %d" % (
#                 orbit_min, orbit_max)
#             plt.errorbar(
#                 x,
#                 y,
#                 ex,
#                 None,
#                 color=colors[i],
#                 marker=markers[i],
#                 linestyle='none',
#                 label=OrbitPhStr,
#                 linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                 markersize=lightCurvePlottingUtilities.getMarkerSize())
#
#     plt.xlabel(
#         lightCurvePlottingUtilities.getOrbitalPhaseAxisString(OrbitalPeriod_BL) )
#     plt.ylabel(
#         lightCurvePlottingUtilities.getFluxAxisString(
#         PlotInstruments[0],PlotVariable) )
#     if PlotInstruments[0].find('Optical') < 0:
#         plt.axhline(y=0, linestyle=':')
#     if getPrintInstrumentName(PlotInstruments).find( "fwhm" ) < 0 and \
# getPrintInstrumentName(PlotInstruments).find( "ew" ) < 0:
#         plt.legend(prop={'size': 10}, framealpha=0.1)
#
#     lightCurvePlottingUtilities.printFigure(
#         getPrintInstrumentName(PlotInstruments) +
#         "-HESSJ0632p057-LC-phaseFolded-%dd-Orbits" %
#         OrbitalPeriod_BL)
#
#
# def plotLightCurve_fluxvsPhase(
#         fDataDict,
#         PlotInstruments,
#         lc_spline_bin_centers,
#         lc_spline_sv,
#         lc_spline_sv_err,
#         OrbitalPeriod_BL=315.,
#         MJDPlotmin=0,
#         MJDPlotMax=90000,
#         PlotVariable=None,
#         icrc2019Plots=False,
#         donotplotaverage=False):
#     """plot flux vs orbital phase
#     """
#     ax = lightCurvePlottingUtilities.paper_figures(4, 4)
#     if PlotInstruments[0].find('Swift XRT') >= 0:
#         colors = lightCurvePlottingUtilities.getColorList(
#                     len(PlotInstruments), "xray" )
#     else:
#         colors = lightCurvePlottingUtilities.getColorList(
#                         len(PlotInstruments) )
#
#     markers = lightCurvePlottingUtilities.getMarkerList(icrc2019Plots)
#
#     if len(PlotInstruments) < 1:
#         return
#
#     c = 0
#     for i in range(len(PlotInstruments)):
#
#         i_plotValue, i_plotError = lightCurvePlottingUtilities.getPlottingVariable(
#             PlotVariable, i)
#         xp = []
#         yp = []
#         ee = []
#         ye = []
#         for p in range(len(fDataDict[PlotInstruments[i]]['MJD'])):
#             if fDataDict[PlotInstruments[i]
#                  ]['MJD'][p] > MJDPlotmin and
#                   fDataDict[PlotInstruments[i]]['MJD'][p] < MJDPlotMax:
#                 xp.append(fDataDict[PlotInstruments[i]]['phase'][p])
#                 yp.append(fDataDict[PlotInstruments[i]][i_plotValue][p])
#                 ee.append(fDataDict[PlotInstruments[i]][i_plotError][p])
#
#                 # plot time periods as 'x-errors'
#                 if 'phase_err' in fDataDict[PlotInstruments[i]] and \
#                         len(fDataDict[PlotInstruments[i]]['phase_err']) == 2 and \
#                         len(fDataDict[PlotInstruments[i]]['phase_err'][0]) \
#                         == len(fDataDict[PlotInstruments[i]]['phase']) and \
#                         len(fDataDict[PlotInstruments[i]]['phase_err'][1]) \
#                         == len(fDataDict[PlotInstruments[i]]['phase']):
#                     ye.append(0.5*(fDataDict[PlotInstruments[i]]['phase_err'][0][p])+\
#                             fDataDict[PlotInstruments[i]]['phase_err'][1][p])
#                 else:
#                     ye.append(0.)
#
#         if len(xp) > 0:
#             ff='full'
#             si=lightCurvePlottingUtilities.getMarkerSize()
#             li=lightCurvePlottingUtilities.getLineWidth()
#             if getPrintInstrumentName(PlotInstruments).find( "XRay" ) >= 0:
#                 ff='none'
#                 ff='full'
#                 si=0.75*lightCurvePlottingUtilities.getMarkerSize()
#                 li=0.5*lightCurvePlottingUtilities.getLineWidth()
#             pLabel=PlotInstruments[i]
#             if pLabel.find( 'XRT' )>=0:
#                 pLabel="$\it{Swift}$-XRT"
#             plt.errorbar(
#                 xp, yp, ee,
#                 ye,
#                 color=colors[c], marker=markers[c], linestyle='none',
#                 label=pLabel,
#                 linewidth=li, fillstyle=ff,
#                 markersize=si )
#         c = c + 1
#
#     glabel = PlotInstruments[0] + " (average)"
#     if lc_spline_bin_centers and len(lc_spline_bin_centers) > 0 \
#             and not icrc2019Plots \
#             and not donotplotaverage:
#         plt.plot(
#             lc_spline_bin_centers,
#             lc_spline_sv,
#             color='tab:gray',
#             linestyle='--',
#             linewidth=lightCurvePlottingUtilities.getLineWidth())
#         plt.fill_between(lc_spline_bin_centers,
#                          lc_spline_sv - lc_spline_sv_err,
#                          lc_spline_sv + lc_spline_sv_err,
#                          color='tab:gray', linestyle='-',
#                          linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                          alpha=0.3)
#     plt.xlabel(
#         lightCurvePlottingUtilities.getOrbitalPhaseAxisString(OrbitalPeriod_BL))
#     plt.ylabel(
#         lightCurvePlottingUtilities.getFluxAxisString(
#             PlotInstruments[0],
#             PlotVariable))
#     if PlotInstruments[0].find('Optical') < 0:
#         plt.axhline(y=0, linestyle=':')
#         plt.legend()
#     #elif PlotVariable and len(PlotVariable) == 1:
#     #    PlotInstruments[0] += "-" + PlotVariable[0]
#
#     lightCurvePlottingUtilities.printFigure(
#         getPrintInstrumentName(PlotInstruments)
#  + "-HESSJ0632p057-LC-phaseFolded-%dd" % OrbitalPeriod_BL)
#
#
# def plotAverageLightCurve_fluxvsPhase(
#         fDataDict,
#         PlotInstruments,
#         OrbitalPeriod_BL=315.,
#         OrbitalPeriodBins=20):
#     """
#     plot average flux vs orbital phase
#
#     calculates also average light curves
#     """
#
#     print("Plot phase binned averaged light curve:")
#     print("\t orbital phase (%.1f d)" % OrbitalPeriod_BL)
#     print("\t number of phase bins (%d)" % OrbitalPeriodBins)
#     print("\t calculating average light curve for ", PlotInstruments)
#     if len(PlotInstruments) < 1:
#         return
#
#     # copy all data into one set of arrays
#     # averaged light curve is calculated from
#     # light curve bins of all data
#     i_MJD = []
#     i_flux = []
#     i_flux_err = []
#     for I in PlotInstruments:
#         i_MJD.extend(fDataDict[I]['MJD'])
#         i_flux.extend(fDataDict[I]['flux'])
#         i_flux_err.extend(fDataDict[I]['flux_err'])
#
#     lightCurvePlottingUtilities.paper_figures(4, 4)
#
#     # calculate phase binned average light curve
#     lc_bincenters, lc_mean, lc_std = lightCurveAverageing.calculateAverageLightCurve(
#         i_MJD, i_flux, i_flux_err, OrbitalPeriod_BL, OrbitalPeriodBins)
#
#     # bin width
#     lc_binw = np.repeat(0.5 / OrbitalPeriodBins, len(lc_bincenters))
#     plt.errorbar(lc_bincenters, lc_mean, lc_std, lc_binw,
#                  color='b', linestyle='none',
#                  linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                  markersize=lightCurvePlottingUtilities.getMarkerSize(),
#                  marker='o', label='average')
#
#     # calculate cubic spline smoothed average light curve
#     F, lc_spline_bin_centers, lc_spline_sv = \
#         lightCurveAverageing.smoothCubeSplineLightCurve(
#             i_MJD, i_flux, i_flux_err,
#             OrbitalPeriod_BL, OrbitalPeriodBins,
#             True, 0.)
#     # get errors on the same
#     F_temp, lc_spline_bin_centers, lc_spline_sv_err = \
#         lightCurveAverageing.smoothCubeSplineLightCurve(
#             i_MJD, i_flux, i_flux_err,
#             OrbitalPeriod_BL, OrbitalPeriodBins,
#             True, 1.)
#     lc_spline_sv_err = list(
#         np.array(lc_spline_sv_err) -
#         np.array(lc_spline_sv))
#
#     plt.plot(lc_spline_bin_centers, lc_spline_sv, color='r', linestyle='--',
#              linewidth=lightCurvePlottingUtilities.getLineWidth(),
#              label='spline average')
#
#     # data set used for averaging
#     wP = []
#     wP = [lightCurveAnalysisOrbitalPeriod.getOrbitalPhase(
#         x, OrbitalPeriod_BL) for x in i_MJD]
#     plt.errorbar(wP, i_flux, i_flux_err, None,
#                  color='g', linestyle='none',
#                  linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                  markersize=lightCurvePlottingUtilities.getMarkerSize(),
#                  marker='o', label='average')
#
#     plt.fill_between(lc_spline_bin_centers,
#                      lc_spline_sv - lc_spline_sv_err,
#                      lc_spline_sv + lc_spline_sv_err,
#                      color='r', linestyle='--',
#                      linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                      alpha=0.3)
#
#     plt.xlabel(
#         lightCurvePlottingUtilities.getOrbitalPhaseAxisString(OrbitalPeriod_BL))
#     plt.ylabel(lightCurvePlottingUtilities.getFluxAxisString(PlotInstruments[0]))
#     plt.legend()
#     plt.axhline(y=lc_spline_sv[0], linestyle=':')
#
#     lightCurvePlottingUtilities.printFigure(
#         getPrintInstrumentName(PlotInstruments) +
#         "-HESSJ0632p057-LC-phaseFolded-%dd-Average" %
#         OrbitalPeriod_BL)
#
#     return F, lc_spline_bin_centers, lc_spline_sv, lc_spline_sv_err
#
#
# def plotLightCurve_fluxvsMJD_XandGray(
#         fDataDict, mjd_min, mjd_max,
#         icrc2019Plots=False,
#         yaxis_min = -0.9, yaxis_max=7.5,
#         fColorDict=None,
#         plot_title=None,
#         convert_erg=False):
#     """plot flux vs MJD with two different y-axis
#
#     note: fixed limits for the y-axis
#     """
#
#     lightCurvePlottingUtilities.paper_figures(4, 4, 1, False)
#     # quick and dirty fix to get consistent colors and markers
#     if len(fColorDict) == 5:
#         colors = lightCurvePlottingUtilities.getColorList(3)
#         colors[-1] = 'black'
#         colors.append( 'gray' )
#     else:
#         colors = lightCurvePlottingUtilities.getColorList(len(fColorDict))
#     markers = lightCurvePlottingUtilities.getMarkerList(icrc2019Plots)
#     markers[3]='+'
#     markers[4]='x'
#
#     fig, ax1 = plt.subplots()
#     ax1_y = ax1.twinx()
#     ax1_x = ax1.twiny()
#     if convert_erg:
#         ax1.set_ylim(ymin=yaxis_min,ymax=yaxis_max*1.8)
#     else:
#         ax1.set_ylim(ymin=yaxis_min,ymax=yaxis_max)
#     ax1_y.set_ylim(ymin=yaxis_min,ymax=yaxis_max)
#
#     c = 0
#     for key, fData in fDataDict.items():
#         print('Plotting %s (%d data points)' % (key, len(fData)))
#
#         # scale all data by 1.e-12
#         # (scale is added to the legend text)
#         y = np.asarray(fData['flux']) / 1.e-12
#         y_err = np.asarray(fData['flux_err']) / 1.e-12
#
#         if key in fColorDict:
#             c = fColorDict[key]
#
#         if len(y) == 0:
#             c += 1
#             continue
#
#         # left labelled data
#         if "VERITAS" in key or \
#                 "HESS" in key or \
#                 "MAGIC" in key:
#             ln1 = ax1.errorbar(
#                 fData['MJD'],
#                 y,
#                 y_err,
#                 fData['MJD_err'],
#                 color=colors[c],
#                 marker=markers[c],
#                 label=key,
#                 linestyle='none',
#                 fillstyle='full',
#                 linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                 markersize=lightCurvePlottingUtilities.getMarkerSize())
#             ax1_y.plot(np.nan, '-r',
#                      color=colors[c], marker=markers[c], label=key,
#                      linestyle='none', fillstyle='full',
#                      linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                      markersize=lightCurvePlottingUtilities.getMarkerSize())
#         # right labelled data
#         else:
#             pLabel=key
#             if key.find('XRT')>=0:
#                 pLabel="$\it{Swift}$-XRT"
#             ln2 = ax1_y.errorbar(
#                 fData['MJD'],
#                 y,
#                 y_err,
#                 fData['MJD_err'],
#                 color=colors[c],
#                 marker=markers[c],
#                 label=pLabel,
#                 linestyle='none',
#                 fillstyle='full',
#                 linewidth=lightCurvePlottingUtilities.getLineWidth(),
#                 markersize=lightCurvePlottingUtilities.getMarkerSize())
#         c += 1
#
#     ax1.axis(xmin=mjd_min,xmax=mjd_max)
#     ax1.locator_params(axis='x', tight=True, nbins=4)
#     ax1.set_xlabel('Modified Julian Day (MJD)')
#     # left label: assume always gamma rays
#     ax1.set_ylabel(lightCurvePlottingUtilities.getFluxAxisString("VERITAS",
#                                                                  None,
#                                                                  "10$^{-12} \\times$",
#                                                                  convert_erg))
#
#     if plot_title:
#         ty=0.9*yaxis_max
#         if convert_erg:
#             ty*=1.8
#         ax1.text(
#               mjd_min+0.75*(mjd_max-mjd_min),
#               ty,
#               plot_title,
#               fontsize=10,
#               weight='bold',
#               )
#
#     # right label: assume always X-rays
#     if 'Swift XRT' in fDataDict or 'NuSTAR' in fDataDict:
#         ax1_y.set_ylabel(lightCurvePlottingUtilities.getFluxAxisString(
#   "Swift XRT", None, "10$^{-12} \\times$"))
#
#     ax1_y.legend(loc=2)
#     if yaxis_min < 0.:
#         plt.axhline(y=0, linestyle=':')
#     lightCurvePlottingUtilities.align_yaxis(ax1, 0, ax1_y, 0)
#
#     # orbital phase axis
#     # (note: use default orbital period here
#     # as defined in lightCurveAnalysisOrbitalPeriod
#     phase_ticks = []
#     label_format = '%.2f'
#     for t in ax1_x.get_xticks():
#         mjd=mjd_min+t*(mjd_max-mjd_min)
#         phase_ticks.append(label_format % (
#             lightCurveAnalysisOrbitalPeriod.getOrbitalPhase(mjd),))
#     # remove first and last label from plotting
#     if len(phase_ticks)>0:
#         phase_ticks[0]=None
#         phase_ticks[-1]=None
#     ax1_x.set_xticklabels(phase_ticks)
#     ax1_x.set_xlabel('orbital phase')
#
#     lightCurvePlottingUtilities.printFigure('XG-HESSJ0632p057-LC')
#
#
