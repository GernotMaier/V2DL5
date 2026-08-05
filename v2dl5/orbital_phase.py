"""Orbital phase calculations."""

import math

import numpy as np
from astropy.time import Time


def get_orbital_phase_from_iso_time(iso_time, orbital_period, mjd_0, phase_reduce=True):
    """
    Calculate orbital phase for a given time.

    Parameters
    ----------
    iso_time: str
        Time in ISO format (e.g., '2023-10-01T00:00:00')
    orbital_period: float
        Orbital period (in units of days)
    mjd_0: float
        Reference MJD (in units of days)
    phase_reduce: bool
        Reduce phase to interval [0,1]

    Returns
    -------
    float or np.ndarray
        Orbital phase

    """
    return get_orbital_phase(
        mjd=Time(iso_time, scale="utc").mjd,
        orbital_period=orbital_period,
        mjd_0=mjd_0,
        phase_reduce=phase_reduce,
    )


def get_orbital_phase(mjd, orbital_period, mjd_0, phase_reduce=True):
    """
    Calculate orbital phase for a given MJD.

    Parameters
    ----------
    mjd: float
        MJD
    orbital_period: float
        Orbital period (in units of days)
    mjd_0: float
        Reference MJD (in units of days)
    phase_reduce: bool
        Reduce phase to interval [0,1]

    Returns
    -------
    float
        Orbital phase

    """
    mjd = np.asarray(mjd)
    scalar_input = False
    if mjd.ndim == 0:
        mjd = mjd[None]  # Makes x 1D
        scalar_input = True

    phase = (mjd - mjd_0) / orbital_period

    if phase_reduce:
        phase -= phase.astype(int)

    for i in range(len(phase)):
        if phase[i] < 0:
            phase[i] = 1.0 + phase[i]

    if scalar_input:
        return np.squeeze(phase)
    return phase


def get_orbital_phase_range(mjd_min, mjd_max, phase_mean, orbital_period, mjd_0, upper_error=True):
    """
    Calculate width of observation bin in orbital phase.

    Parameters
    ----------
    mjd_min: float
        Start time of observation
    mjd_max: float
        End time of observation
    phase_mean: float
        Mean phase of observation
    orbital_period: float (in units of days)
        Orbital period
    mjd_0: float
        Reference MJD
    upper_error: bool
        Upper error

    Returns
    -------
    float
        Width of observation bin in orbital phase
    """
    ph_min = get_orbital_phase(mjd=mjd_min, orbital_period=orbital_period, mjd_0=mjd_0)
    ph_max = get_orbital_phase(mjd=mjd_max, orbital_period=orbital_period, mjd_0=mjd_0)
    if abs(ph_max - ph_min) < 1.0e-3:
        ph_err = 0.0
    elif ph_max > ph_min:
        ph_err = 0.5 * (ph_max - ph_min)
    else:
        ph_err = 0.5 * (1.0 + ph_max - ph_min)

    # check if error reaches over '0' or '1' -> simply cut
    # (not entirely correct)
    if upper_error:
        if phase_mean + ph_err > 1.0:
            ph_err = 1.0 - phase_mean
    else:
        if phase_mean - ph_err < 0.0:
            ph_err = phase_mean
    return ph_err


def get_orbit_number(mjd, orbital_period, mjd_orbit_count):
    """
    Return orbit number for a given MJD since MJD0.

    Parameters
    ----------
    mjd: float
        MJD
    orbital_period: float
        Orbital period (in units of days)
    mjd_orbit_count: float
        Reference MJD for counting orbits

    Returns
    -------
    int
        Orbit number
    """
    return math.ceil((mjd - mjd_orbit_count) / orbital_period)


# def getMJDOrbitZeroPhase(mjd, object="HESS J0632+057", orbital_period=317.3):
#     """
#     Return MJD of phase 0 for the given orbital period period.
#
#     Parameters:
#     -----------
#     mjd: float
#         MJD
#     object: str
#         name of the object
#     orbital_period: float
#         orbital period
#
#     Returns:
#     --------
#     MJD0: float
#         MJD0 for the object
#
#     """
#     norbit = math.floor((mjd - get_mjd0()) / orbital_period)
#
#     return get_mjd0() + norbit * orbital_period
