"""Binary properties."""

def binary_properties():
    """
    List of properties of selected binary systems.

    Orbital periods is given in days.

    Returns
    -------
    dict
        Dictionary with binary properties.
    """
    return {
        "HESS J0632+057": {
            "name": "HESS J0632+057",
            "orbital_period": 317.3,
            "mjd_0": 54857.0,  # Bongiorno et al 2011
            "mjd_orbit_count": 52953.0,  # Adams et al 2021
        },
        "LS I +61 303": {
            "name": "LS I +61 303",
            "orbital_period": 26.496,  # Gregory et al. 2002
            "mjd_0": 2443366.775 - 2400000.5,
        },
        "LS I +61 303 superorbital": {
            "name": "LS I +61 303",
            "orbital_period": 1667.0,  # Gregory et al. 2002
            "mjd_0": 2443366.775 - 2400000.5,
        },
        "LS 5039": {
            "name": "LS 5039",
            "orbital_period": 3.90603,
            "mjd_0": 51942.59,
        },
    }
