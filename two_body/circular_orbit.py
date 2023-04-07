# Section 2.6 in "Orbital Mechanics for Engineering Students, Howard D. Curtis, 3rd edition"

import numpy as np


def calculate_speed(mu, R, altitude):
    """Calculate the speed of a satellite in circular orbit around a planet

    Inputs:
    mu - Gravitational parameter [km^3/s^2]
    R - Radius of planet [km]
    altitude - Altitude of satellite above planet surface [km]

    Outputs:
    velcity - velocity of satellite in orbit [km/s]
    """

    velocity = np.sqrt(mu / (R + altitude))

    return velocity


def calculate_period(mu, R, altitude):
    """Calculate the period of a satellite in circular orbit around a planet

    Inputs:
    mu - Gravitational parameter [km^3/s^2]
    R - Radius of planet [km]
    altitude - Altitude of satellite above planet surface [km]

    Outputs:
    period - Time taken by the satellite to complete one orbit [s]
    """

    period = (2 * np.pi) / np.sqrt(mu) * pow((R + altitude), 3 / 2)

    return period

def calculate_coverage(R, altitude):
    """Calculate the percent of surface that the satellite can see
    
    Inputs:
    R_e - Radius of planet [km]
    altitude - Altitude of satellite above planet surface [km]
    
    Output:
    phi - Maximum latitude [degrees]
    coverage - percentage of earth's surface visible to satellite
    """

    phi = np.rad2deg(np.arccos(R / (R + altitude)))
    coverage = 100 * (1 - np.cos(np.deg2rad(phi))) / 2

    return phi, coverage