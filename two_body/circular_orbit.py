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
