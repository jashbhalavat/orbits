import numpy as np


def rf_from_r0v0_ta(r0, v0, dt, mu):
    """
    This function computes the state vector (r,v)
    from the initial state vector (r0, v0) and the
    change in true anomaly.

    mu - gravitational parameter (km^3/s^2)
    r0 - initial position vector (km)
    v0 - initial velocity vector (km/s)
    dt - change in true anomaly (degrees)
    r - final position vector (km)
    v - final velocity vector (km/s)
    """

    