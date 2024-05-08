import numpy as np
import lagrange_functions


def rv_from_r0v0_ta(r0, v0, dt, mu):
    """Algorithm 2.3 in Orbital Mechanics for Engineering Students (3rd edition)

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
    f, g, fdot, gdot = lagrange_functions.f_and_g_ta(r0, v0, dt, mu)

    r = f * r0 + g * v0
    v = fdot * r0 + gdot * v0

    return r, v
