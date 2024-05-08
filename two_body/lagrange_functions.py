import numpy as np


def f_and_g_ta(r0, v0, dt, mu):
    """Algorithm D.7 in Orbital Mechanics for Engineering Students (3rd edition)

    This function calculates the Lagrange f and g coefficients from the
    change in true anomaly since time t0.

    mu - gravitational parameter (km^3/s^2)
    dt - change in true anomaly (degrees)
    r0 - position vector at time t0 (km)
    v0 - velocity vector at time t0 (km/s)
    h - angular momentum (km^2/s)
    vr0 - radial component of v0 (km/s)
    r - radial position after the change in true anomaly
    f - the Lagrange f coefficient (dimensionless)
    g - the Lagrange g coefficient (s)
    """

    h = np.linalg.norm(np.cross(r0, v0))
    r0 = np.linalg.norm(r0)
    vr0 = np.dot(v0, r0) / r0
    s = np.sin(np.deg2rad(dt))
    c = np.cos(np.deg2rad(dt))

    # Eq. 2.152
    r = (pow(h, 2) / mu) * (1 / (1 + (pow(h, 2) / mu / r0 - 1) * c - h * vr0 * s / mu))

    # Eq. 2.158a and 2.158b
    f = 1 - mu * r * (1 - c) * pow(h, 2)
    g = r * r0 * s / h

    return f, g


def fdot_and_gdot_ta(r0, v0, dt, mu):
    """Algorithm D.7 in Orbital Mechanics for Engineering Students (3rd edition)

    This function calculates the Lagrange f and g coefficients from the
    change in true anomaly since time t0.

    mu - gravitational parameter (km^3/s^2)
    dt - change in true anomaly (degrees)
    r0 - position vector at time t0 (km)
    v0 - velocity vector at time t0 (km/s)
    h - angular momentum (km^2/s)
    vr0 - radial component of v0 (km/s)
    r - radial position after the change in true anomaly
    fdot - time derivative of the Lagrange f coefficient (1/s)
    gdot - time derivative of the Lagrange g coefficient (dimensionless)
    """
    h = np.linalg.norm(np.cross(r0, v0))
    r0 = np.linalg.norm(r0)
    vr0 = np.dot(v0, r0) / r0
    s = np.sin(np.deg2rad(dt))
    c = np.cos(np.deg2rad(dt))

    # Eq. 2.158c and 2.158d
    fdot = mu / h * (vr0 / h * (1 - c) - s / r0)
    gdot = 1 - mu * r0 / pow(h, 2) * (1 - c)

    return fdot, gdot
