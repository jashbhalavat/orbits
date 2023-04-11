# Section 2.7 in "Orbital Mechanics for Engineering Students, Howard D. Curtis, 3rd edition"

import numpy as np
def calculate_eccentricity(R, r_p, r_a):
    """Calculate the eccentricity of an elliptical orbit
    
    Inputs:
    R - Radius of planet [km]
    r_p - Orbit Perigee from planet surface [km]
    r_a - Orbit Apogee from planet surface [km]
    
    Output:
    e - Eccentricity of orbit
    """

    if r_a > r_p:
        e = (r_a - r_p) / (r_a + r_p + 2*R)
    else:
        print("Orbit perigee is greater than orbit apogee")

    return e

def calculate_angularmomentum(R, mu, r_p, e):
    """Calculate angular momentum
    
    Inputs:
    R - Radius of planet [km]
    mu - Gravitational parameter [km^3/s^2]
    r_p - Orbit Perigee from planet surface [km]
    e - Eccentricity of orbit
    
    Outputs:
    h - Angular momentum [km^2/s]
    """

    h = np.sqrt((R + r_p) * mu * (1 + e))

    return h

def calculate_apvel(h, r):
    """Calculate the velocity at apogee or perigee
    
    Inputs:
    R - Radius of planet [km]
    h - Angular momentum [km^2/s]
    r - Orbit apogee/perigee [km]
    
    Outputs:
    v_ap = Velocity at apogee/perigee [km/s]
    """

    v_ap = h / (R + r)

    return v_ap

def calculate_semimajoraxis(R, r_a, r_p):
    """Calculate the semimajor axis of the orbit
    
    Inputs:
    R - Radius of planet [km]
    r_a - Orbit apogee [km]
    r_p - Orbit perigee [km]
    
    Outputs:
    a - Semimajor axis [km]
    """

    a = (2 * R + r_a + r_p) / 2

    return a

def calculate_period(mu, a):
    """Calculate the period of an orbit
    
    Inputs:
    mu - Gravitational parameter [km^3/s^2]
    a - Semimajor axis [km]
    
    Outputs:
    T - Orbit period [s]
    """

    T = 2*np.pi / (np.sqrt(mu)) * pow(a, 3/2)

    return T
