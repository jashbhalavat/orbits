import numpy as np
from typing import List
import click
import json

import circular_orbit


@click.command()
@click.option(
    "-e", "--eccentricity", default=None, type=float, help="Orbit eccentricity"
)
@click.option(
    "-r_p",
    "--radius_of_perigee",
    default=None,
    type=float,
    help="Perigee altitude from planet surface [km]",
)
@click.option(
    "-r_a",
    "--radius_of_apogee",
    default=None,
    type=float,
    help="Apogee altitude from planet surface [km]",
)
def main(
    eccentricity,
    r_perigee,
    r_apogee
):
    """Decides which type of orbit function to call depending on the eccentricity

    TODO

    Inputs:
    Eccentricity
    Perigee Radius
    Apogee Radius

    Outputs:
    """

    # TODO - Parameters.json, library
    # TODO - Find a more generic way to compute elements

    # Variables
    G = 6.67430e-11  # Nm^2 / kg^s
    mu_earth = 398600  # km^3 / s^2
    R_earth = 6378  # km

    if eccentricity == 0:
        # Circular orbit
        if r_perigee != r_apogee:
            print("Eccentricity is zero (indicates circular orbit), but apogee and perigee radii are different")
            return
        speed = circular_orbit.calculate_speed(mu_earth, R_earth, r_perigee)
        period = circular_orbit.calculate_period(mu_earth, R_earth, r_perigee)
        phi, coverage = circular_orbit.calculate_coverage(R_earth, r_apogee)

        print("Speed:", speed, "km/s")
        print("Period:", period, "s")
        print("Maximum latitude:", phi, "degrees")
        print("Percentage of surface area visible to spacecraft:", coverage, "%")
    elif eccentricity > 0 and eccentricity < 1:
        # Elliptical orbit



if __name__ == "__main__":
    main()
