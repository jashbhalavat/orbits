import numpy as np
from typing import List
import click
import json

import circular_orbit


@click.command()
@click.option(
    "-e", "--eccentricity", default=None, type=float, help="Orbit eccentricity"
)
@click.option("-a", "--semimajor_axis", default=None, type=float, help="Semimajor axis")
@click.option("-i", "--inclination", default=None, type=float, help="Orbit inclination")
@click.option(
    "-o",
    "--longitude_of_ascending_node",
    default=None,
    type=float,
    help="Longitude of ascending node",
)
@click.option(
    "-w",
    "--argument_of_periapsis",
    default=None,
    type=float,
    help="Argument of periapsis",
)
@click.option("-t", "--true_anomaly", default=None, type=float, help="True anomaly")
def main(
    eccentricity,
    semimajor_axis,
    inclination,
    longitude_of_ascending_node,
    argument_of_periapsis,
    true_anomaly,
):
    """Decides which type of orbit function to call depending on the eccentricity

    TODO

    Inputs:
    Eccentricity
    Semimajor Axis
    Inclination
    Longitude of Ascending Node
    Argument of Periapsis
    True Anomaly

    Outputs:
    """

    # TODO - Parameters.json, library

    # Variables
    G = 6.67430e-11  # Nm^2 / kg^s
    mu_earth = 398600  # km^3 / s^2
    R_earth = 6378  # km

    if eccentricity == 0:
        # Circular orbit
        speed = circular_orbit.calculate_speed(mu_earth, R_earth, semimajor_axis)
        period = circular_orbit.calculate_period(mu_earth, R_earth, semimajor_axis)
        phi, coverage = circular_orbit.calculate_coverage(R_earth, semimajor_axis)

        print("Speed:", speed, "km/s")
        print("Period:", period, "s")
        print("Maximum latitude:", phi, "degrees")
        print("Percentage of surface area visible to spacecraft:", coverage, "%")
    elif eccentricity > 0 and eccentricity < 1:
        # Elliptical orbit


if __name__ == "__main__":
    main()
