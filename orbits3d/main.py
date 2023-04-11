import numpy as np
from typing import List
import click
import json


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