# Homework 1
# ASTR 8080

## Author

Serat Saad

## Installation

Ensure you have Python installed on your system. You also need to download and install the required libraries using pip:

pip install astropy matplotlib numpy

## Running the Script

To run this script, navigate to the directory containing `homework1.py` and execute it with Python. Use the following command in your terminal:

python homework1.py

## Features

- **Plot Planet Positions**: Generates plots for the positions of Mercury, Venus, Mars, Jupiter, and Saturn in ecliptic and equatorial coordinates for the year 2020 to 2030.
- **Lowest Airmass Calculations**: Find out the observations of these planets at the lowest airmass from APO.
- **Quasar Observability**: Determines which of the 1,111 quasars can be observed at the lowest airmass at 11PM MST from APO for any given month.

## Output

- Two PNG files showing the plotted positions of the planets:
  - `planet_position_ecliptic.png` for ecliptic coordinates.
  - `planet_position_equator.png` for equatorial coordinates.
- Output indicating the time and airmass of the lowest airmass observations for each planet.
- Output identifying the best quasar to observe at the lowest airmass for the specified month.
