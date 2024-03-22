# M. Kaldor
# v1 1/25/2024
# ASTR 8080 HW0

# to import this from another directory:
# import sys
# sys.path.insert(0, '../week1/')
# import hw0
# from hw0 import function

# IMPORT BLOCK
###############################
###############################
import math

import numpy as np
import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt
from matplotlib import rc
import time
import pdb
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time


# FUNCTIONS
###############################
###############################
def cart_conv(ra, dec):
    """
    This is a function that converts a right ascension and declination from spherical to Cartesian coordinates.

    Args: ra [hours, minutes, seconds]
        dec [degrees, minutes, seconds]

    Returns: x [float]
        y [float]
        z [float]

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: converting into Cartesian units
    #
    #
    # CALLING SEQUENCE: cart_conv("17h18m00s", "49:04:02")
    #
    #
    # INPUTS: ra in hms, dec in degrees minutes seconds
    #
    #
    # -------------------------------------------------------------
    c = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg))
    c.representation_type = "cartesian"
    return c

###############################
def xyz_eq(ra, dec):
    """
    This is a function converts the ra and dec from spherical into xyz Cartesian coords via the class equations.

    Args:

    Returns:

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE:
    #
    #
    # CALLING SEQUENCE:
    #
    #
    # INPUTS:
    #
    # -------------------------------------------------------------
    #
    c = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg))
    ra_eq = c.ra.deg
    dec_eq = c.dec.deg
    x = math.cos(math.radians(ra_eq))*math.cos(math.radians(dec_eq))
    y = math.sin(math.radians(ra_eq))*math.cos(math.radians(dec_eq))
    z = math.sin(math.radians(dec_eq))
    return x, y, z

###############################
def zenith_year(dec):
    """
    This is a function plots the zenith at night at a particular declination throughout the year.

    Args:

    Returns:

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE:
    #
    #
    # CALLING SEQUENCE:
    #
    #
    # INPUTS:
    #
    # -------------------------------------------------------------
    #
    ralist = [i*10 for i in range(0,36)]
    eqlist = [SkyCoord(ra, dec, frame='icrs', unit=(u.deg, u.deg)) for ra in ralist]
    gallist = [eq.transform_to("galactic") for eq in eqlist]
    points = [(gal.l.deg,gal.b.deg) for gal in gallist]
    for point in points:
        plt.scatter(point[0], point[1], color="blue")
    plt.title("b vs. l")
    plt.xlabel("l")
    plt.ylabel("b")
    plt.show()
    return gallist, points

###############################
def planet_pos_now(planets):
    """
    This is a function that returns the position of planets in ecliptic coordinates at the current time.

    Args:

    Returns:

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE:
    #
    #
    # CALLING SEQUENCE:
    #
    #
    # INPUTS:
    #
    # -------------------------------------------------------------
    #
    planet_list = [astropy.coordinates.get_body(planet, astropy.time.Time.now()) for planet in planets]
    planet_list_ec = [p.transform_to("heliocentrictrueecliptic") for p in planet_list]
    for i in np.arange(0, len(planets)):
        plt.scatter(planet_list_ec[i].lon, planet_list_ec[i].lat, label=planets[i])
    plt.legend()
    plt.show()
    return planet_list_ec



# MAIN
###############################
###############################
if __name__ == '__main__':
    # MEK set the timer
    time0 = time.time()

    # MEK call coordinate conversion function
    coord = cart_conv("04h14m01s", "02:12:01")
    print("from function", coord, "\n")

    # MEK calculate from equations in class notes in order to compare results
    eq_coord = xyz_eq("04h14m01s", "02:12:01")
    print("from class equatinons", eq_coord, "\n")
    print("differences: x =", coord.x - eq_coord[0], "y =", coord.y - eq_coord[1], "z =", coord.z - eq_coord[2], "\n")

    # MEK establish Galactic center coordinates and convert to equatorial
    gal_cen = SkyCoord(0, 0, frame='galactic', unit=(u.hourangle, u.deg))
    icrs_gal_cen = gal_cen.transform_to("icrs")
    print("Galactic center in equatorial frame is", icrs_gal_cen, "\n")

    # MEK find out what constellation this coordinate is in
    con = gal_cen.get_constellation()
    print("Constellation at Galactic center is", con, "! This is right at the edge of the constellation.\n")

    # MEK plot (l,b) overhead at night throughout the year
    print(zenith_year(36))

    # MEK get locations for Mercury, Venus, Mars
    p = planet_pos_now(["mercury", "venus", "mars"])
    print(p)

    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))