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
import sfdmap
import extinction
from numpy.random import random


# FUNCTIONS
###############################
###############################
def func1(x):
    """
    This is a function that

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

###############################
def func2():
    """
    This is a function

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


# MAIN
###############################
###############################
if __name__ == '__main__':
    # MEK set the timer
    time0 = time.time()

    c1 = astropy.coordinates.SkyCoord(263.75, -12.9, unit=(u.deg, u.deg), frame="icrs")
    c2 = astropy.coordinates.SkyCoord("20h24m59.9s", "10:06:00", unit=(u.hourangle, u.deg), frame="icrs")
    cart1 = c1.cartesian.xyz
    cart2 = c2.cartesian.xyz
    dot = cart1[0]*cart2[0]+cart1[1]*cart2[1]+cart1[2]*cart2[2]
    mag_cart1 = np.sqrt(cart1[0]**2+cart1[1]**2+cart1[2]**2)
    mag_cart2 = np.sqrt(cart2[0]**2+cart2[1]**2+cart2[2]**2)
    angle = np.arccos(dot/(mag_cart1*mag_cart2))
    print("dot product angle=", angle)

    sep = c1.separation(c2)
    print("separation angle=", sep.rad, "rad")

    # Multiply by range, shift to minimum value from zero
    ra1 = 15*random(100)+30
    dec1 = 4*random(100)-2
    cat1 = SkyCoord(ra1, dec1, unit=(u.deg, u.deg))
    ra2 = 15 * random(100) + 30
    dec2 = 4 * random(100) - 2
    cat2 = SkyCoord(ra2, dec2, unit=(u.deg, u.deg))
    plt.scatter(ra1, dec1, color="magenta", marker="x")
    plt.scatter(ra2, dec2, color="cyan", marker="o")
    plt.title("DEC vs. RA")
    plt.xlabel("RA")
    plt.ylabel("DEC")

    id1, id2, d2, d3 = cat2.search_around_sky(cat1, 0.166*u.deg)
    match1 = cat1[id1]
    match2 = cat2[id2]
    plt.scatter(match1.ra, match1.dec, color="orange", marker="x")
    plt.scatter(match2.ra, match2.dec, color="orange", marker="o")
    plt.show()

    ra = np.append(ra1, ra2)
    dec = np.append(dec1, dec2)

    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))
