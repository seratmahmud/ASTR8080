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
def ra_cap(ra):
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
    # MEK you add 6 because for RA caps, the center of the cap is actually 90 degrees offset from where the edge of the
    # cap is - think about the bobble and the brim of the hat! :)
    ra_pt = ra+6
    dec_pt = 0
    pt = SkyCoord(ra_pt, dec_pt, frame="icrs", unit=(u.hourangle, u.deg))
    v = [float(pt.cartesian.x), float(pt.cartesian.y), float(pt.cartesian.z), 1]

    return(v)

###############################
def dec_cap(dec):
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
    ra_pt = 0
    dec_pt = 90
    pt = SkyCoord(ra_pt, dec_pt, frame="icrs", unit=(u.hourangle, u.deg))
    v = [float(pt.cartesian.x), float(pt.cartesian.y), float(pt.cartesian.z), 1-np.sin(np.radians(dec))]

    return (v)

###############################
def gen_cap(ra, dec, r):
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
    pt = SkyCoord(ra, dec, frame="icrs", unit=(u.deg, u.deg))
    v = [float(pt.cartesian.x), float(pt.cartesian.y), float(pt.cartesian.z), 1-np.cos(np.radians(r))]

    return (v)


# MAIN
###############################
###############################
if __name__ == '__main__':
    # MEK set the timer
    time0 = time.time()

    print(ra_cap(5))
    print(dec_cap(36))
    print(gen_cap(5, 36, 1))

    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))
