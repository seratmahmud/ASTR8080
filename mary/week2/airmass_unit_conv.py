# M. Kaldor
# v1 1/18/2024
# ASTR 8080 Airmass and Unit Conversion work

# to import this from another directory:
# import sys
# sys.path.insert(0, '../week1/')
# import airmass_unit_conv
# from airmass_unit_conv import function

# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt
from matplotlib import rc
import pdb
import time
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time


# FUNCTIONS
###############################
###############################
def coord_conv(ra, dec):
    """
    This is a function that converts a right ascension from hours/minutes/seconds to decimal degrees and
    a declination from degrees/minutes/seconds to decimal degrees.

    Args: ra [hours, minutes, seconds]
        dec [degrees, minutes, seconds]

    Returns: ra [degrees]
        dec [degrees]

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: converting into the same units
    #
    #
    # CALLING SEQUENCE: coord_conv("17h18m00s", "49:04:02")
    #
    #
    # INPUTS: ra in hms, dec in degrees minutes seconds
    #
    #
    # -------------------------------------------------------------
    c = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg))
    return c


# MAIN
###############################
###############################
if __name__ == '__main__':
    # MEK set the timer
    time0 = time.time()

    # MEK call coordinate conversion function
    coord = coord_conv("17h18m03s", "49:04:02")
    print("from function", coord)
    ra_func = coord.ra.deg
    dec_func = coord.dec.deg

    # MEK check conversion function against equations from class notes
    ra_eq = 15*(17+(18/60)+(3/3600))
    dec_eq = 49+(4/60)+(2/3600)
    print("my calculations", ra_eq, dec_eq)
    print("difference between", ra_func-ra_eq, ",", dec_func-dec_eq)

    # MEK look at all the options for info you can pull out of SkyCoord object, its RA, and its DEC
    print(dir(coord))
    print(dir(coord.ra))
    print(dir(coord.dec))

    #MEK use time.now() to get today's MJD and JD, check against class notes
    t = Time.now()
    MJD_now = t.mjd
    JD_now = t.jd
    print(MJD_now)
    print(JD_now)
    print("difference between JD and MJD should be 2400000.5 and is", JD_now - MJD_now)


    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))