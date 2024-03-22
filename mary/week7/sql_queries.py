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
from astropy.io import fits


# FUNCTIONS
###############################
###############################
def func1():
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

    hdu = fits.open("/Users/marykaldor/ASTR8080/mary/week7/sql_test_query.fits")
    hdr = hdu[1].header
    data = hdu[1].data
    ra = data['ra']
    dec = data['dec']
    g = data['g']
    gbins = np.histogram(g, bins=4)
    print(gbins)

    bin1 = np.where((gbins[1][0]<=g) & (g<=gbins[1][1]))[0]
    bin2 = np.where((gbins[1][1]<=g) & (g<=gbins[1][2]))[0]
    bin3 = np.where((gbins[1][2]<=g) & (g<=gbins[1][3]))[0]
    bin4 = np.where((gbins[1][3]<=g) & (g<=gbins[1][4]))[0]

    plt.scatter(ra[bin1], dec[bin1], color="black", marker="o", s=120, alpha=0.5, label="Brightest")
    plt.scatter(ra[bin2], dec[bin2], color="black", marker="o", s=80, alpha=0.5, label="2nd Brightest")
    plt.scatter(ra[bin3], dec[bin3], color="black", marker="o", s=20, alpha=0.5, label="3rd Brightest")
    plt.scatter(ra[bin4], dec[bin4], color="black", marker="o", s=2, alpha=0.5, label="Dimmest")

    plt.title("Dec vs. RA of SDSS Objects with G Mag Weighted Size")
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.legend(fontsize=8)
    plt.gca().invert_xaxis()

    plt.show()

    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))
