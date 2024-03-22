# M. Kaldor
# v1 2/1/2024
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

    ra = 2 * np.pi * (random(10000) - 0.5)
    dec = np.arcsin(1. - random(10000) * 2.)
    plt.scatter(ra, dec, marker=",", color="black", s=0.3, alpha=0.5)
    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.title("DEC vs. RA")
    plt.show()
    # There are fewer points at the poles and more at the equator.

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="aitoff")
    ax.scatter(ra, dec, marker="o", color="red", s=0.3, alpha=0.8)
    ax.grid(color="blue", linestyle="dashed", linewidth=2.5)
    xlab = ["14h", "16h", "18h", "20h", "22h", "0h", "2h", "4h", "6h", "8h", "10h"]
    ax.set_xticklabels(xlab, weight=800)
    plt.title("Aitoff Projection")
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="lambert")
    ax.scatter(ra, dec, marker="o", color="red", s=0.3, alpha=0.8)
    ax.grid(color="blue", linestyle="dashed", linewidth=2.5)
    xlab = ["14h", "16h", "18h", "20h", "22h", "0h", "2h", "4h", "6h", "8h", "10h"]
    ax.set_xticklabels(xlab, weight=800)
    plt.title("Lambert Projection")
    plt.show()

    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))
