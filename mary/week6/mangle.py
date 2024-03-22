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
import pymangle

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
from matplotlib.patches import Rectangle
import math
import sys
sys.path.insert(0, '/Users/marykaldor/ASTR8080/mary/week6/spher_cap.py')
import spher_cap
from spher_cap import gen_cap


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

    cap1 = gen_cap(76, 36, 5)
    cap2 = gen_cap(75, 35, 5)
    print(cap1)
    print(cap2)

    minter = pymangle.Mangle("intersection.ply")
    mbothcaps = pymangle.Mangle("bothcaps.ply")
    mflip1 = pymangle.Mangle("intflip.ply")

    ra_rand_int, dec_rand_int = minter.genrand(10000)
    ra_rand_b, dec_rand_b = mbothcaps.genrand(10000)
    ra_rand_flip, dec_rand_flip = mflip1.genrand(10000)

    plt.scatter(ra_rand_b, dec_rand_b, color="red")
    plt.scatter(ra_rand_int, dec_rand_int, color="blue")
    plt.savefig("/Users/marykaldor/ASTR8080/mary/week6/bothcaps.png")
    plt.show()

    plt.scatter(ra_rand_int, dec_rand_int, color="blue")
    plt.scatter(ra_rand_flip, dec_rand_flip, color="yellow")
    plt.savefig("/Users/marykaldor/ASTR8080/mary/week6/flip.png")
    plt.show()

    '''
    f = open("intersection.ply", "a")
    f.write("1 polygons\n")
    f.write("polygon 1 ( 2 caps, 1 weight, 0 pixel ")
    f.write(str(cap1))
    f.close()
    '''

    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))
