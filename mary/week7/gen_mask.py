# M. Kaldor
# v1 2/20/2024
# ASTR 8080 gen_mask

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
import pymangle
import sys
sys.path.insert(0, '/Users/marykaldor/ASTR8080/mary/week6/')
from spher_cap import ra_cap
from spher_cap import dec_cap
sys.path.insert(0, '/Users/marykaldor/ASTR8080/mary/hw/')
from hw2 import square_area
from hw2 import random_pop

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

    '''
    mask = pymangle.Mangle(“file.ply”)
    ra_rand, dec_rand = mask.genrand(10000)
    '''

    ra5 = ra_cap(5)
    ra6 = ra_cap(6)
    dec30 = dec_cap(30)
    dec40 = dec_cap(40)

    area_deg = square_area(5*15, 6*15, 30, 40)[0]
    area_str = area_deg*((np.pi/180)**2)
    print(area_str, "Sr")

    f = open("mask.ply", "w")
    f.write("2 polygons\n")
    f.write(" polygon 1 ( 4 caps, 0.9 weight, 0 pixel, "+str(area_str)+" str):\n")
    f.write("  "+str(ra5[0])+" "+str(ra5[1])+" "+str(ra5[2])+" "+str(ra5[3])+"\n")
    f.write("  "+str(ra6[0])+" "+str(ra6[1])+" "+str(ra6[2])+" "+str(ra6[3]*-1)+"\n")
    f.write("  "+str(dec30[0])+" "+str(dec30[1])+" "+str(dec30[2])+" "+str(dec30[3])+"\n")
    f.write("  "+str(dec40[0])+" "+str(dec40[1])+" "+str(dec40[2])+" "+str(dec40[3]*-1)+"\n")
    f.close()

    ra10 = ra_cap(10)
    ra12 = ra_cap(12)
    dec60 = dec_cap(60)
    dec70 = dec_cap(70)

    area_deg_2 = square_area(10 * 15, 12 * 15, 60, 70)[0]
    area_str_2  = area_deg_2 * ((np.pi / 180) ** 2)
    print(area_str_2, "Sr")

    f = open("mask.ply", "a")
    f.write(" polygon 2 ( 4 caps, 0.2 weight, 0 pixel, " + str(area_str_2) + " str):\n")
    f.write("  " + str(ra10[0]) + " " + str(ra10[1]) + " " + str(ra10[2]) + " " + str(ra10[3]) + "\n")
    f.write("  " + str(ra12[0]) + " " + str(ra12[1]) + " " + str(ra12[2]) + " " + str(ra12[3] * -1) + "\n")
    f.write("  " + str(dec60[0]) + " " + str(dec60[1]) + " " + str(dec60[2]) + " " + str(dec60[3]) + "\n")
    f.write("  " + str(dec70[0]) + " " + str(dec70[1]) + " " + str(dec70[2]) + " " + str(dec70[3] * -1) + "\n")
    f.close()

    mask = pymangle.Mangle("mask.ply")
    print(mask)

    fullsky = random_pop(0,360,-90,90,1000000)
    good = mask.contains(fullsky[0], fullsky[1])
    inmask = np.where(good==True)
    print(inmask)

    plt.scatter(fullsky[0], fullsky[1], color="black")
    plt.scatter(fullsky[0][inmask], fullsky[1][inmask], color="magenta")
    plt.show()

    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))
