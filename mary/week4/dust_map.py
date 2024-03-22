# M. Kaldor
# v1 1/30/2024
# ASTR 8080 Dust Map

# to import this from another directory:
# import sys
# sys.path.insert(0, '../week1/')
# import dust_map
# from dust_map import function

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


# FUNCTIONS
###############################
###############################
def dust_correct(ra, dec):
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
    z = 1.08
    # c = SkyCoord(ra, dec).galactic
    c = SkyCoord(ra, dec, frame="icrs", unit=u.deg)
    m = sfdmap.SFDMap(dustdir, scaling=1)
    # ebv = m.ebv(c.l.value, c.b.value, frame="galactic")
    ebv = m.ebv(c.ra.value, c.dec.value)
    wave = np.array([3543., 4770., 6231, 7625., 9134.])
    A = extinction.fitzpatrick99(wave, 3.1 * ebv)
    return A


###############################
def make_grid(xcenter, ycenter, step, num):
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
    range = step*num
    x = np.arange(xcenter-(range/2), xcenter+(range/2)-1, step)
    y = np.arange(ycenter-(range/2), ycenter+(range/2)-1, step)
    print(xcenter)
    print(np.median(x))
    print(ycenter)
    print(np.median(y))
    print(x, "\n")
    print(y)
    xmap, ymap = np.meshgrid(x, y)
    plt.plot(xmap, ymap, marker=",", color='k', linestyle='none')
    plt.show()


# MAIN
###############################
###############################
if __name__ == '__main__':
    # MEK set the timer
    time0 = time.time()

    dustdir = "/Users/marykaldor/ASTR8080/mary/week4/sfddata/"

    ra1, dec1 = 246.933, 40.795
    ra2, dec2 = 236.562, 2.440
    gri1 = [18.81, 18.74, 18.81]
    gri2 = [19.10, 18.79, 18.72]
    gr1 = gri1[0]-gri1[1]
    gr2 = gri2[0]-gri2[1]
    ri1 = gri1[1]-gri1[2]
    ri2 = gri2[1]-gri2[2]
    plt.scatter((ri1, ri2), (gr1, gr2))
    plt.xlabel("r-i")
    plt.ylabel("g-r")
    plt.title("g-r vs. r-i")
    # plt.show()

    # These quasars have relatively similar colors and they should!

    A1 = dust_correct(ra1, dec1)
    A2 = dust_correct(ra2, dec2)
    z = 1.08
    gri1corr = gri1 - A1[1:4]
    gri2corr = gri2 - A2[1:4]
    gr1corr = gri1corr[0] - gri1corr[1]
    gr2corr = gri2corr[0] - gri2corr[1]
    ri1corr = gri1corr[1] - gri1corr[2]
    ri2corr = gri2corr[1] - gri2corr[2]
    plt.scatter((ri1corr, ri2corr), (gr1corr, gr2corr))
    plt.xlabel("r-i")
    plt.ylabel("g-r")
    plt.title("g-r vs. r-i")
    # plt.show()
    plt.close()

    # Their colors are more similar!

    make_grid(236.6, 2.4, 1, 100)
    make_grid(246.9, 40.8, 1.3, 100)


    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))
