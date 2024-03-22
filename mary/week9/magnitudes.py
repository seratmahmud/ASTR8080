# M. Kaldor
# v1 3/5/2024
# ASTR 8080 magnitudes

# to import this from another directory:
# import sys
# sys.path.insert(0, '../week9/')
# import magnitudes
# from magnitudes import function

# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib
#matplotlib.use("MacOSX")
import matplotlib.pyplot as plt
from matplotlib import rc
import time
import pdb
import math
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
#import sfdmap
#import extinction
from numpy.random import random
import astropy.io.fits as fits
import sys
sys.path.insert(0, '../../runnoe/week8')
from sdss_sweep_data_index import sdss_sweep_data_index


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

    v = 15.256
    bv = 0.873
    g = v+0.60*(bv)-0.12
    print("g =", g, "SDSS says g = 15.70")

    ra = 248.85833
    dec = 9.79806

    swfiles = sdss_sweep_data_index(ra, dec, 5/60/60, sweepdir="/astro/astr8020/dr15/eboss/sweeps/dr13_final")

    hdus = [fits.open(swfile) for swfile in swfiles]
    objs = [hdu[1].data for hdu in hdus]
    objs = np.hstack(objs)

    primaryflag = 2**8
    id = np.where((objs["RESOLVE_STATUS"] & primaryflag) != 0)
    objs_match = objs[id]    

    star = SkyCoord(ra, dec, unit=(u.degree, u.degree))
    coords = SkyCoord(objs["RA"], objs["DEC"], unit=(u.degree, u.degree))
    sep = star.separation(coords)
    small_sep = np.where(sep<2/3600*u.degree)
    sdssmatch = objs[small_sep]
    flux = sdssmatch["PSFFLUX"][0]
    print("Flux", flux)
    mags = 22.5-2.5*np.log10(flux)
    print("Mags", mags)
    print("u =", mags[0], "SDSS says u = 17.28")
    print("g =", mags[1], "SDSS says g = 15.70")
    print("r =", mags[2], "SDSS says r = 15.18")
    print("i =", mags[3], "SDSS says i = 14.71")
    print("z =", mags[4], "SDSS says z = 14.55")

    ra_faint = 248.90813
    dec_faint = 9.81937

    swfiles = sdss_sweep_data_index(ra_faint, dec_faint, 5/60/60, sweepdir="/astro/astr8020/dr15/eboss/sweeps/dr13_final")

    hdus = [fits.open(swfile) for swfile in swfiles]
    objs = [hdu[1].data for hdu in hdus]
    objs = np.hstack(objs)

    primaryflag = 2**8
    id = np.where((objs["RESOLVE_STATUS"] & primaryflag) != 0)
    objs_match = objs[id]

    star_faint = SkyCoord(ra_faint, dec_faint, unit=(u.degree, u.degree))
    sep_faint = star_faint.separation(coords)
    small_sep_faint = np.where(sep_faint<2/3600*u.degree)
    sdssmatch = objs[small_sep_faint]
    flux = sdssmatch["PSFFLUX"][0]
    print("Flux", flux)
    mags = 22.5-2.5*np.log10(flux)
    print("Mags", mags)
    print("u =", mags[0], "SDSS says u = 24.49")
    print("g =", mags[1], "SDSS says g = 25.12")
    print("r =", mags[2], "SDSS says r = 23.09")
    print("i =", mags[3], "SDSS says i = 20.66")
    print("z =", mags[4], "SDSS says z = 19.38")


    
    # MEK time
    time1 = time.time()
    print(" Time [s]: {0:.3f}".format(time1 - time0))
