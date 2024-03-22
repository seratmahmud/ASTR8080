# M. Kaldor
# v1 3/7/2024
# ASTR 8080 forced photometry

# to import this from another directory:
# import sys
# sys.path.insert(0, '../week9/')
# import forced_photometry
# from forced_photometry import function

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
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
#import sfdmap
#import extinction
from numpy.random import random
from astropy.io import fits
import os
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

    dr15path = "/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-008162-6-stargal-primary.fits.gz"
    wisepath = "/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-008162-6-wise-stargal-primary.fits.gz"
    dr15hdu = fits.open(dr15path)
    wisehdu = fits.open(wisepath)
    dr15hdr = dr15hdu[0].header
    wisehdr = wisehdu[0].header
    dr15data = dr15hdu[1].data
    wisedata = wisehdu[1].data

    print("Number of objects in DR15 is", np.shape(dr15data)[0])
    print("Number of objects in WISE is", np.shape(wisedata)[0])

    ra = 143.209
    dec = 36.701
    swfiles = sdss_sweep_data_index(ra, dec, 2/60/60, sweepdir="/astro/astr8020/dr15/eboss/sweeps/dr13_final")
    swdr15 = [file.replace("star", "stargal-primary") for file in swfiles]
    swwise = [file.replace("star", "wise-stargal-primary") for file in swfiles]

    dr15hdus = [fits.open(file) for file in swdr15]
    dr15objs = [hdu[1].data for hdu in dr15hdus]
    dr15objs = np.hstack(dr15objs)

    wisehdus = [fits.open(file) for file in swwise]
    wiseobjs = [hdu[1].data for hdu in wisehdus]
    wiseobjs = np.hstack(wiseobjs)

    primaryflag = 2**8
    id = np.where(((dr15objs["RESOLVE_STATUS"] & primaryflag) != 0) & dr15objs["OBJC_TYPE"]==6)
    objs_match = dr15objs[id]

    ras = dr15objs["RA"]
    decs = dr15objs["DEC"]
    coords = SkyCoord(ras, decs, unit=(u.degree, u.degree))

    star = SkyCoord(ra, dec, unit=(u.degree, u.degree))

    sep = star.separation(coords)
    small_sep = np.where(sep<2/3600*u.degree)
    dr15match = dr15objs[small_sep][0]
    wisematch = wiseobjs[small_sep]
    dr15flux = dr15match["PSFFLUX"]
    wiseflux = [wisematch["W1_NANOMAGGIES"][0], wisematch["W2_NANOMAGGIES"][0]]

    print("DR15 flux", dr15flux)
    print("WISE flux", wiseflux)

    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))
