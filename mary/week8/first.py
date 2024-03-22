# M. Kaldor
# v1 2/27/2024
# ASTR 8080 first

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

    hdu = fits.open("/astro/astr8020/FIRST/first_08jul16.fits")
    hdr = hdu[1].header
    data = hdu[1].data
    ra = data['ra']
    dec = data['dec']
    plt.scatter(ra, dec, color="blue", alpha=0.5)
    '''
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.title("FIRST 08 July 2016")
    '''
    #plt.savefig("/Users/marykaldor/ASTR8080/mary/week8/first.png")

    ra100 = ra[:100]
    dec100 = dec[:100]

    '''
    for i in range(0, len(ra100)):
        os.system(f"python ../../runnoe/week8/sdssDR15query.py {ra100[i]} {dec100[i]} >> first_100.txt")
    '''

    hdu = fits.open("/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-000094-4-gal.fits.gz")
    hdr = hdu[1].header
    data = hdu[1].data
    ra = data['ra']
    dec = data['dec']
    plt.scatter(ra, dec, color="black", alpha=0.5)

    hdu = fits.open("/astro/astr8020/dr15/eboss/sweeps/dr13_final/datasweep-index-gal.fits")
    hdr = hdu[1].header
    data = hdu[1].data
    ra = data['ra']
    dec = data['dec']
    plt.scatter(ra, dec, color="yellow", alpha=0.01)
    plt.title("Sweeps and Quasars")
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.show()

    swfiles = sdss_sweep_data_index(ra100, dec100, 5/60/60, sweepdir="/astro/astr8020/dr15/eboss/sweeps/dr13_final")
    objs_struct = [fits.open(file) for file in swfiles]
    objs = [obj[1].data for obj in objs_struct]
    objs = np.hstack(objs)

    primaryflag = 2 ** 8
    w = np.where((objs["RESOLVE_STATUS"] & primaryflag) != 0)
    objs = objs[w]

    cat1 = SkyCoord(ra100, dec100, unit=(u.deg, u.deg))
    csweeps = SkyCoord(ra=objs["RA"]*u.degree, dec=objs["DEC"]*u.degree)
    id1, id2, d2, d3 = csweeps.search_around_sky(cat1, 5/60/60 * u.deg)
    match1 = cat1[id1]
    match2 = csweeps[id2]
    plt.scatter(ra100, dec100, color="black")
    plt.scatter(match2.ra, match2.dec, color="blue", marker="o", label="Sweeps")
    plt.scatter(match1.ra, match1.dec, color="magenta", marker="x", label="First 100")
    plt.title("Quasars Sweep Match")
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.legend()
    plt.show()

    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))
