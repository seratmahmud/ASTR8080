# L. Schult
# v1 
# ASTR 8080 
# to import this from another directory:
# import sys
# sys.path.insert(0, '../week1/')
# import polycalc
# from polycalc import get_poly_o3

# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import astropy
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import pdb
import time
import sys
sys.path.insert(0, '../../runnoe/week8/')
from sdss_sweep_data_index import sdss_sweep_data_index


# FUNCTIONS
###############################
###############################
def testname(x):
    """
    This is a function that 

    Args:
        x (datatype): A...
    Returns:
        y (datatype): A ...
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   You pass butter
#
# CALLING SEQUENCE:
#   variable = testname(x) 
#
# INPUTS:
#   x
#-
#-------------------------------------------------------------  
    y     = x**2. 
    return y 


# MAIN
###############################
###############################
if __name__ == '__main__':
    sdss_sweep = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-004002-1-stargal-primary.fits.gz')
    wise_sweep = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-004002-1-wise-stargal-primary.fits.gz')
    sdss_data = sdss_sweep[1].data
    wise_data = wise_sweep[1].data
    print(sdss_data.shape)
    print(wise_data.shape)

    PG_skycoord = SkyCoord(143.209, 36.701, unit=u.degree) 
    sweeps = sdss_sweep_data_index(float(PG_skycoord.ra.value), float(PG_skycoord.dec.value), 2./60., objtype='star', sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final')
    sweeps = [s.replace('star', 'stargal-primary') for s in sweeps]
    fitsweeps = [fits.open(fi) for fi in sweeps]
    objsweeps = [g[1].data for g in fitsweeps]
    objsweeps = np.hstack(objsweeps)
    objskycoord = SkyCoord(objsweeps['RA'], objsweeps['DEC'], unit=u.degree)
    objectseps = PG_skycoord.separation(objskycoord)
    objectind = np.where(objectseps.arcsecond < 2 )
    print(objskycoord[objectind]) 
    print('ugriz mags', objsweeps[objectind]['PSFFLUX']) 
    wisesweeps = [f.replace('stargal', 'wise-stargal') for f in sweeps]
    wisefits = [fits.open(fi) for fi in wisesweeps]
    wisedata = [fi[1].data for fi in wisefits]
    wisedata = np.hstack(wisedata)
    wiseobj = wisedata[objectind]
    print('W1 flux:', wiseobj['W1_NANOMAGGIES'])
    print('W2 flux:', wiseobj['W2_NANOMAGGIES'])

    print('hello world XD')

