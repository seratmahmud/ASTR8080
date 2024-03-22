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
import os
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
    # LSS getting FIRST data
    file = fits.open('/astro/astr8020/FIRST/first_08jul16.fits')
    data = file[1].data
    hdr = file[1].header
    ras = data['ra'][:100] # FIRST RAs
    decs = data['dec'][:100]  # FIRST DECs
    
    # LSS getting a single sdss fits file
    filesdss = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-004136-2-star.fits.gz')
    datasdss = filesdss[1].data
    
    # LSS getting sweeps
    file_sweeps = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/datasweep-index-gal.fits')
    data_sweeps = file_sweeps[1].data
    ra_sweeps = data_sweeps['ra']
    dec_sweeps = data_sweeps['dec'] 


    # for i in range(ras.shape[0]):
    #     cmdstr = f'python sdssDR15query.py {ras[i]} {decs[i]} >> file.txt'
    #     os.system(cmdstr)

    #filesdss = fits.open('/astro/astr8020/calib0bj-004136-2-star.fits.gz')
    #data = filesdss[1].data
    sweeps = sdss_sweep_data_index(ras, decs, 0.5/60./60., sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final')
    fitsweeps = [fits.open(fi) for fi in sweeps]
    fitobjs = [g[1].data for g in fitsweeps]
    fitobjs = np.hstack(fitobjs)
    ra_sdss = fitobjs['RA']
    dec_sdss = fitobjs['DEC']
    
    
    # LSS making object catalogs
    cat1 = SkyCoord(ras, decs, unit=u.degree, frame='icrs')
    cat2 = SkyCoord(ra_sweeps, dec_sweeps, unit=u.degree, frame='icrs')

    id1, id2, d2, d3 = cat2.search_around_sky(cat1, 1/6*u.degree)

    # LSS plotting
    plt.scatter(ra_sweeps, dec_sweeps, marker='.', label='SWEEPS')
    plt.scatter(ras, decs, label='FIRST')
    plt.scatter(ra_sdss, dec_sdss, marker='.', label='random sdss')
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.legend(loc=2)
    plt.show()

    plt.scatter(ras[id1], decs[id1], marker='.', label='FIRST')
    plt.scatter(ra_sweeps[id2], dec_sweeps[id2], marker='.', label='SDSS objs')
    plt.legend()
    plt.show()

    # /astro/astr8020/dr15/eboss/sweeps/dr13_final/calibObj-004136-2-star.fits.gz
    print('hello world XD')

