# M. Rizzo Smith
# v1 3/18/24
# Homework 4
# ASTR 8080


# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
import astropy 
from matplotlib import rc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.io import fits
import sys
import pymangle
from numpy.random import random
import time

sys.path.insert(0, '../../week6/')
import spherical_caps
from spherical_caps import ra_cap, dec_cap, gen_cap

# FUNCTIONS
###############################
###############################

def write_ply(caps, num_poly):
    """
    Function to take in an arbitray number of caps. Assumes the first 4 caps correspond to a lat-lon rectangle bounded region, and the follwoing caps are individual spherical caps.

    Args:
        caps (numpy array) - Array of strings that correspond to caps that we wish to write to a ply file.

    Returns:
        Function does not return any variables but creates a .ply file called 'aal_caps.ply' in the current directory
    """
#-----------------------------------------------------------------------
#+
# PURPOSE:
#   Take in an array of caps to write to a single .ply file
#
# CALLING SEQUENCE:
#   write_ply(caps)
#
# INPUTS:
#   caps - numpy array of caps written as strings
#-
#-----------------------------------------------------------------------
    # MRS Create the ply file, specify lat lon as first 4 caps
    f = open('all_caps.ply', 'w')
    
    # MRS Write PLY file assumes 1 for area
    f.write(f'{num_poly} polygons\n')
    for i in range(0, num_poly):
        f.write(f'polygon {i+1}, 1 caps, 1 weight, 0 pixel, 1 str):\n')
        f.write(f' {caps[i][0]:.5f} {caps[i][1]:.5f} {caps[i][2]:.5f} {caps[i][3]:.5f}\n')
    f.close()

        #for circles in caps:
            #f.write(f' {caps[i][0]:.5f} {caps[i][1]:.5f} {caps[i][2]:.5f} {caps[i][3]:.5f}\n')
    #f.close()

    return


###############################
###############################
if __name__ == '__main__':
    spher_caps = np.array([gen_cap(163, 50, 2), gen_cap(167, 50, 2)])
    write_ply(spher_caps, 2)
    
    hdul = fits.open('../../week8/first_08jul16.fits')
    #hdul = fits.open('/astro/astr8020/FIRST/first_08jul16.fits')
    first_data = hdul[1].data
    
    first_objs = SkyCoord(first_data['ra'], first_data['dec'], unit=(u.deg, u.deg))
    m = pymangle.Mangle('all_caps.ply')
    ind_match = m.contains(first_objs.ra.degree, first_objs.dec.degree)
    tp = np.dtype([('ra', 'f8'), ('dec', 'f8')]) 
    first_data = np.zeros(len(first_objs[ind_match]), dtype=tp)
    first_data['ra'] =first_objs[ind_match].ra.degree
    first_data['dec']=first_objs[ind_match].dec.degree

    hdu = fits.BinTableHDU(first_data)
    hdu.writeto('FIRST_Match.fits', overwrite=True)


    test_hdul = fits.open('FIRST_Match.fits')
    test_data = test_hdul[1].data

    plt.plot(first_objs.ra.degree, first_objs.dec.degree, 'o', color='tomato', ms=1, alpha=0.5)
    plt.plot(first_objs[ind_match].ra.degree, first_objs[ind_match].dec.degree, 'o', color='blue')
    plt.plot(test_data['ra'], test_data['dec'], 'x', color='green')
    plt.show()
