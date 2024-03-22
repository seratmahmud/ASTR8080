# M. Rizzo Smith
# v1 3//7/24
# Week 9 Lecture 2 Tasks
# ASTR 8080 forced photometry


# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
import astropy 
from matplotlib import rc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic
from astropy.coordinates import get_body
from astropy.time import Time
from astropy.io import fits
import sys
sys.path.insert(0, '.')
import sdss_sweep_data_index
from sdss_sweep_data_index import sdss_sweep_data_index

# FUNCTIONS
###############################
###############################

# None for now

# MAIN
###############################
###############################
if __name__ == '__main__':
    hdul = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-008162-3-stargal-primary.fits.gz')
    sdss_dat = hdul[1].data

    hdul2 = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-008162-3-wise-stargal-primary.fits.gz')
    wise_data = hdul2[1].data
    print(len(sdss_dat), len(wise_data))

    test_star = SkyCoord(143.209, 36.701, unit=(u.deg, u.deg))
    swfiles = sdss_sweep_data_index(float(test_star.ra.deg), float(test_star.dec.deg), 2./60., objtype='star', sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final')
    
    for i in range(0, len(swfiles)):
        swfiles[i] = swfiles[i].replace('star', 'stargal-primary')
    
    objs_struct = [ fits.open(file) for file in (swfiles) ]
    objs        = [obj[1].data for obj in objs_struct ]
    objs = np.hstack(objs)

    objs_pos = SkyCoord(objs['RA'], objs['DEC'], unit=(u.deg, u.deg))
    
    # MRS Calculate sep from star to get index 
    z_ang_sep = test_star.separation(objs_pos)
    good = np.where(z_ang_sep.arcsecond<2)
     
    
    star_match = objs[good]

    star_sdss_flux = star_match["PSFFLUX"][0]
    print(star_sdss_flux)

    wise_file = [ line.replace('star', 'wise-star') for line in swfiles]
    
    wise_struct = [fits.open(file) for file in wise_file]
    wise_objs = [obj[1].data for obj in wise_struct]
    wise_objs = np.hstack(wise_objs)

    wise_star = wise_objs[good]
    star_wise_flux = [wise_star["W1_NANOMAGGIES"][0], wise_star["W2_NANOMAGGIES"][0]]
    print(star_wise_flux)
