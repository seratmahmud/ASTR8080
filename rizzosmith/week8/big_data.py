# M. Rizzo Smith
# v1 1/25/24
# Week 8 Lecture 1 Tasks
# ASTR 8080 Survery Matching and Big Data


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
import os

sys.path.insert(0, '../../runnoe/week8/')
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
    #hdul = fits.open('first_08jul16.fits')
    hdul = fits.open('/astro/astr8020/FIRST/first_08jul16.fits')
    data = hdul[1].data

    num_objs = 100
    ra_100 = data['ra'][:num_objs]
    dec_100 = data['dec'][:num_objs]

    #for i in range(0, num_objs):
        #os.system(f"python ../../runnoe/week8/sdssDR15query.py {ra_100[i]} {dec_100[i]}>>output.txt")
    curr_dir = '/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/'
    gal_fits = fits.open(curr_dir+'calibObj-000094-3-gal.fits.gz')
    gal_data = gal_fits[1].data
    
    sweep_ind = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/datasweep-index-gal.fits')
    sweep_dat = sweep_ind[1].data
    
    plt.plot(sweep_dat['ra'], sweep_dat['dec'], 'o', color='lightblue')
    plt.plot(data['ra'], data['dec'], 'o', color='tomato')
    plt.plot(gal_data['ra'], gal_data['dec'], 'o', color='lightgreen')
    plt.show()

    swfiles = sdss_sweep_data_index(180, 45, 0.5, objtype='gal',sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final')
    print('Here', swfiles) 
    match_first = sdss_sweep_data_index(ra_100, dec_100, 2/60/60, objtype='star', sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final')
    objs_struct = [ fits.open(file) for file in (match_first) ]
    objs        = [obj[1].data for obj in objs_struct ]
    objs = np.hstack(objs)
    
    # MRS Pull only primary objects
    primaryflag = 2**8
    w = np.where( (objs["RESOLVE_STATUS"] & primaryflag) != 0)
    objs = objs[w]
    
    csweeps = SkyCoord(ra=objs["RA"]*u.degree, dec=objs["DEC"]*u.degree)
    cat1 = SkyCoord(ra_100, dec_100, unit=(u.deg, u.deg))
    ind1, ind2, d2d, d3d = csweeps.search_around_sky(cat1, 2/60/60*u.degree)

    plt.plot(ra_100, dec_100, 'o', color='tomato')
    plt.plot(ra_100[ind1], dec_100[ind1], 'o', color='blue')
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.show()
    
