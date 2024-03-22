# M. Rizzo Smith
# v1 3/5/24
# Week 9 Lecture 1 Tasks
# ASTR 8080 magnitude system transformations


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

    V = 15.256	
    B_V = 0.873	
    g = V + 0.64*(B_V) - 0.13
    g_Nav = 15.70
    print(g)

    star = SkyCoord('16:35:26 +09:47:53', unit=(u.hourangle, u.degree))
    print(star.ra.deg, star.dec.deg)
		     

    swfiles = sdss_sweep_data_index(float(star.ra.deg), float(star.dec.deg), 2./60., objtype='star', sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final')
    #swfiles = sdss_sweep_data_index(248.85833, 9.79806, 2./60., objtype='star', sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final')
    objs_struct = [ fits.open(file) for file in (swfiles) ]
    objs        = [obj[1].data for obj in objs_struct ]
    objs = np.hstack(objs)
    
    objs_pos = SkyCoord(objs['RA'], objs['DEC'], unit=(u.deg, u.deg))
    
    z_ang_sep = star.separation(objs_pos)
    print(len(z_ang_sep))
    good = np.where(z_ang_sep.arcsecond<2)
    
    star_match = objs[good]
    
    star_mags = 22.5 - 2.5*np.log10(star_match["PSFFLUX"][0])
    print(star_mags) 

