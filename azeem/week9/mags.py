# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import os

c2 = SkyCoord(ra='16h35m26', dec='+09d47m53', frame='icrs')
#sep = c2.separation(c1)
#print(sep.rad)


V = 15.256
B_V = 0.873
g = V + 0.64*(B_V) - 0.13
print("g mag from UBVRI to ugriz: ",g) #UBVRI mag, ugriz = 15.70 in g band

#sweeping

import sys
sys.path.insert(0,'../../runnoe/week8/')
from sdss_sweep_data_index import sdss_sweep_data_index
#swfiles = sdss_sweep_data_index(248.85833, 9.79806, 0.36, objtype='star',sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final/')


#hdus = [fits.open(swfile) for swfile in swfiles]
#objs = [hdu[1].data for hdu in hdus]
#objs = np.hstack(objs)
#ra = objs["RA"]
#dec = objs["DEC"]
#obj_id = objs["ID"]
#flux = objs["PSFFLUX"]

def getmags(ra_obj, dec_obj):
	swfiles = sdss_sweep_data_index(ra_obj, dec_obj, 0.36, objtype='star',sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final/')
	hdus = [fits.open(swfile) for swfile in swfiles]
	objs = [hdu[1].data for hdu in hdus]
	objs = np.hstack(objs)
	ra = objs["RA"]
	dec = objs["DEC"]
	obj_id = objs["ID"]
	flux = objs["PSFFLUX"]
	c1 = SkyCoord(ra=ra, dec=dec, frame='icrs', unit='deg')
	c2 = SkyCoord(ra=ra_obj, dec=dec_obj, frame='icrs', unit='deg')
	sep = c2.separation(c1)
	good = np.where(sep.arcsecond<2)
	sep = sep[good]

	fluxgood =flux[good]
	print('PSFFLUX',fluxgood)
	m = 22.5-2.5*np.log10(fluxgood)
	print("sweep mags: ",m)
getmags(248.85833, 9.79806)
mags_sdss = [17.28, 15.70, 15.18, 14.71, 14.55]
print("SDSS mags: ", mags_sdss)
#lets try for a faint STAR
getmags(248.91054, 9.76041)
