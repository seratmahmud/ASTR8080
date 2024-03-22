# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../../runnoe/week8/')
from sdss_sweep_data_index import sdss_sweep_data_index
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
import os
#import sys
#sys.path.insert(0, '.')
#from mags import getmags

stargal = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-008162-1-stargal-primary.fits.gz')

wise = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-008162-1-wise-stargal-primary.fits.gz')

hdrstargal = stargal[1].header
datastargal = stargal[1].data
stargal_ra = datastargal["ra"]

hdrwise = wise[1].header
datawise = wise[1].data
wise_ra = datawise["ra"]
print(len(stargal_ra))
print(len(wise_ra)) #THEY MATCH!

#get PSFFLUX of object
swfiles = sdss_sweep_data_index(143.209, 36.701, 0.36, objtype='star',sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final/')
for i in range(0, len(swfiles)):
	swfiles[i] = swfiles[i].replace("star","stargal-primary")
hdus = [fits.open(swfile) for swfile in swfiles]
objs = [hdu[1].data for hdu in hdus]
objs = np.hstack(objs)
ra = objs["RA"]
dec = objs["DEC"]
obj_id = objs["ID"]
flux = objs["PSFFLUX"]
c1 = SkyCoord(ra=ra, dec=dec, frame='icrs', unit='deg')
#c2 = SkyCoord(ra='16h35m26', dec='+09d47m53', frame='icrs')
c2 = SkyCoord(ra=143.209, dec=36.701, frame='icrs', unit='deg')
sep = c2.separation(c1)
good = np.where(sep.arcsecond<2)
sep = sep[good]
fluxgood =flux[good]
mags_sdss = [17.28, 15.70, 15.18, 14.71, 14.55]
goodobj = obj_id[good]
print('PSFFLUX: ',fluxgood)

#now for WISE object
for i in range(0, len(swfiles)):
        swfiles[i] = swfiles[i].replace("star","wise-star")
wisehdus = [fits.open(swfile) for swfile in swfiles]
wiseobjs = [hdu[1].data for hdu in wisehdus]
wiseobjs = np.hstack(wiseobjs)
wise_star = wiseobjs[good]
wiseflux = [wise_star["W1_NANOMAGGIES"][0],wise_star["W2_NANOMAGGIES"][0]]
print("W1, W2: ",wiseflux)
#getmags(143.209, 36.701)

