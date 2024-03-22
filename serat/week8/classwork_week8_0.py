

#Import
##############
##############
#############

from astropy.io import fits
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import sdss_sweep_data_index
from astropy.coordinates import search_around_sky
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u


# Part 1 & 2
hdul = fits.open("/astro/astr8020/FIRST/first_08jul16.fits")
data_first = hdul[1].data
hdr = hdul[1].header
#print(hdr)


ra = data_first['RA']
dec = data_first['DEC']

plt.figure(figsize=(10,6))
plt.scatter(ra,dec,s=1)
plt.xlabel("Right Ascnension (degrees)")
plt.ylabel("Declination (degrees)")
plt.savefig("First_plot.png")
plt.show()

# Part 3
#for i in tqdm(range(100)):
#	os.system(f"python sdssDR15query.py {ra[i]} {dec[i]}>>file.txt")
	 
# Part 4 & 5
hdul = fits.open("/astro/astr8020/dr15/eboss/sweeps/dr13_final/datasweep-index-star.fits")
data_sdss = hdul[1].data
hdr_sdss = hdul[1].header

ra_new = data_sdss["RA"]
dec_new = data_sdss["DEC"]

plt.figure(figsize=(10,6))
plt.scatter(ra_new,dec_new,s=1)
plt.xlabel("Right Ascension (degrees)")
plt.ylabel("Declination (degrees)")
plt.savefig("sdss_plot.png")
plt.show()


# Part 6
sweepdir="/astro/astr8020/dr15/eboss/sweeps/dr13_final"
file_name = sdss_sweep_data_index.sdss_sweep_data_index(180, 45, 0.5, sweepdir=sweepdir)
#print("Files to be read:", file_name) 


# Part 7
ra_radio = ra[0:100]
dec_radio = dec[0:100]
swfiles = sdss_sweep_data_index.sdss_sweep_data_index(ra_radio, dec_radio,\
														  0, sweepdir=sweepdir)
objs_struct = [ fits.open(file) for file in swfiles ]
objs        = [obj[1].data for obj in objs_struct ] 
objs = np.hstack(objs)
csweeps = SkyCoord(ra=objs["RA"]*u.degree, dec=objs["DEC"]*u.degree)
cin = SkyCoord(ra=ra_new*u.degree, dec=dec_new*u.degree)



idx1, idx2, _, _ = search_around_sky(csweeps, cin, (2/3600)*u.deg)
ra_match = objs["RA"][idx1]
dec_match = objs["DEC"][idx1]

print(ra_match)
print(dec_match)
#plt.figure(figsize=(10,6))
#plt.scatter(ra_match,dec_match,s=1)
#plt.xlabel("Right Ascension (degrees)")
#plt.ylabel("Declination (degrees)")
#plt.savefig("sweep_plot.png")
#plt.show()













