# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os


if __name__ == '__main__':
	#task1
	#Read in .fits file
	SDSS = fits.open('/astro/astr8020/FIRST/first_08jul16.fits')
	hdr = SDSS[1].header
	data = SDSS[1].data
	#plot the ra and decs
	ra = data['ra']
	dec = data['dec']
	plt.plot(ra,dec)
	plt.show()

	#task3
	#grab first 100 points
	radiora = data['ra'][:100]
	radiodec = data['dec'][:100]
	#crossmatch radio to optical in SDSS
	for i in np.arange(0,len(radiora)):
		os.system("python ../../runnoe/week8/sdssDR15query.py" + " " + str(radiora[i]) + " " + str(radiodec[i]) + " >> file.txt")

	#task4
	#open an SDSS dr15 file
	sweep = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-008162-1-star.fits.gz')
	hdrsweep = sweep[1].header
	datasweep = sweep[1].data
	#task6
	#plot arrays of RA and dec on the sky
	#shows where files overlap in ra and dec
	index = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/datasweep-index-star.fits')
	hdrindex = index[1].header
	dataindex = index[1].data
	ra_star = dataindex['ra']
	dec_star = dataindex['dec']
	plt.scatter(ra_star, dec_star)
	plt.xlabel('RA')
	plt.ylabel('DEC')
	plt.show()

	#calculate which files would be read to find objects within a certain radius at a given ra and dec (can be arrays)
	import sys
	sys.path.insert(0,'../../runnoe/week8/')
	from sdss_sweep_data_index import sdss_sweep_data_index
	swfiles = sdss_sweep_data_index(180, 45, 0.5, objtype='star',sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final/')
	print(swfiles)
