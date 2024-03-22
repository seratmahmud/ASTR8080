# J. Runnoe
# 10/7/2020
# This code searches the sweeps within 0.5 deg
# of a test RA/dec position.  It returns the 
# sweep files that it finds.  Then it plots the
# RA dec values in all of the sweeps, with the RA
# dec from the index file overplotted.
# The goal is to show what the index is doing.
# This doesn't query the sweeps so you have to 
# manually download the data first.

# this code assumes you have manually downloaded:
# datasweep-index-star.fits
# and in 301/ 
# calibObj-003647-1-star.fits.gz	calibObj-003647-4-star.fits.gz	calibObj-003704-3-star.fits.gz	calibObj-004335-6-star.fits.gz
# calibObj-003647-2-star.fits.gz	calibObj-003704-1-star.fits.gz	calibObj-003704-4-star.fits.gz	calibObj-004516-6-star.fits.gz
# calibObj-003647-3-star.fits.gz	calibObj-003704-2-star.fits.gz	calibObj-003900-6-star.fits.gz
# from
# https://dr15.sdss.org/sas/dr15/eboss/sweeps/dr13_final/

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from time import time
from astropy.io import fits 
from sdss_sweep_data_index import sdss_sweep_data_index
from matplotlib.pyplot import cm
import matplotlib
import matplotlib.pyplot as plt

if __name__=="__main__":

    # read in the index file
    index_struct = fits.open('datasweep-index-star.fits')
    index     = index_struct[1].data
    index_hdr = index_struct[1].header

    # grab sweeps for a random position
    swfiles = sdss_sweep_data_index(145.285854,34.741254,0.5)

    # get some colors
    clr = cm.plasma(np.linspace(0,1,np.shape(swfiles)[0]))

    # loop over the sweeps file and plot RA,dec 
    # from each one in a different color
    for i in range(np.shape(swfiles)[0]):
        #import pdb;pdb.set_trace()
        f   = fits.open(swfiles[i])
        dat = f[1].data
        ra  = dat['RA']
        dec = dat['DEC']

        plt.plot(ra,dec,color=clr[i],marker=',')
        plt.xlabel('RA')
        plt.ylabel('DEC')

    # now find RA,decs in the index file that matched
    # to our point, just recreate the sdss_sweep_data_index

    cin    = SkyCoord(ra=145.285854*u.degree, dec=34.741254*u.degree) 
    cindex = SkyCoord(ra=index["RA"]*u.degree, dec=index["DEC"]*u.degree)
    sep    = cin.separation(cindex)
    m2     = np.where(sep < (0.5+0.36)*u.deg)[0]

    ra = index['RA'][m2]
    dec = index['DEC'][m2]
    plt.plot(ra,dec,'k+')
    plt.plot(145.285854,34.741254,'cP')
    #plt.xlim(142,150)
    #plt.ylim(33,36)
    plt.tight_layout()
    plt.savefig('./sweep.png', format='png')
    plt.show()
    import pdb;pdb.set_trace()
