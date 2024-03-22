# M. Rizzo Smith
# v1 2/22/24
# Week 7 Lecture 2 Tasks
# ASTR 8080 SQL Learning tasks


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
# FUNCTIONS
###############################
###############################

# None for now

# MAIN
###############################
###############################
if __name__ == '__main__':
    hdul = fits.open('SQL_query.fits')
    data = hdul[1].data
    binned = np.histogram(data['g'], bins = 4)
    for i in range(0, len(binned[1])-1):
        bins = (np.where((binned[1][i] <= data['g']) & (data['g'] <= binned[1][i+1]))[0])
        plt.scatter(data['ra'][bins], data['dec'][bins], s=200/((i+1)*2), color = 'black', alpha=0.7)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.xlim(np.max(data['ra']) + 0.5/60, np.min(data['ra'])-0.5/60)
    plt.gca().invert_xaxis()
    plt.show()
