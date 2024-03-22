# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


if __name__ == '__main__':
    #task2
    #Read in .fits file from SQL
    SDSS = fits.open('sql_test.fits')
    hdr = SDSS[1].header
    data = SDSS[1].data
    #plot the ra and decs
    ra = data['ra']
    dec = data['dec']
    plt.scatter(ra,dec)
    plt.xlabel('RA [degrees]')
    plt.ylabel('dec [degrees]')
    plt.gca().invert_xaxis() #invert for ra
    plt.show()

    #task3
    #plot using bins to scale brightness by point size
    gmag = data['g']
    hist, bin_edges = np.histogram(gmag, bins=4)
    for i in range(0,len(bin_edges)-1):
        bins = (np.where((bin_edges[i] <= gmag) & (gmag <= bin_edges[i+1]))[0])
        plt.scatter(ra[bins],dec[bins],s=200/((i+1)*2), color = 'black', alpha=0.7)
    plt.xlabel('RA [degrees]')
    plt.ylabel('dec [degrees]')
    plt.gca().invert_xaxis() #invert for ra
    plt.show()
