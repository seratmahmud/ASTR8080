# S. Saad
# ASTR 8080 week3, classwork0

# IMPORT BLOCK
###############################
###############################
from numpy.random import random
import numpy as np
import matplotlib.pyplot as plt


from astropy.coordinates import SkyCoord
import sfdmap
from astropy.io import fits
from astropy import wcs
from extinction import fitzpatrick99


# FUNCTIONS
###############################
###############################


def main():
    # Defining RA and Dec using np.random
    ra = 2 * np.pi * (random(10000) - 0.5) 
    dec = np.arcsin(1. - random(10000) * 2.)
    
    # Plotting in 2D plane
    plt.scatter(ra, dec, alpha=0.5, s=1)
    plt.xlabel('Right Ascension in radians')
    plt.ylabel('Declination in radians')
    plt.grid()
    plt.show()
    
    # Aitoff plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="aitoff") 
    ax.scatter(ra, dec, alpha=0.5, s=1)
    ax.grid(color='blue', linestyle='solid', linewidth=0.5)
    fig.show()
    
    # Lambert plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="lambert") 
    ax.scatter(ra, dec, alpha=0.5, s=1)
    ax.grid(color='blue', linestyle='solid', linewidth=0.5)
    fig.show()
    
    # Taking RA and Dec value to form a meshgrid
    ra_center_1, dec_center_1 = 0, 0
    ra_bins_1, dec_bins_1 = 1, 1
    size = 100
    ra_range_1 = np.linspace(ra_center_1 - 180 * ra_bins_1, ra_center_1 + 180 * ra_bins_1, size)
    dec_range_1 = np.linspace(dec_center_1 - 90 * dec_bins_1, dec_center_1 + 90 * dec_bins_1, size)
    
    ra_grid, dec_grid = np.meshgrid(ra_range_1, dec_range_1)
    sc = SkyCoord(ra_grid, dec_grid, unit='deg')
    sc.galactic
    
    # Getting the redenning
    dustdir = 'sfddata-master/'
    m = sfdmap.SFDMap(dustdir, scaling=1)
    ebv = m.ebv(sc.ra.value,sc.dec.value, frame='galactic')
    
    ebv_max = np.max(ebv)
    
    clevels = np.array([0.999, 0.99, 0.9])
    clevels = (1-clevels)*ebv_max
    
    # Contour plotting to see the galactic plane in 2D
    plt.contour(sc.ra.value, sc.dec.value, ebv, levels=clevels)
    plt.xlabel("Galactic Longitude (l)")
    plt.ylabel("Galactic Latitude (b)")
    plt.show()
    
    # Using wcs to plot it in aitoff
    w = wcs.WCS(naxis=2) 
    w.wcs.cdelt=[1,1] 
    w.wcs.crval=[0.5,-89.5] 
    w.wcs.crpix=[1,1] 
    w.wcs.ctype = ["RA---AIT", "DEC--AIT"] 
    xmap,ymap = w.wcs_world2pix(ra_grid, dec_grid, 0)
    lev = np.arange(50)*0.03 
    plt.contourf(xmap,ymap,ebv,levels=lev) 
    plt.show()
    
    

    
# MAIN
###############################
###############################
if __name__ == '__main__':
    main()