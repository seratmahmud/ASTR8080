# M. Rizzo Smith
# v1 1/25/24
# Week 3 Lecture 2 Tasks
# ASTR 8080 coordinate system transformations


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
from numpy.random import random
import healpy


# FUNCTIONS
###############################
###############################

# None for now

# MAIN
###############################
###############################
if __name__ == '__main__':

    ra = 360.* (random(1000000))
    dec = (180./np.pi) * (np.arcsin(1-random(1000000)*2))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='aitoff')
    ax.scatter(ra, dec, s=0.1, color='tomato', alpha=0.1)
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=800)
    ax.grid(color='black', linestyle='solid', linewidth=1.5)
    plt.show()


    phi = ra * (np.pi)/180.
    theta = (np.pi / 2) - (dec * (np.pi/180))
    
    pix = healpy.ang2pix(1, theta, phi)

    # Area of pixel is going to be 4pi steradians divided by 12*Nside**2
    area = (4 * np.pi) / 12.
    area *= (180./np.pi)**2
    hp = healpy.nside2pixarea(1, degrees=True)
    print(hp, area)
    
    points = np.histogram(pix, 12)
    plt.hist(pix, 12)
    plt.show()

    w = np.where(pix == 2)
    plt.plot((phi/np.pi), np.cos(theta), 'o', ms=0.1, color='tomato', alpha=0.2)
    plt.plot((phi[w]/np.pi), np.cos(theta[w]), 'o', ms=0.1, color='saddlebrown', alpha=0.2)
    plt.show()
    
    ra_rad = (ra - 180) * (np.pi / 180) 
    dec_rad = dec * (np.pi / 180)
    print(np.min(ra_rad), np.max(ra_rad))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='aitoff')
    ax.scatter(ra_rad, dec_rad, s=0.1, color='tomato', alpha=0.2)
    ax.scatter(ra_rad[w], dec_rad[w], s=0.1, color='saddlebrown', alpha=0.2)
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=800)
    ax.grid(color='black', linestyle='solid', linewidth=1.5)
    plt.show()

    pix2 = healpy.ang2pix(2, theta, phi) 
    p5 = np.where(pix == 5)
    #p2in5 = np.where(pix2 == (np.where(pix == 5)))
    plt.plot((phi/np.pi), np.cos(theta), 'o', ms=0.1, color='tomato', alpha=0.2)
    plt.plot((phi[p5]/np.pi), np.cos(theta[p5]), 'o', ms=0.1, color='saddlebrown', alpha=0.2)
    plt.show()

