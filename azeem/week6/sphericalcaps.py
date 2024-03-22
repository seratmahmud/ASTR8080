# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pdb
import time
import math

from astropy import units as u
from astropy.coordinates import SkyCoord

#task1
def ra_bound(ra,dec):
    #AB all in degrees
    #add 90 to ra
    #make SkyCoord object
    c = SkyCoord(ra=ra+90, dec=dec, unit='deg')

    #make 4-array
    four_array = np.array([c.cartesian.x, c.cartesian.y, c.cartesian.z, 1])
    return(four_array)

#task2
def dec_bound(ra1,dec1):
    #AB convert dec to radians
    h = np.radians(dec1)
    #take sin
    h = np.sin(h)

    #AB make SkyCoord object
    c1 = SkyCoord(ra=ra1, dec=90, unit='deg')
    #print(dec1)
    #make 4-array
    four_array1 = np.array([c1.cartesian.x, c1.cartesian.y,c1.cartesian.z, 1-h])
    #print(four_array1)
    return(four_array1)
#task3
def circle_cap(ra,dec,theta):
    #AB for more general case
    #convert theta to rads and take cosine
    h1 = np.radians(theta)
    h1 = np.cos(h1)
    #make SkyCoord object
    c2 = SkyCoord(ra=ra, dec=dec, unit='deg')
    #make 4-array
    four_array2 = np.array([c2.cartesian.x, c2.cartesian.y, c2.cartesian.z, 1-h1])
    return(four_array2)


if __name__ == '__main__':
    ra = 5*15 #converts from hms to dec
    dec = 0 #in degrees
    ra_bound(ra, dec)

    ra2 = 0
    dec2 = 30
    dec_bound(ra2,dec2)

    ra3 = 5*15
    dec3 = 36
    theta=1
    circle_cap(ra3,dec3,theta)


    # #for mangle
    # #cap 1
    # theta3 = 5
    # h3 = np.radians(theta3)
    # h3 = np.cos(h3)
    # c3 = SkyCoord(ra=76, dec=36, unit='deg')
    # four_array3 = np.array([c3.cartesian.x, c3.cartesian.y, c3.cartesian.z, 1-h3])
    # print(four_array3)
    #
    # #cap 2
    # theta4 = 5
    # h4 = np.radians(theta4)
    # h4 = np.cos(h4)
    # c4 = SkyCoord(ra=75, dec=35, unit='deg')
    # four_array4 = np.array([c4.cartesian.x, c4.cartesian.y, c4.cartesian.z, 1-h4])
    # print(four_array4)
    #
    # #flip cap
    # theta5 = 5
    # h5 = np.radians(theta5)
    # h5 = np.cos(h5)
    # c5 = SkyCoord(ra=76, dec=36, unit='deg')
    # four_array5 = np.array([c5.cartesian.x, c5.cartesian.y, c5.cartesian.z, -1*(1-h5)])
    # print(four_array5)
