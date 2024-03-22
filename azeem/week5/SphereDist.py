# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pdb
import time

from astropy import units as u
from astropy.coordinates import SkyCoord

if __name__ == '__main__':
    #task 1
    #AB get SkyCoord objs and convert to cartesian
    c1 = SkyCoord(ra=263.75, dec=-12.9, frame='icrs', unit='deg')
    c2 = SkyCoord(ra='20h24m59.9', dec='+10d06m', frame='icrs')

    dot = np.dot(c1.cartesian.xyz,c2.cartesian.xyz)
    print(dot)

    x1 = c1.cartesian.x
    y1 = c1.cartesian.y
    z1 = c1.cartesian.z

    x2 = c2.cartesian.x
    y2 = c2.cartesian.y
    z2 = c2.cartesian.z

    #AB find magnitude of vectors
    mag1 = np.sqrt(x1**2 + y1**2 + z1**2)
    mag2 = np.sqrt(x2**2 + y2**2 + z2**2)

    #find the angle
    ang = np.arccos(np.divide(dot,mag1*mag2))
    print(ang)

    #check for separation
    sep = c2.separation(c1)
    print(sep.rad)

    #task 2
    #AB populate area of sky with 100 random points
    rabox = (np.random.random(100)*15)+30
    decbox = (np.random.random(100)*4)-2
    rabox2 = (np.random.random(100)*15)+30
    decbox2 = (np.random.random(100)*4)-2
    ##for hw
    print(rabox)
    print(decbox)

    #task 3
    #AB cross match catalogs with search_around_sky
    cat1 = SkyCoord(ra=rabox,dec=decbox,unit='deg')
    cat2 = SkyCoord(ra=rabox2,dec=decbox2, unit='deg')

    id1,id2,d2,d3 = cat2.search_around_sky(cat1,(1/6)*u.deg)

    #plot
    plt.scatter(rabox,decbox)
    plt.scatter(rabox2,decbox2)
    plt.scatter(rabox[id1],decbox[id1],marker="x", color='green')
    plt.scatter(rabox2[id2],decbox2[id2],marker="x",color='green')
    plt.show()
