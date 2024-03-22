# L. Schult
# v1 Tues 6 Feb 2024
# ASTR 8080 
# to import this from another directory:
# import sys
# sys.path.insert(0, '../week1/')
# import polycalc
# from polycalc import get_poly_o3

# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import random
from matplotlib import rc
import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
import pdb
import time
import healpy as hp

# FUNCTIONS
###############################
###############################
def testname(x):
    """
    This is a function that 

    Args:
        x (datatype): A...
    Returns:
        y (datatype): A ...
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   You pass butter
#
# CALLING SEQUENCE:
#   variable = testname(x) 
#
# INPUTS:
#   x
#-
#-------------------------------------------------------------  
    y     = x**2. 
    return y 

def randomglobe(rarange=[0, 2*np.pi], decrange=[-np.pi/2, np.pi/2], 
                projection='xy', n=10000, plot=True):
    """
    This is a function that generates random points and plots them.
    projections available are xy, aitoff, and lambert.

    Args:
        projection (string): the projection the map will be in. default is xy
        n (int): the number of random points to generate. default is 10000
    Returns:
        y (datatype): A ...
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   plot random points on different maps
#
# CALLING SEQUENCE:
#   randomglobe(projection='aitoff') 
#
# INPUTS:
#   projection, n
#-
#-------------------------------------------------------------  
    ra = ((rarange[1]-rarange[0])*(random(n)))+rarange[0]
    dec = np.arcsin(((np.sin(decrange[1])-np.sin(decrange[0]))*random(n))+np.sin(decrange[0]))
    if plot:
        if projection=='xy':
            plt.scatter(ra, dec, marker=',', s=1)
            plt.xlabel('ra')
            plt.ylabel('dec')
            plt.show()
        elif projection == 'aitoff' or 'lambert':
            fig = plt.figure()
            xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
            ax = fig.add_subplot(111, projection=projection)
            ax.set_xticklabels(xlab, weight=800)
            ax.grid(color='b', linestyle='dashed', linewidth=3)
            ax.scatter(ra-np.pi, dec, marker=',', s=1)
            #ax.set_xticks()
            plt.show()
    
    return ra, dec

# MAIN
###############################
###############################
if __name__ == '__main__':
    c1 = SkyCoord(263.75, -12.9, unit='deg', frame='icrs')
    c2 = SkyCoord('20h24m59.9s','+10d6m0s', frame='icrs')

    c1cart = [c1.cartesian.x, c1.cartesian.y, c1.cartesian.z]
    c2cart = [c2.cartesian.x, c2.cartesian.y, c2.cartesian.z]

    # LSS getting angle from dotprod
    dotprod = np.dot(c1cart, c2cart)
    mag = lambda x: np.sqrt(x[0]**2+x[1]**2+x[2]**2)
    mag1 = mag(c1cart)
    mag2 = mag(c2cart)
    angle = np.arccos(dotprod / (mag1*mag2))
    print(f'separation from dotprod {angle}')

    # LSS separation angle using astropy separation
    print(f'separation angle from separation: {np.deg2rad(c1.separation(c2))}')
    
    # LSS generating small area random points
    smallrandra1, smallranddec1 = randomglobe(rarange=[np.deg2rad(30), np.deg2rad(45)], 
                decrange=[np.deg2rad(-2), np.deg2rad(2)], n=100, plot=False)
    smallrandra2, smallranddec2 = randomglobe(rarange=[np.deg2rad(30), np.deg2rad(45)], 
                decrange=[np.deg2rad(-2), np.deg2rad(2)], n=100, plot=False)

    # LSS making object catalogs
    cat1 = SkyCoord(np.rad2deg(smallrandra1), np.rad2deg(smallranddec1), unit=u.degree, frame='icrs')
    cat2 = SkyCoord(np.rad2deg(smallrandra2), np.rad2deg(smallranddec2), unit=u.degree, frame='icrs')

    id1, id2, d2, d3 = cat2.search_around_sky(cat1, 1/6*u.degree)

    # LSS making skymap
    randra1, randdec1 = randomglobe(n=10000, plot=True, projection='aitoff')
    randra2, randdec2 = randomglobe(n=10000, plot=False)

    # LSS healpy stuff
    ns = 1 # LSS nside

    # LSS convert RA/Dec to Phi/Theta
    randphi1 = randra1
    randphi2 = randra2
    randthet1 = np.pi/2 - randdec1
    randthet2 = np.pi/2 - randdec2

    # LSS get pixels of all the points
    randpx1 = hp.ang2pix(ns, randthet1, randphi1)
    randpx2 = hp.ang2pix(ns, randthet2, randphi2)
    
    # LSS print the pixel size
    print(f'pixel size = {hp.nside2pixarea(ns, degrees=True)} sqdeg')

    # LSS histogram of pixel
    plt.hist(randpx1, histtype='step', density=True, bins=np.arange(13))
    plt.hist(randpx2, histtype='step', density=True, bins=np.arange(13))
    plt.show()

    # LSS Plotting    
    plt.scatter(np.rad2deg(smallrandra1), np.rad2deg(smallranddec1), marker='+', s=10)
    plt.scatter(np.rad2deg(smallrandra2), np.rad2deg(smallranddec2), marker='x', s=10)
    plt.scatter(np.rad2deg(smallrandra1)[id1], np.rad2deg(smallranddec1)[id1], marker='v', color='k', s=10)
    plt.scatter(np.rad2deg(smallrandra2)[id2], np.rad2deg(smallranddec2)[id2], marker='v', color='k', s=10)
    
    plt.xlabel('ra')
    plt.ylabel('dec')
    plt.show()
    print('hello world XD')

