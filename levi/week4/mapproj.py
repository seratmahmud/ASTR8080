# L. Schult
# v1 1 Feb 2024
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
from numpy.random import random
import matplotlib.pyplot as plt
from matplotlib import rc
import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
import pdb
import time

# FUNCTIONS
###############################
###############################
def randomglobe(projection='xy', n=10000):
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
    ra = 2*np.pi*(random(10000)-0.5)
    dec = np.arcsin(1.-random(10000)*2)
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
        ax.scatter(ra, dec, marker=',', s=1)
        #ax.set_xticks()
        plt.show()


# MAIN
###############################
###############################
if __name__ == '__main__':
    randomglobe(projection='lambert')
    print('hello world XD')

