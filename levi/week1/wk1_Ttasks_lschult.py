# L. Schult
# v1 
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
from matplotlib import rc
import astropy
import pdb
import time
from astropy.io import fits

# FUNCTIONS
###############################
###############################
# def radecplot(ras, decs):
#     """
#     This is a function that plots sky positions according to input ra/dec arrays

#     Args:
#         ras (ndarray): An array of ras
#         decs (ndarray): An array of decs
#     Returns:
#         nothing - ends by showing matplotlib
#     """
# #-------------------------------------------------------------
# #+
# # PURPOSE:                                                 
# #   basic sky plotting function
# #
# # CALLING SEQUENCE:
# #   radecplot(ras, decs) 
# #
# # INPUTS:
# #   ra, dec - plotted using plt.plot()
# #-
# #-------------------------------------------------------------  
#     plt.plot(ras, decs, 'bx')
#     plt.show()  

# MAIN
###############################
###############################
if __name__ == '__main__':
    fx = fits.open('~/Desktop/classes/ASTR8080TEMP/week1/struc.fits')
    objs = fx[1].data # LSS reading in files
    hdr = fx[1].header
    #plt.plot(objs['RA'], objs['DEC'], 'bx')
    #plt.show()
    ra = objs['RA']
    dec = objs['DEC']
    extinct = objs['EXTINCTION'][:,0]
    #mask = np.where(extinct > 0.22)[0] # masking for extinction > 0.22
    mask = extinct > 0.22
    ranotex = ra[~mask] # LSS making two separate RA+DEC to plot objects with 
    raex = ra[mask] # LSS extinction > 0.22 as a different color on the plot. 
    decnotex = dec[~mask]
    decex = dec[mask]
    plt.plot(ranotex, decnotex, 'bx', label='ext<0.22')
    plt.plot(raex, decex, 'rx', label='ext>0.22')
    plt.legend()
    plt.title('object locations')
    plt.xlabel('Right Ascension (º)')
    plt.ylabel('Declination (º)')
    plt.show()


