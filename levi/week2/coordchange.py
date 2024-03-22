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
from astropy.coordinates import SkyCoord
import astropy.units as u
import pdb
import time

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


# MAIN
###############################
###############################
if __name__ == '__main__':
    #print('hello world XD')
    ra = '12:30:45.3' # random coordinates
    dec = '+29:50:44'
    obj = SkyCoord(ra, dec, frame='icrs', unit=(u.hourangle, u.deg))
    print(f'RA: {obj.ra.hms} DEC: {obj.dec.deg}')

    print('SkyCoord conversion')
    obj_cart = obj.representation_type='cartesian' # LSS converting to cartesian
    print(obj)
    print('manual conversion')
    obj.representation_type='spherical' # LSS converting back for manual conversion
    x = np.cos(obj.ra.rad)*np.cos(obj.dec.rad)
    y = np.sin(obj.ra.rad)*np.cos(obj.dec.rad)
    z = np.sin(obj.dec.rad)
    print(x, y, z)

    # ra/dec of center of galaxy + what constellation is it in.
    galcenter = SkyCoord(0, 0, frame='galactic', unit=u.deg)
    galctr_radec = galcenter.transform_to('icrs')
    print(f'gal center in ra/dec: {galctr_radec}')
    print(f'gal center in constellation: {galcenter.get_constellation()}')

    # LSS for nashville (lat=36º N) plot (l,b) changes
    ra_arr = np.arange(0, 360, 10)
    dec_arr = np.zeros_like(ra_arr)
    dec_arr.fill(36)

    # LSS make SkyCoord objects + transform to galactic
    c = SkyCoord(ra_arr, dec_arr, frame='icrs', unit='deg')
    galcoord = c.transform_to('galactic')
    
    # LSS Plotting
    plt.scatter(galcoord.l.deg, galcoord.b.deg)
    plt.xlabel('l')
    plt.ylabel('b')
    plt.title('Precession of Zenith in Galactic Coordinates')
    plt.show()

