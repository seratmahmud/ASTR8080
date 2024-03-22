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
from astropy.coordinates import SkyCoord, Angle
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


def racap(rabound):
    """
    This is a function that gets the 4 vector in xyz(1-h)
    rabound should be '5h' string that is in hourangle
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
    racenter = Angle(rabound) + 6*u.hourangle
    racenter = SkyCoord(racenter, 0, unit=u.hourangle, frame='icrs')
    centervector = [racenter.cartesian.x.value, racenter.cartesian.y.value, 
                    racenter.cartesian.z.value, 1]
    # racenter = racenter
    return centervector

def deccap(decbound):
    """
    This is a function that gets the 4 vector in xyz(1-h)
    decbound should be in '30d' in degrees
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
    decbound_ang = Angle(decbound) # LSS convert to a degree value
    deccenter = SkyCoord(0, 90, unit=[u.hourangle, u.degree], frame='icrs')
    centervector = [deccenter.cartesian.x.value, deccenter.cartesian.y.value, 
                     deccenter.cartesian.z.value, 1-np.sin(decbound_ang.radian)]
    return centervector

def radeccap(ra, dec, theta):
    """
    This is a function that gets the 4 vector in xyz(1-h)
    format = '02h35m15s', '25d51m20s'
    theta in degrees
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
    capcenter = SkyCoord(ra, dec, unit=[u.hourangle, u.degree], frame='icrs')
    centervector = [capcenter.cartesian.x.value, capcenter.cartesian.y.value, 
                      capcenter.cartesian.z.value, 1-np.cos(np.deg2rad(theta))]
    return centervector

# MAIN
###############################
###############################
if __name__ == '__main__':
    print(radeccap('05h00m00s','36d00m00s', 1))
    print('hello world XD')

