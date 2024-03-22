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
import sys
sys.path.insert(0, '../week6/')
import cap
from cap import racap, deccap, radeccap
from manglefun import twoply_writer
sys.path.insert(0, '../hw/hw2/')
from hw2_levischult import skyarea


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
    ra5cap = racap('5h')
    ra6cap = racap('6h')
    dec30cap = deccap('30d')
    dec40cap = deccap('40d')
    fourcaps1 = [ra5cap, ra6cap, dec30cap, dec40cap]

    # LSS find area1 in sqdeg
    latlonrect1 = skyarea(ramin=5*15, ramax=6*15, decmin=30, decmax=40)
    latlonrect1 = latlonrect1 * ((np.pi / 180)**2)
    print(f'{latlonrect1=}')

    ra10cap = racap('10h')
    ra12cap = racap('12h')
    dec60cap = deccap('60d')
    dec70cap = deccap('70d')
    fourcaps2 = [ra10cap, ra12cap, dec60cap, dec70cap]

    # LSS find area2 in sqdeg
    latlonrect2 = skyarea(ramin=10*15, ramax=12*15, decmin=60, decmax=70)
    latlonrect2 = latlonrect1 * ((np.pi / 180)**2)
    print(f'{latlonrect2=}')

    # LSS write caps to ply file
    twoply_writer(2, [fourcaps1, fourcaps2], [[0, 1, 2, 3], [0, 1, 2, 3]], 'latlonrects')


    print('hello world XD')

