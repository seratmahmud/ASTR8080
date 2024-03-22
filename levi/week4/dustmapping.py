# L. Schult
# v1 30 Jan 2024
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
import sfdmap
import extinction

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
    dustdir = './sfddata/'
    m = sfdmap.SFDMap(dustdir, scaling=1)
    wave = np.array([3543., 4770., 6231., 7625., 9134.])

    # LSS getting cords + mags
    ra1, dec1 = 246.933, 40.795
    c1 = SkyCoord(ra1, dec1, unit=[u.deg, u.deg]).galactic
    ugriz1 = np.array([18.82, 18.81, 18.74, 18.81, 18.89]) # LSS ugriz magnitudes
    ra2, dec2 = 236.562, 2.440
    c2 = SkyCoord(ra2, dec2, unit=[u.deg, u.deg]).galactic
    ugriz2 = np.array([19.37, 19.10, 18.79, 18.72, 18.62])

    # LSS extinction correcting
    ebv1 = m.ebv(c1.l.value, c1.b.value, frame='galactic') 
    # LSS getting dust map for certain skycoord
    A1 = extinction.fitzpatrick99(wave, 3.1*ebv1)
    # LSS calculating extinction ^ then correcing below
    ugriz1_corr = ugriz1 - A1
    # LSS doing the second star
    ebv2 = m.ebv(c2.l.value, c2.b.value, frame='galactic')
    A2 = extinction.fitzpatrick99(wave, 3.1*ebv2)
    ugriz2_corr = ugriz2 - A2

    # LSS meshgridding
    xlinspace = np.linspace(int(236.6-50), int(236.6+50), 100) # LSS: RA
    ylinspace = np.linspace(int(2.4-50), int(2.4+50), 100) # LSS: DEC
    xmap1, ymap1 = np.meshgrid(xlinspace, ylinspace)
    xlinspace = np.linspace(int(246.9-50), int(246.9+50), int(100/1.3)) # LSS making second grid
    ylinspace = np.linspace(int(40.8-50), int(40.8-50), 100)
    xmap2, ymap2 = np.meshgrid(xlinspace, ylinspace)

    # LSS convert to galactic coords
    map1 = SkyCoord(xmap1, ymap1, unit=(u.degree, u.degree), frame='galactic')
    map2 = SkyCoord(xmap2, ymap2, unit=(u.degree, u.degree), frame='galactic')
    ebvmap1 = 4
    # LSS contour plotting
    

    # LSS plotting color + corrections
    colorcorrectplot = False # LSS change to True to see the plot! I have it 
    # LSS bypassed rn so that I can work on the contour plot.
    if colorcorrectplot == True:
        plt.plot(ugriz1[1]-ugriz1[2], ugriz1[2]-ugriz1[3], marker='.', 
                label='q1 uncorr', color='r')
        plt.plot(ugriz2[1]-ugriz2[2], ugriz2[2]-ugriz2[3], marker='.', 
                label='q2 uncorr', color='r')
        plt.plot(ugriz1_corr[1]-ugriz1_corr[2], ugriz1_corr[2]-ugriz1_corr[3], marker='x', 
                label='q1 corr', color='b')
        plt.plot(ugriz2_corr[1]-ugriz2_corr[2], ugriz2_corr[2]-ugriz2_corr[3], marker='x', 
                label='q2 corr', color='b')
        plt.xlabel('g-r')
        plt.ylabel('r-i')
        plt.legend()
        plt.show()

    zc12 = 1.08
    
    print('hello world XD')

