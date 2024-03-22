# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord
import pymangle
import sys
sys.path.insert(0,'/home/baria/')
from sphericalcaps import ra_bound, dec_bound

def four_poly(ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4):

    #-------------------------------------------------------------
    #+
    # PURPOSE:
    #   Gives area for a "lat-long" rectangle given 4 caps that are RA and dec
    #   bounded.
    #
    #
    #
    # CALLING SEQUENCE:
    #   four_poly(ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4)
    #
    # INPUTS:
    #   ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4 - RAs and decs for each of the 4 caps
    #-
    #-------------------------------------------------------------

    #call spherical cap functions
    ra_bound(ra1, dec1)
    ra_bound(ra2, dec2)
    dec_bound(ra3, dec3)
    dec_bound(ra4,dec4)

    #calculate area
    area = (ra2*(np.pi/180)-(np.pi/180)*ra1)*(np.sin(dec4*(np.pi/180))-np.sin((np.pi/180)*dec3))
    print(area)

if __name__ == '__main__':

    #prints (x,y,z, 1-h) for each cap and gives RA. Polygon files made by copying
    #and pasting these values
    four_poly(5*15,0,6*15,0,0,30,0,40)
    four_poly(10*15,0,12*15,0,0,60,0,70)
