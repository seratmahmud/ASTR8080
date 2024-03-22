# L. Schult
# v1 18 jan 2024
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
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time


# Task:
# Use Skycoord from astropy.coordinates to convert a dec
# in (o,’,”) format to decimal degrees. Do the same for an
# RA in hms format

# Use Time.now() from astropy.time to obtain today’s
# MJD and today’s JD

# FUNCTIONS
###############################
###############################
def getradecdecimal(ra, dec):
    """
    This is a function that converts RA in hms to decimal degrees

    Args:
        ra (string): input ra as a string in hms coordinates e.g. '00h42m30s'
        dec (string): input dec as a string in d'" coordinates e.g. '+41d12m00s'
        uses astropy.SkyCoord -- see documentation for input options.
    Returns:
        radecimal, decdecimal (float): output ra and dec in decimal degrees
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   convert ra(dec) from hours(degrees)minutes(')seconds(") to decimal degrees
#
# CALLING SEQUENCE:
#   radecimal, decdecimal = getradecdecimal(ra, dec)
#
# INPUTS:
#   ra, dec
#-
#-------------------------------------------------------------  
    c = SkyCoord(ra, dec)
    return c.ra.degree, c.dec.degree


def dateto_mjdjd(date):
    """
    A function that prints today's MJD and JD. 

    Args:
        date (string): date in the form '2000-08-15' year-mo-da

    Returns:
        mjd, jd (float): MJD and JD of input date 
    
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   convert date from year-mo-da to MJD and JD
#
# CALLING SEQUENCE:
#   datemjd, datejd = dateto_mjdjd(date)
#
# INPUTS:
#   date
#-
#-------------------------------------------------------------  
    t = Time(date)
    return t.mjd, t.jd

# MAIN
###############################
###############################
if __name__ == '__main__':
    # LSS example coordinate conversion using random coords.
    rahms = '18h36m24s' 
    decdms = '23d45m18s'
    radecimal, decdecimal = getradecdecimal(rahms, decdms)
    print(f'coordinates: {rahms}, {decdms}')
    print(f'converted {radecimal}º, {decdecimal}º')

    # LSS example get mjd/jd from today's date
    today = Time.now()
    todaymjd, todayjd = dateto_mjdjd(today)
    print(f'today is {today}')
    print(f'today\'s MJD is {todaymjd}')
    print(f'today\'s JD is {todayjd}')


