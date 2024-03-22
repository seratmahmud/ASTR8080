# M. Rizzo Smith
# v1 2/13/24
# Week 6 Lecture 1 Tasks
# ASTR 8080 spherical caps


# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
import astropy 
from matplotlib import rc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic
from astropy.coordinates import get_body
from astropy.time import Time

# FUNCTIONS
###############################
###############################

def ra_cap(ra):
    # Dont forget the bobble
    ra = ra + 6 # Adding 6h because ra cap defined by +90 degrees from RA bound
   # print(str(int(ra)).zfill(2))
    c = SkyCoord(ra, 0.0, unit=(u.hourangle, u.degree))
    c_4array = [c.cartesian.x.value, c.cartesian.y.value, c.cartesian.z.value, 1]
 
    return c_4array

def dec_cap(dec):
    dec_cap = SkyCoord('00 00 00.00 90 00 00', unit=(u.hourangle, u.degree))
    d_cap_4 = [dec_cap.cartesian.x.value, dec_cap.cartesian.y.value, dec_cap.cartesian.z.value, 1-np.sin(dec * np.pi/180)]
    return d_cap_4

def gen_cap(ra, dec, r):
    c = SkyCoord(ra, dec, unit=(u.degree, u.degree))
    cap = [c.cartesian.x.value, c.cartesian.y.value, c.cartesian.z.value, 1-np.cos(np.radians(r))]
    return cap
# MAIN
###############################
###############################
if __name__ == '__main__':
    print(ra_cap(5))
    print(dec_cap(36))
    print(gen_cap(5, 36, 1))
