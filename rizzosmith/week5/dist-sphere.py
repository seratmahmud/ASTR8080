# M. Rizzo Smith
# v1 1/25/24
# Week 3 Lecture 2 Tasks
# ASTR 8080 coordinate system transformations


# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
import astropy 
from matplotlib import rc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from numpy.random import random
# FUNCTIONS
###############################
###############################

# None for now

# MAIN
###############################
###############################
if __name__ == '__main__':
    obj1 = SkyCoord(263.75, -12.9, frame='icrs', unit='deg')
    obj2 = SkyCoord('20 24 59.9 10 06 00.0' ,frame='icrs', unit=(u.hourangle, u.degree))
    a = [obj1.cartesian.x, obj1.cartesian.y, obj1.cartesian.z]
    b = [obj2.cartesian.x, obj2.cartesian.y, obj2.cartesian.z]
    dot = np.dot(a, b)
    mag1 = np.sqrt(a[0]**2 + a[1]**2 + a[2]**2)
    mag2 = np.sqrt(b[0]**2 + b[1]**2 + b[2]**2)

    z_ang = np.arccos(dot / (mag1 * mag2))
    print(z_ang)
   
    z_ang_sep = obj1.separation(obj2)
    print(z_ang_sep.rad)

    # Generate random points between 2 and 3 hours, will do this better for homework
    ra1 = (random(100) * 15) + 30
    ra2 = (random(100) * 15) + 30

    dec1 = (random(100)*4) - 2
    dec2 = (random(100)*4) - 2
    cat1 = SkyCoord(ra1, dec1, unit=(u.deg, u.deg))
    cat2 = SkyCoord(ra2, dec2, unit=(u.deg, u.deg))
    rad= (1/6) * u.degree
    ind1, ind2, d2d, d3d = cat2.search_around_sky(cat1, rad)

    
    plt.plot(ra1, dec1, 'o', color='red')
    plt.plot(ra2, dec2, 'x', color='orange')
    plt.plot(ra1[ind1], dec1[ind1], 'v', color='blue')
    plt.plot(ra2[ind2], dec2[ind2], 'v', color='blue')
    plt.show()

    big_ra = np.append(ra1, ra2)
    big_dec = np.append(dec1, dec2)
    PC = SkyCoord('02 20 05 -00 06 12', frame='icrs', unit=(u.hourangle))

