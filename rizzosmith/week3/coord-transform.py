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
from astropy.coordinates import Galactic
from astropy.coordinates import get_body
from astropy.time import Time

# FUNCTIONS
###############################
###############################

# None for now

# MAIN
###############################
###############################
if __name__ == '__main__':

    # MRS initialize test object
    obj = SkyCoord('05 46 00 +28 56 00', unit=(u.hourangle, u.degree))
    print(obj.cartesian)


    # MRS Check cartesian transform is right according to math in slides 
    x_math = np.cos(obj.ra.radian)*np.cos(obj.dec.radian)
    print(x_math, obj.cartesian.x)

    y_math = np.sin(obj.ra.radian)*np.cos(obj.dec.radian)
    print(y_math, obj.cartesian.y)

    z_math = np.sin(obj.dec.radian)
    print("z: ", z_math, obj.cartesian.z )

    # MRS set the galactic center coordinate in coordinates
    GC = SkyCoord(0, 0, frame=Galactic, unit=(u.degree, u.degree))
    print(GC.l, GC.b)

    # MRS transform to ICRS to get ra and dec of GC
    GC = GC.transform_to('icrs')
    print(GC.ra.hms, GC.dec.degree)
    print(GC.get_constellation())
    
    # MRS Create a list of RA from 0 to 360 
    ra_coords = np.arange(0, 360)
    
    # MRS Create matching Dec for the Zenith at 36N from Nashville 
    dec_coords = (ra_coords * 0) + 36.0
    
    # MRS Initialize sky coord object in ICRS for list of RA and Dec throughout the year
    objs_year = SkyCoord(ra_coords, dec_coords, unit=(u.degree, u.degree))
    
    # MRS Transform to Galactic coordinates
    objs_year = objs_year.galactic
    
    # MRS Create and display a plot of galactic lattiude and longitude of zenith
    plt.plot(objs_year.l.degree, objs_year.b.degree, 'o')
    plt.xlabel('l')
    plt.ylabel('b')
    plt.show()
    
    # MRS This is gonna be such a nice optimized function for HW1 but alas for now I just want to plot it quick.

    
    # MRS Initializing planet positions
    now = Time.now()
    
    # MRS Pull each planet sky coord position for now
    merc = get_body('mercury', now)
    venus = get_body('venus', now)
    mars = get_body('mars', now)

    # MRS Transform into the ecliptic coordinate system
    mars = mars.transform_to('heliocentrictrueecliptic')
    venus = venus.transform_to('heliocentrictrueecliptic')
    merc = merc.transform_to('heliocentrictrueecliptic')

    # Plot ecliptic lat and lon of each planet
    plt.plot(mars.lat.degree, mars.lon.degree, 'o', color = 'red', label = 'Mars')
    plt.plot(venus.lat.degree, venus.lon.degree, 'o', color = 'orange', label = 'Venus')
    plt.plot(merc.lat.degree, merc.lon.degree, 'o', color = 'black', label='Mercury')
    plt.legend()
    plt.xlabel('Ecliptic Lattitude')
    plt.ylabel('Ecliptic Longitude')
    plt.show()
