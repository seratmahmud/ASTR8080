# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pdb
import time

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates import AltAz
from astropy.coordinates import Galactic, Longitude



if __name__ == '__main__':
    #1
    #AB get our coordinates
    c = SkyCoord(ra='10h05m02s', dec='+10d05m02s', frame='icrs', unit=(u.hourangle, u.deg))
    print("skycoord: ",c)
    #AB convert them to cartesian
    c.representation_type = 'cartesian'
    print("cartesian: ", c)

    #2
    #AB input coordinates of galactic center using galactic coordinates
    gc = SkyCoord(l="0h0m0s", b='0h0m0s', frame=Galactic)
    print("galactic center: ", gc)
    #AB convert them to equatorial coordinates
    icrs_gc = gc.transform_to("icrs")
    print("galactic center in icrs: ",icrs_gc)

    #AB find constellation gc is in
    con = gc.get_constellation()
    print("what constellation is the gc in? ",con)
    #AB the galactic center is very near the edge of the constellation

    #3
    #AB make list of ra's from 0 to 360 degrees
    ras = np.arange(0,360)
    ra= Longitude(ras, unit=u.deg)
    #AB make list of declination of Nashville throughout the year (stays the same)
    nashdec = [36] * 360
    nashcoords = SkyCoord(ra, nashdec, frame = 'icrs', unit='deg')
    #AB convert to galactic coordinates
    nashgal = nashcoords.galactic
    nashgal.l.degree
    nashgal.b.degree
    #AB plot
    plt.scatter(nashgal.l.degree, nashgal.b.degree)
    plt.xlabel("Galactic Longitude")
    plt.ylabel("Galactic Latitude")
    plt.show()
