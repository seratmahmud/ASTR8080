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

import healpy

if __name__ == '__main__':
    #task1
    #AB Generate random points
    ra = 360*(np.random.random(1000000))
    dec = (180/np.pi)*np.arcsin(1-np.random.random(1000000)*2)
    #AB plot them
    fig = plt.figure()
    ax = fig.add_subplot(111,projection="aitoff")
    ax.scatter(ra,dec, s=1)
    ax.grid(color='blue', linestyle='dashed', linewidth=3)
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'] #convert to hms
    ax.set_xticklabels(xlab, weight=800)
    fig.show()
    pdb.set_trace()

    #task2
    #AB define your parameters for ang2pix
    nside = 1
    phi = np.radians(ra)
    theta = np.pi/2 - np.radians(dec)
    pix = healpy.ang2pix(nside,theta,phi)
    #calculate area of a pixel
    area = (4*np.pi)/12
    area *= (180/np.pi)**2

    print(area)

    #task3
    #AB make histogram to compare
    plt.hist(pix, 12)
    plt.show()
