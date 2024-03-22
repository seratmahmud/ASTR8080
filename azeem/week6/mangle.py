# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u
from astropy.coordinates import SkyCoord
import pymangle

if __name__ == '__main__':
    #task3
    #AB read in ply file and generate random points
    minter = pymangle.Mangle("intersection.ply")
    mboth = pymangle.Mangle("bothcaps.ply")
    ra_rand, dec_rand = minter.genrand(10000)
    ra_rand2, dec_rand2 = mboth.genrand(10000)

    #plot!
    fig = plt.figure()
    ax1 = fig.add_subplot()
    ax1.scatter(ra_rand2, dec_rand2)
    ax1.scatter(ra_rand, dec_rand)
    fig.savefig('intersect.png')

    #AB flip the sign and read in ply file
    mflip1 = pymangle.Mangle("mflip.ply")
    ra_rand3, dec_rand3 = mflip1.genrand(10000)

    #plot!
    fig2 = plt.figure()
    ax2 = fig.add_subplot()
    ax2.scatter(ra_rand3, dec_rand3)
    ax2.scatter(ra_rand, dec_rand)
    fig2.savefig('flip1.png')
