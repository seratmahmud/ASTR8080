# M. Rizzo Smith
# v1 1/25/24
# Week 6 Lecture 2 Tasks
# ASTR 8080 Mangle

# ASTR 8080
# to import this from another directory:
# import sys
# sys.path.insert(0, ‘../week1/’)
# import polycalc
# from polycalc import get_poly_o3


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
import sys 
sys.path.insert(0, '../week6/')
import spherical_caps
from spherical_caps import gen_cap
import pymangle

# FUNCTIONS
###############################
###############################

# MAIN
###############################
###############################
if __name__ == '__main__':
    cap1 = gen_cap(76, 36, 5)
    cap2 = gen_cap(75, 35, 5)
    
    #print(cap1)
    #print(cap2)
    # np.savetext('intersection.ply')
    minter = pymangle.Mangle('intersection.ply')
    ra_rand, dec_rand = minter.genrand(10000)
    
    mflip = pymangle.Mangle('flip.ply')
    ra_flip, dec_flip = mflip.genrand(10000)
    
    mflip2 = pymangle.Mangle('flip2.ply')
    ra_flip2, dec_flip2 = mflip2.genrand(10000)

    bof = pymangle.Mangle('bothcaps.ply')
    ra_bof, dec_bof = bof.genrand(10000)
    plt.plot(ra_bof, dec_bof, 'o', color='lightblue')
    plt.plot(ra_rand, dec_rand, 'o', color='green')
    plt.plot(ra_flip, dec_flip, 'o', color='teal')
    plt.plot(ra_flip2, dec_flip2, 'o', color ='lightblue')
    plt.savefig('int-flip1-flip2.png', dpi=300)
    plt.show()
    
    mdubflip = pymangle.Mangle('dubflip.ply')
    ra_dub, dec_dub = mdubflip.genrand(1000000)
    
    plt.plot(ra_dub, dec_dub, 'o', color='tomato')
    plt.savefig('double-flip.png', dpi=300)
    plt.show()
