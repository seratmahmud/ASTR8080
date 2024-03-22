# L. Schult
# v1 
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
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
import pdb
import time
import sys
sys.path.insert(0, '../week6/')
import cap
from cap import radeccap
import pymangle


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

def twoply_writer(numpolygon, capvectors, polygonvects, filename):
    """
    This is a function that writes mangle polygon files (ending in .ply)

    Args:
        numpolygon (int): number of polygons
        capvectors (ndarray): should be n x 4 ndarray with cap 4-vectors
        polygonvects (list): len(polygonvects) = numpolygon. each element of
        the list is itself a list that has the index of the cap each polygon will take
        filename (str): polygon filename no ply needed
    Returns:
        
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
    with open(filename+'.ply', 'w') as outf:
        first_line = str(numpolygon) + ' polygons\n'
        outf.write(first_line)
        for pidx in range(numpolygon): # LSS taking care of 0 indexing
            polystr = f'polygon {pidx+1} ( {len(polygonvects[pidx])} caps, 1 weight, 0 pixel, 0 str):\n'
            outf.write(polystr) # LSS print first line with number of caps based
            # LSS on length of list of vector indices from polygonvects
            for cap_idx in polygonvects[pidx]:
                vectstr = str(capvectors[pidx][cap_idx]).replace('[', '')
                vectstr = vectstr.replace(']', '')
                vectstr = vectstr.replace(',', '')
                outf.write(vectstr+'\n')
        
        outf.close()
    return

# MAIN
###############################
###############################
if __name__ == '__main__':
    #rahms = Angle('76d').hms
    #print(str(rahms[0])+'h'+str(rahms[1])+'m'+str(rahms[2])+'s')
    a = SkyCoord(76, 36, unit=[u.degree, u.degree], frame='icrs')
    b = SkyCoord(75, 35, unit=u.degree, frame='icrs')

    cap1 = radeccap(a.ra.to_string(u.hour), a.dec.to_string(u.degree), 5)
    print(cap1)
    cap2 = radeccap(b.ra.to_string(u.hour), b.dec.to_string(u.degree), 5)
    print(cap2)
    minter = pymangle.Mangle('intersection.ply')
    mboth = pymangle.Mangle('bothcaps.ply')
    mnotint = pymangle.Mangle('NOTintersection.ply')

    minter_ra_rand, minter_dec_rand = minter.genrand(10000)
    mboth_ra_rand, mboth_dec_rand = mboth.genrand(10000)
    mnoti_ra_rand, mnoti_dec_rand = mnotint.genrand(10000)

    plt.plot(mboth_ra_rand, mboth_dec_rand, color='C0', marker=',')
    plt.plot(minter_ra_rand, minter_dec_rand, color='g', marker=',')
    plt.savefig('./mboth_minter_rand.png')
    plt.show()

    plt.plot(minter_ra_rand, minter_dec_rand, color='g', marker=',')
    plt.scatter(mnoti_ra_rand, mnoti_dec_rand, color='orange', marker=',')
    plt.savefig('./minter_mflip_rand.png')
    plt.show()
    print('hello world XD')

