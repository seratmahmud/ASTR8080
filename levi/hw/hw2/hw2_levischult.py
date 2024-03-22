# L. Schult
# v1 12 Feb 2024
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
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from numpy.random import random
from matplotlib import rc
import astropy
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
import pdb
import time

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

def skyarea(ramin, ramax, decmin, decmax):
    """
    This is a function that calculates area of given coords(deg) in sq. deg.'

    Args:
        ramin, ramax (float): A min(max) RA in degrees for which to calculate area
        decmin, decmax (float): a min(max) dec in degrees for which to calc. area
    Returns:
        area (float): the area of the rectangle specified in square degrees
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   calculate area of a rectangle bounded by ra + dec lines
#
# CALLING SEQUENCE:
#   area = skyarea(20, 30, 10, 15) 
#
# INPUTS:
#   ramin, ramax, decmin, decmax
#-
#-------------------------------------------------------------  
    # LSS calculate area in sqdeg - formula from slides
    areasqdeg = ((180 / np.pi)**2)*(np.deg2rad(ramax)-np.deg2rad(ramin))*(
        np.sin(np.deg2rad(decmax)) - np.sin(np.deg2rad(decmin)))
    return areasqdeg

def randomglobe(rarange=[0, 360], decrange=[-90, 90], 
                projection='xy', n=10000, plot=True):
    """
    This is a function that generates random points and plots them.
    projections available are xy, aitoff, and lambert.

    Args:
        projection (string): the projection the map will be in. default is xy
        n (int): the number of random points to generate. default is 10000
    Returns:
        y (datatype): A ...
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   plot random points on different maps
#
# CALLING SEQUENCE:
#   randomglobe(projection='aitoff') 
#
# INPUTS:
#   projection, n
#-
#-------------------------------------------------------------  
    rarange = np.deg2rad(rarange)
    decrange = np.deg2rad(decrange)
    ra = np.rad2deg(((rarange[1]-rarange[0])*(random(n)))+rarange[0])
    dec = np.rad2deg(np.arcsin(((np.sin(decrange[1])-np.sin(decrange[0]))*random(n))+np.sin(decrange[0])))
    if plot:
        if projection=='xy':
            plt.scatter(ra, dec, marker=',', s=1)
            plt.title('random points')
            plt.xlabel('ra')
            plt.ylabel('dec')
            plt.show()
        elif projection == 'aitoff' or 'lambert':
            fig = plt.figure()
            xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
            ax = fig.add_subplot(111, projection=projection)
            ax.set_xticklabels(xlab, weight=800)
            ax.grid(color='b', linestyle='dashed', linewidth=3)
            ax.scatter(ra-np.pi, dec, marker=',', s=1)
            ax.set_title('random points')
            #ax.set_xticks()
            plt.show()
    
    return ra, dec

def densepix(filein, nsi=4, plot=True):
    """
    This is a function that finds the healpy pixels (nside specified) that have
    the most objects based on a fits file that is read in. Will plot all objects
    and those in the denser pixels in a different color/marker 
    NOTE: assumes fits file is one array with first axis being RA, DEC, 4nsidepx, 8nsidepix, 16nsidepix
    Args:
        filein (string): fits filename
        nsi (int): must be power of 2. note that function assumes the structure
        of different nside configurations in the array
        plot (bool): whether or not to plot the points and those in the densest 
        pixels
    Returns:
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   find healpy pixels with most objects from a fitsfile with ra,dec,nside4pix, nside8pix, nside16pix
#
# CALLING SEQUENCE:
#   densepix('starpositions_pixels.fits')
#
# INPUTS:
#   filein, nsi
#-
#-------------------------------------------------------------  
    # LSS reading in file
    with fits.open(filein) as fitsin:
        data = fitsin[0].data
    # LSS getting index of nside column
    nsindex = int(np.emath.logn(2, nsi))
    # LSS counting quasars/bin
    qspix = data[nsindex,:]
    bins, counts = np.unique(qspix, return_counts=True)
    # LSS finding the bins with the most quasars
    countsorder = np.argsort(counts)[::-1]
    overdense_bins = bins[countsorder]
    od5bins = overdense_bins[:5]
    
    # LSS find quasars that are in 5 overdense pixels
    od5bins_mask = np.isin(qspix, od5bins)
    ramasked = data[0,:][od5bins_mask]
    decmasked = data[1,:][od5bins_mask]

    if plot:
        plt.scatter(data[0,:], data[1,:], marker=',', label='all points')
        plt.scatter(ramasked, decmasked, color='red', marker='x', label='points in 5 densest pixels')
        plt.xlabel('Right Ascension (º)')
        plt.ylabel('Declination (º)')
        plt.legend()
        plt.show()
    return



# MAIN
###############################
###############################
if __name__ == '__main__':

    ###### LSS Q1 ########
    # LSS verify area works for spherical caps:
    print(f'my sky area function for half sky = {skyarea(0, 360, 0, 90)}')
    print(f'All sky area = 41252.96 sqdeg / 2 = {41252.96/2}')
   
    # LSS problem 1 declarations
    a0 = skyarea(ramin=0, ramax=90, decmin=0, decmax=20)
    a1 = skyarea(ramin=0, ramax=90, decmin=0, decmax=40)
    a2 = skyarea(ramin=0, ramax=90, decmin=0, decmax=60)
    a3 = skyarea(ramin=0, ramax=90, decmin=0, decmax=80)
    boxareas = [a0, a1, a2, a3]
    
    xy4rect = [(0,0), (90, 0), (180, 0), (270, 0)]
    rectlist = []
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for idx, ht in enumerate(range(20, 100, 20)):
        ax.add_patch(Rectangle(xy4rect[idx], width=90, height=ht, color=f'C{idx}', 
                               label=f'height = {ht}deg, area={boxareas[idx]}',
                               alpha=0.2))
    
    #ax.add_collection(pc)
    ax.set_ylim(-90, 90)
    ax.set_xlim(0, 360)
    ax.legend()
    ax.set_xlabel('Right Ascension (º)')
    ax.set_ylabel('Declination (º)')
    ax.set_title('areas of increasing declination. Offset for ease of viewing')
    plt.show()

    ####### LSS Q2 ########
    # LSS generate random points for whole globe
    npts = 100000
    rand_ra, rand_dec = randomglobe(n=npts)

    # LSS calculate area of rectangle and sky
    # LSS rectangle will be (0,0) with width=90 and height=20 º 
    rectarea = skyarea(0, 90, 0, 20)
    allskyarea = skyarea(0, 360, -90, 90)

    # LSS count points in rectangle
    ptsinrect = 0

    for idx, point in enumerate(zip(rand_ra, rand_dec)):
        if point[0] >= 0 and point[0] <= 90 and point[1] >= 0 and point[1] <= 20:
            ptsinrect += 1 # LSS if point is within height/width of rectangle add 1

    print(f'Points in rect/sky area percent    = {ptsinrect / npts}')
    print(f'Percentage of rect area / sky area = {rectarea / allskyarea}')
    


    ####### LSS Q3 ########
    # LSS below is from hw1
    # LSS load in data + parse raw text into skycoords
    with open('./HW1quasarfile.dat', 'r') as infile:

        rawtxt = np.loadtxt(infile, dtype=str)
    coordstrings = [] 
    strtim = time.time()
    for line in rawtxt:
        coordstrings.append(f'{line[:2]}h{line[2:4]}m{line[4:9]} {line[9:12]}d{line[12:14]}m{line[14:]}s')
    quasarcoords = SkyCoord(coordstrings)

    # LSS change to phi/theta coords for hp ang2pix
    quasarcoords.representation_type = 'physicsspherical'
    quasarpix = np.zeros((3, quasarcoords.shape[0]))
    nsides = [4, 8, 16]
    # LSS getting pix num for 3 different nsides
    for idx, ns in enumerate(nsides):
        quasarpix[idx, :] = hp.ang2pix(ns, quasarcoords.theta.rad, 
                                       quasarcoords.phi.rad)
    
    # LSS creating combined array 
    # LSS RA, DEC, nside4 pix, nside8 pix, nside16 pix
    quasarcoords.representation_type = 'spherical'
    combined_coordspix = np.zeros((5, quasarcoords.shape[0]))
    combined_coordspix[0, :] = quasarcoords.ra.deg
    combined_coordspix[1, :] = quasarcoords.dec.deg
    combined_coordspix[2:, :] = quasarpix
    
    # LSS write to fits file
    hdu_out = fits.PrimaryHDU(combined_coordspix)
    hdu_out.writeto('quasarcoords_pix.fits', overwrite=True)

    ####### LSS Q4 #########
    densepix('quasarcoords_pix.fits')

    #print('hello world XD')

