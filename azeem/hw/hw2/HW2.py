# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pdb
import time
import math
from matplotlib.patches import Rectangle
from numpy.random import random
import astropy.io.fits as fits

from astropy import units as u
from astropy.coordinates import SkyCoord
from collections import Counter
import healpy




def square(ra_min,ra_max,dec_min,dec_max):
    """
    This is a function that finds the area of a square field in square degrees,
    and returns the positions of the input coordiantes.

    Args:
        ra_min (float): A number that represents the lowest RA of your box.
        ra_max (float): A number that represents the highest RA of your box.
        dec_min (float): A number that represents the lowest dec of your box.
        dec_max (float): A number that represents the highest dec of your box.
    Returns:
        xs (list): List of RA's to make plotting a square easier
        ys (list): List of decs to make plotting a square easier
        area (float): The area of your square as a function of the inputs
    """

    #-------------------------------------------------------------
    #+
    # PURPOSE:
    #   Find area of field on the sky in square degrees
    #
    #
    #
    # CALLING SEQUENCE:
    #   square(ra_min, ra_max, dec_min, dec_max)
    #
    # INPUTS:
    #   ra_min, ra_max, dec_min, dec_max - borders of field
    #-
    #-------------------------------------------------------------

    #AB formula for area of lat-lon rectangle bounded by ra and dec, in square degrees
    area = (ra_max*(np.pi/180)-(np.pi/180)*ra_min)*(np.sin(dec_max*(np.pi/180))-np.sin((np.pi/180)*dec_min))*(180/np.pi)**2
    print(area)

    #AB put positions into an array, helps for plotting
    coord = [[ra_min, dec_min], [ra_min, dec_max],[ra_max, dec_max], [ra_max, dec_min]]

    #AB repeat the first point to create a 'closed loop'
    coord.append(coord[0])

    #AB create lists of x and y values
    xs, ys = zip(*coord)

    #return x and y values, as well as the area of the square
    return(xs, ys, area)

def random_pop(ra_min,ra_max,dec_min,dec_max):
    """
    This is a function that randomly populates a square field that is drawn on
    the surface of a sphere.

    Args:
        ra_min (float): A number that represents the lowest RA of your box.
        ra_max (float): A number that represents the highest RA of your box.
        dec_min (float): A number that represents the lowest dec of your box.
        dec_max (float): A number that represents the highest dec of your box.
    Returns:
        ra_good (array): An output numpy array of all the random RA's in the box.
        dec_good (array): An output numpy array of all the random decs in the box.
    """

    #-------------------------------------------------------------
    #+
    # PURPOSE:
    #   Randomly generate RA's and decs in a square field
    #
    #
    #
    # CALLING SEQUENCE:
    #   random_pop(ra_min, ra_max, dec_min, dec_max)
    #
    # INPUTS:
    #   ra_min, ra_max, dec_min, dec_max - borders of field
    #-
    #-------------------------------------------------------------

    #AB generate random points for ra and dec
    ra_rand = (random(10000))*360 #puts into degrees
    dec_rand = (180/np.pi)*np.arcsin(1.-random(10000)*2.)

    #AB make points bounded in the rectangle
    good = np.where((ra_rand>=ra_min) & (ra_rand<=ra_max) & (dec_rand >= dec_min) & (dec_rand <= dec_max))

    #AB grab the random points that are in bounded area
    ra_good = ra_rand[good]
    dec_good = dec_rand[good]

    #return randomly generated ras and decs that were bounded in rectangle
    return ra_good, dec_good

def pix_num(quasar_file):
    """
    This is a function that determines the pixel number at Nside = 4, 8, and 16
    using the HEALpix hierarchy for quasars. It also stores the pixel numbers, RA,
    and dec in a fits file.

    Args:
        quasar_file (data file): .dat file of quasars
    Returns:
        N/A
    """

    #-------------------------------------------------------------
    #+
    # PURPOSE:
    #   Determine pixel number for quasars in quasar_file.
    #
    #
    #
    # CALLING SEQUENCE:
    #   pix_num(quasar_file)
    #
    # INPUTS:
    #   quasar_file - a file containing the RA's and decs of quasars.
    #-
    #-------------------------------------------------------------

    #AB makes empty array to store quasar positions
    all_pos = []

    #AB To open file and make an array:
    with open(quasar_file, 'r') as file:
        for line in file:
            qpos = line
            qpos = f'{qpos[:2]} {qpos[2:4]} {qpos[4:9]} {qpos[9:12]} {qpos[12:14]} {qpos[14:]}'
            cquasar = SkyCoord(qpos, unit=(u.hourangle, u.deg)) #converts to SkyCoord object
            all_pos.append(cquasar) #make array
    all_pos = np.array(all_pos)

    #AB make RAs and decs of quasars into SkyCoord object
    c = SkyCoord(all_pos, unit=(u.hourangle, u.deg))

    #AB define parameters for healpy
    phi = np.radians(c.ra.degree) #want in radians
    theta = np.pi/2 - np.radians(c.dec.degree)
    #AB want pixel numbrt at Nside = 4, 8, and 16
    pix4 = healpy.ang2pix(4,theta,phi)
    pix8 = healpy.ang2pix(8,theta,phi)
    pix16 = healpy.ang2pix(16,theta,phi)

    #AB Make recarray with tags ra, dec, and pixnum
    type = np.dtype([('ra', 'f8'), ('dec', 'f8'), ('pixnum', 'f8', 3)])
    objs = np.zeros(len(c), dtype=type) #empty array to add to later
    #AB convert ra and dec back to degrees
    objs['ra'] = phi * 180./np.pi
    objs['dec'] = c.dec.degree
    #AB make pixnum a 3-array
    objs['pixnum'] = [px for px in zip(pix4, pix8, pix16)]
    #AB convert into fits file format and save
    #maybe...? for x,y,z in zip(pix4, pix8,pix16):
    rec_array = fits.BinTableHDU(objs)
    rec_array.writeto('HW2quasars.fits', overwrite=True)

def quasar_plot(file_name):
    """
    This is a function that plots the RAs and decs of the quasars, as well as
    plotting the location of the five most over-dense regions at Nside = 4.

    Args:
        file_name (data file): .fits file of quasars
    Returns:
        N/A
    """

    #-------------------------------------------------------------
    #+
    # PURPOSE:
    #   Plot quasar positions and show the 5 most over-dense regions
    #
    #
    #
    # CALLING SEQUENCE:
    #   quasar_plot(file_name)
    #
    # INPUTS:
    #   file_name - fits file containing the ra, dec, and pixel numbers of the
    #               quasars.
    #-
    #-------------------------------------------------------------
    #AB set the timer
    time0  = time.time()
    #AB open fits file that was produced with pix_num()
    file = fits.open(file_name)
    hdr = file[0].header
    data = file[1].data

    #AB grab RA, dec and pixel numbers
    ra = data['ra']
    dec = data['dec']
    pixnum = data['pixnum']
    #plot!
    plt.scatter(ra,dec, label='Quasars', color='black')
    #AB make colors and labels for the over densities plot
    colors = ['red', 'orange', 'yellow', 'green', 'blue']
    labels= ['Most Dense', '2nd', '3rd', '4th', 'Least Dense']


    #AB want only the Nside = 4 pixels
    pix4 = [pixels[0] for pixels in pixnum]
    #AB grabs most 5 most populated pixels
    over_dense = Counter(pix4).most_common(5)
    #AB loop and plot the 5 most dense pixels
    for i in range (0, len(over_dense)):
        good = np.where(pix4 == over_dense[i][0])
        plt.scatter(ra[good], dec[good], color=colors[i], label=labels[i])
    plt.xlabel("RA [degrees]")
    plt.ylabel("Dec [degrees]")
    plt.title("Quasar positions")
    plt.legend(fontsize=8)
    #AB time
    time1  = time.time()
    print("    opening of fits file and grabbing 5 most over-dense pixels".format(1e6*(time1-time0)))
    plt.show()





if __name__ == '__main__':
    #1
    #AB set the timer
    time0  = time.time()
    #AB loop over 4 different max declinations
    for i in range(0,4):
        decs = [10,20,30,40] #4 max decs
        #grab points and area from the function square
        xs, ys, area = square(0,60,0,decs[i])
        #AB time
        time1  = time.time()
        print("    calculation of the area and making list of postions".format(1e6*(time1-time0)))
        #plot areas
        plt.figure()
        plt.plot(xs,ys)
        plt.title('Square {}'.format(i+1))
        plt.xlabel('RA [degrees]')
        plt.ylabel('dec [degrees]')
        plt.legend(["Area= {}".format(area)])
        #AB time before interactive plot
        #so it doesn't time plot being open
        time2  = time.time()
        print("    plot [s]: {0:.3f}".format((time2-time1)))
        plt.show()

    #2
    #AB grab points and area from the function square
    xs, ys, area = square(0,30,0,30)
    time3  = time.time()
    print("    calculation of the area and making list of postions".format(1e6*(time3-time2)))
    #AB get ratio of area of square to area of sphere
    area_sphere = (4 * np.pi * (180/np.pi)**2) #area of a sphere
    area_ratio = area/area_sphere
    #AB grab ra and dec random bounded points from random_pop
    ra_good, dec_good = random_pop(0,30,0,30)
    time4  = time.time()
    print("   randomly generating points and calculating area of sphere".format(1e6*(time4-time3)))
    #AB calculate ratio of points in rectangle to points in the sky
    rect_ratio = len(ra_good) / 10000
    print("area ratio= ", area_ratio)
    print("rectangle ratio= ", rect_ratio)

    #3
    #AB make fits file with ra's, decs and pixel numbers of quasars
    time5  = time.time()
    print("    making fits quasar file".format(1e6*(time5-time4)))
    pix_num('HW1quasarfile.dat')

    #4
    #AB plot quasar positons and where the 5 most dense regions are located
    quasar_plot('HW2quasars.fits')
