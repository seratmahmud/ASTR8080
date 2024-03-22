# M. Kaldor
# v1 2/13/2024
# ASTR 8080 HW2

# to import this from another directory:
# import sys
# sys.path.insert(0, '../week1/')
# import hw2
# from hw2 import function

# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt
from matplotlib import rc
import time
import pdb
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import sfdmap
import extinction
from numpy.random import random
from matplotlib.patches import Rectangle
import math
import healpy as hp
import astropy.io.fits as fits
import collections

# FUNCTIONS
###############################
###############################
def square_area(ra_min, ra_max, dec_min, dec_max):
    """
    This is a function that takes in right ascension and declination borders and calculates the area that that square
    covers. It returns the area and the RA and dec converted into radians so that you can plot the square once you've
    run the function.

    Args: ra_min = [float]
        ra_max = [float]
        dec_min = [float]
        dec_max = [float]

    Returns: area = [float]
        ra_min = [float]
        ra_max = [float]
        dec_min = [float]
        dec_max = [float]

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: Calculating the area and providing the necessary values to plot the square on an Aitoff projection
    #   (or another projection that requires radian values to plot)
    #
    #
    # CALLING SEQUENCE: square_area(0, 360, -90, 90)
    #
    #
    # INPUTS: ra_min = minimum right ascension value in degrees
    #   ra_max = maximum right ascension value in degrees
    #   dec_min = declination value in degrees
    #   dec_max = maximum declination value in degrees
    #
    # -------------------------------------------------------------
    #
    # MEK calculate area of square from equation in the class notes
    area = (180/np.pi)*(ra_max-ra_min)*(np.sin(np.radians(dec_max))-np.sin(np.radians(dec_min)))

    # MEK convert RAs and DECs to radians
    ra_min = np.radians(ra_min)
    ra_max = np.radians(ra_max)
    dec_min = np.radians(dec_min)
    dec_max = np.radians(dec_max)

    return(area, ra_min, ra_max, dec_min, dec_max)


###############################
def random_pop(ra_min, ra_max, dec_min, dec_max, n):
    """
    This is a function randomly (by area) populates a square of given bounds with a number of points and returns the
    coordinates of these random points for plotting.

    Args: ra_min = [float]
        ra_max = [float]
        dec_min = [float]
        dec_max = [float]
        n = [int]

    Returns: ra = [list of floats]
        dec = [list of floats]

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: randomly distributing points across the sky by area
    #
    #
    # CALLING SEQUENCE: random_pop(0, 360, -90, 90, 1000)
    #
    #
    # INPUTS: ra_min = minimum right ascension value in degrees
    #   ra_max = maximum right ascension value in degrees
    #   dec_min = declination value in degrees
    #   dec_max = maximum declination value in degrees
    #   n = number of points to put in the specified area
    #
    # -------------------------------------------------------------
    #
    # MEK convert RAs and DECs to radians
    ra_min = np.radians(ra_min)
    ra_max = np.radians(ra_max)
    dec_min = np.radians(dec_min)
    dec_max = np.radians(dec_max)

    # MEK calculate the ranges to populate, create random points within that range (multiply by range, add minimum
    # value of the range to scale random numbers from 0 to 1 to be random numbers within RA and dec ranges)
    ra_range = ra_max-ra_min
    dec_range = dec_max-dec_min
    ra = random(n)*ra_range+ra_min
    dec = random(n)*dec_range+dec_min

    return np.degrees(ra), np.degrees(dec)

###############################
def pixnum():
    """
    This is a function that reads in a SDSS quasar file and produces a fits file that contains the RA, dec, and
    pixel number (in the HEALpix ring hierarchy) for a 4-side, 8-side, and 16-side pixel structure in a recarray.

    Args: none - takes in data from external file

    Returns: none - writes a fits file to wherever you would like the newly formatted data to go

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: format the SDSS data in a more useful way (easier to grab RA and dec and associates the pixel numbers
    #   with the coordinates), frame the SDSS coordinates in the HEALpix ring hierarchy
    #
    #
    # CALLING SEQUENCE: pixnum()
    #
    #
    # INPUTS: none - make sure the two paths go to somewhere that you would like them to go on your local device
    #
    # -------------------------------------------------------------
    #
    # MEK pull data from where quasar file is located (confirm location on your computer)
    datafile = "/Users/marykaldor/ASTR8080/runnoe/week3/HW1quasarfile.dat"
    data = np.genfromtxt(datafile, unpack=True, dtype=str)

    # MEK initialize lists to keep track of quasars
    coordlist = []
    ra = []
    dec = []

    # MEK pull RA and dec from the quasar file (the formatting is consistent so you can index through), use radians
    for q in data:
        coord = q[0:2] + " " + q[2:4] + " " + q[4:9] + " " + q[9:12] + " " + q[12:14] + " " + q[14:18]
        coord = SkyCoord(coord, unit=(u.hourangle, u.deg))
        coordlist.append(coord)
        ra.append(coord.ra.rad)
        dec.append(coord.dec.rad)

    # MEK convert lists to arrays
    ra = np.array(ra)
    dec = np.array(dec)

    # MEK convert RA and dec to phi and theta for the HEALpix ring hierarchy, pull pixel numbers for 4-, 8-, and 16-side
    phi = ra
    theta = (np.pi / 2) - dec
    pix4 = hp.ang2pix(4, theta, phi)
    pix8 = hp.ang2pix(8, theta, phi)
    pix16 = hp.ang2pix(16, theta, phi)

    # MEK create recarray and populate it
    type = np.dtype([('ra', 'f8'), ('dec', 'f8'), ('pixnum', 'f8', 3)])
    X = np.zeros(len(coordlist), dtype=type)
    X['ra'] = ra
    X['dec'] = dec
    X['pixnum'] = [px for px in zip(pix4, pix8, pix16)]

    # MEK write rec array to fits file somewhere local
    rec = fits.BinTableHDU(X)
    rec.writeto("/Users/marykaldor/ASTR8080/mary/week7/hw2_q3.fits")

###############################
def dense_pix():
    """
    This is a function that plots all SDSS quasars from the fits file in task 3 (pixnum() function) and highlights the
    5 densest pixels on top of those points

    Args: none - make sure you pull from the fits file wherever it is on your local device

    Returns: none - plots the quasars and densest pixels

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: visualize the SDSS sky and where the most quasars the Survey has are
    #
    #
    # CALLING SEQUENCE: dense_pix()
    #
    #
    # INPUTS: none - make sure you pull from the fits file that you wrote in task 3 (pixnum() function)
    #
    # -------------------------------------------------------------
    #
    # MEK open and read fits file that you just created
    f = fits.open("/Users/marykaldor/ASTR8080/mary/week7/hw2_q3.fits")
    hdr = f[0].header
    data = f[1].data

    # MEK pull RA, dec, and pixel number from it, convert from radians to degrees, plot
    ra_fits = data['ra']*180/np.pi
    dec_fits = data['dec']*180/np.pi
    numpix_fits = data['pixnum']
    plt.scatter(ra_fits, dec_fits, color='black', label="All Quasars")

    # MEK locate the pixel values for the 4-side HEALpix ring hierarchy for each quasar, find the 5 most common pixels,
    # plot them
    pix4 = [pix[0] for pix in numpix_fits]
    maxes = collections.Counter(pix4).most_common(5)
    colors = ["red", "yellow", "green", "blue", "purple"]
    labels = ["Densest", "2nd Densest", "3rd Densest", "4th Densest", "5th Dense"]
    i = 0
    for max in maxes:
        m = np.where(pix4==max[0])
        plt.scatter(ra_fits[m], dec_fits[m], color=colors[i], label=labels[i]+", pix "+str(math.trunc(max[0])))
        i += 1
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.title("Dec vs. RA of SDSS Quasars")
    plt.legend(fontsize=8)


# MAIN
###############################
###############################
if __name__ == '__main__':
    # MEK set the timer, start task 1
    time0 = time.time()

    # MEK calculate all 4 areas with same RA bounds and varying dec bounds (all dec bounds of same height)
    sq1 = square_area(0, 40, 0, 10)
    sq2 = square_area(0, 40, 20, 30)
    sq3 = square_area(0, 40, 40, 50)
    sq4 = square_area(0, 40, 60, 70)
    sq5 = square_area(0, 360, 0, 90)

    # MEK plot the 4 squares with their areas as their legend labels
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="aitoff")
    ax.grid(color="blue", linestyle="dashed", linewidth=2.5)
    xlab = ["14h", "16h", "18h", "20h", "22h", "0h", "2h", "4h", "6h", "8h", "10h"]
    ax.set_xticklabels(xlab, weight=800)
    ax.add_patch(Rectangle((sq1[1], sq1[3]), sq1[2] - sq1[1], sq1[4] - sq1[3], color="yellow", linewidth=1,
                           label="Area=" + str(round(sq1[0], 3))))
    ax.add_patch(Rectangle((sq2[1], sq2[3]), sq2[2] - sq2[1], sq2[4] - sq2[3], color="blue", linewidth=1,
                           label="Area=" + str(round(sq2[0], 3))))
    ax.add_patch(Rectangle((sq3[1], sq3[3]), sq3[2] - sq3[1], sq3[4] - sq3[3], color="green", linewidth=1,
                           label="Area=" + str(round(sq3[0], 3))))
    ax.add_patch(Rectangle((sq4[1], sq4[3]), sq4[2] - sq4[1], sq4[4] - sq4[3], color="orange", linewidth=1,
                           label="Area=" + str(round(sq4[0], 3))))
    plt.title("Same RA bounds with Equal (and Increasing) Dec Bounds", pad=20)
    plt.legend()

    # MEK time before opening plots (so that you can look for as long as you'd like!)
    time1 = time.time()
    print(" Time to complete task 1 first plot [mus]: {0:.3f}".format(1e6 * (time1 - time0)))

    plt.show()
    plt.close()
    
    # MEK reset time prior to starting second part of task 1 (so that viewing plots from part 1 doesn't impact timing)
    time2 = time.time()

    # MEK plot the full upper hemisphere - subtract pi from the RA bound so that it plots normally (Aitoff projection
    # goes from -pi to pi instead of 0 to 2pi, which is what the homework instructions said to input)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="aitoff")
    ax.grid(color="blue", linestyle="dashed", linewidth=2.5)
    xlab = ["14h", "16h", "18h", "20h", "22h", "0h", "2h", "4h", "6h", "8h", "10h"]
    ax.set_xticklabels(xlab, weight=800)
    ax.add_patch(Rectangle((sq5[1]-np.pi, sq5[3]), sq5[2] - sq5[1], sq5[4] - sq5[3], color="pink", linewidth=1,
                           label="Area=" + str(round(sq5[0], 3))))
    plt.title("Full Hemisphere Cap", pad=20)
    plt.legend()
    
    # MEK time before opening plots
    time3 = time.time()
    print(" Time to complete task 1 second plot [mus]: {0:.3f}".format(1e6 * (time3 - time2)))
    
    plt.show()
    plt.close()

    # MEK reset time, start task 2
    time4 = time.time()

    # MEK randomly populate the entire sky with 100,000 points, pull out RA and dec coordinates
    whole = random_pop(0,360,-90,90,100000)
    ra_rand = whole[0]
    dec_rand = whole[1]

    # MEK take a smaller section of this sky to see if the points/area ratio is consistent across smaller scales
    part = np.where((ra_rand > 0) & (ra_rand<100) & (dec_rand>0) & (dec_rand<100))

    # MEK plot all points and subgroup of points on top of them
    plt.plot(ra_rand, dec_rand, '.k', label="Full Sky")
    plt.plot(ra_rand[part], dec_rand[part], '.r', label="Small Part of Sky")
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.title("Randomly Populated Sky and Smaller Section")
    plt.legend()

    # MEK time before opening plots
    time5 = time.time()
    print(" Time to complete task 2 part 1 [mus]: {0:.3f}".format(1e6 * (time5 - time4)))

    plt.show()
    plt.close()

    # MEK reset time
    time6 = time.time()

    # Calculate ratio of points per unit area by dividing the total number of points (for the full sky and for the
    # smaller section) by the area of that section (either full sky area or area of smaller section)
    part_area = square_area(0,100,0,100)[0]
    whole_area = square_area(0,360,-90,90)[0]
    part_ratio = np.shape(part)[1]/part_area
    whole_ratio = np.shape(whole)[1]/whole_area
    print("Points per area in small section =", str(part_ratio))
    print("Points per area in whole sky =", str(whole_ratio))

    # MEK time before opening plots
    time7 = time.time()
    print(" Time to complete task 2 part 2 [mus]: {0:.3f}".format(1e6 * (time7 - time6)))
    # MEK they seem relatively similar, but won't be perfect until we get infinite points on the sky

    # MEK reset time, start task 3
    time8 = time.time()

    pixnum()

    # MEK time, start task 4
    time9 = time.time()
    print(" Time to complete task 3 [mus]: {0:.3f}".format(1e6 * (time9 - time8)))

    dense_pix()

    # MEK time
    time10 = time.time()
    print(" Time to complete task 4 [mus]: {0:.3f}".format(1e6 * (time10 - time9)))
    plt.show()
    plt.close()



