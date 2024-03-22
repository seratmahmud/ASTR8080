# M. Kaldor
# v1 2/23/2024
# ASTR 8080 HW3

# to import this from another directory:
# import sys
# sys.path.insert(0, '../hw/')
# import hw3
# from hw3 import function

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
import sys
sys.path.insert(0, '../week6/')
from spher_cap import ra_cap
from spher_cap import dec_cap
from spher_cap import gen_cap
sys.path.insert(0, '../hw/')
from hw2 import square_area
from hw2 import random_pop
import pymangle


# FUNCTIONS
###############################
###############################
def write_ply(caps, filename):
    """
    This is a function that writes a .ply file named "filename.ply" in the same directory that you call it in. This
    file is used later in the homework.

    Args: caps = [nx4 array]
        filename = [str]

    Returns: writes a file named "filename.ply"

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: Sort caps of the sky into a readable format for mangle to interpret
    #
    #
    # CALLING SEQUENCE: write_ply(rectangle, "hw3_mask")
    #
    #
    # INPUTS: caps = array containing any number of rows with 4 columns where each row denotes one cap on the sky,
    #       ** NOTE! The first 4 rows should be the input to the latitude-longitude rectangle that contains the survey
    #       for this homework, specifically.
    #       filename = what you would like to call the file that you write to (this will be referenced again in other
    #       functions)
    #
    # -------------------------------------------------------------
    #
    # MEK create and open a file in write mode
    f = open(filename+".ply", "w")

    # MEK set areas equal to an arbitrary 1 (area cannot be set for these caps ahead of time, as the HW questions
    # asks for the area to be calculated
    area = 1

    # MEK sort the caps into those that describe the latitude-longitude rectangle and those that describe circles on
    # the sky
    latlon = caps[0:4]
    circles = caps[4:]

    # MEK write any number of polygons depending on how many circles make up the survey
    f.write(f"{len(circles)} polygons\n")

    # MEK iterate through all circular observations, 4 caps for lat-lon rectangle and 1 cap for the observation
    for i in range(0, len(circles)):
        f.write(f" polygon {i+1} ( 5 caps, 1.0 weight, 0 pixel, {area} str):\n")
        # MEK increase number of spaces between each row to conform to mangle formatting
        for cap in latlon:
            f.write("  " + str(cap[0]) + " " + str(cap[1]) + " " + str(cap[2]) + " " + str(cap[3]) + "\n")
        f.write("  " + str(circles[i][0]) + " " + str(circles[i][1]) + " " + str(circles[i][2]) + " " + str(circles[i][3]) + "\n")
    # MEK close the written file
    f.close()

###############################
def latlon(ra1, ra2, dec1, dec2):
    """
    This is a function that produces the set of length=4 lists that will define the latitude-longitude rectangle
    for the survey

    Args: ra1 = [float]
        ra2 = [float]
        dec1 = [float]
        dec2 = [float]

    Returns: racap1 = [list]
        racap2 = [list]
        deccap1 = [list]
        deccap2 = [list]

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: Creating the borders of a latitude-longitude rectangle that properly indicate the inner area of the
    # 4 borders of the cap, which involves flipping the 4th (h) value in each list to be negative, which indicates to
    # mangle to reverse the part of the sky that it looks in.
    #
    #
    # CALLING SEQUENCE: latlon(10.25, 11.25, 30, 40)
    #
    #
    # INPUTS: ra1 = degrees of right ascension (lower limit of rectangle)
    #       ra2 = degrees of right ascension (upper limit of rectangle)
    #       dec1 = degrees of declination (lower limit of rectangle)
    #       dec2 = degrees of declination (upper limit of rectangle)
    #
    # -------------------------------------------------------------
    #
    # MEK call RA cap function that produces a length=4 list of x,y,z,h coordinates that outline the RA cap
    racap1 = ra_cap(ra1)
    racap2 = ra_cap(ra2)

    # MEK flip the sign on the last term of the list to indicate to mangle to use the area below/left of the line
    # instead of above/right of the line
    racap2[-1] = -1*racap2[-1]

    # MEK call dec cap function that produces a length=4 list of x,y,z,h coordinates that outline the dec cap
    deccap1 = dec_cap(dec1)
    deccap2 = dec_cap(dec2)

    # MEK flip the sign on the last term of the list to indicate to mangle to use the area below the line instead of
    # above the line
    deccap2[-1] = -1*deccap2[-1]

    # MEK return the 4 lists of length=4 that outline the area of the lat-lon rectangle
    return(racap1, racap2, deccap1, deccap2)

###############################
def random_area(ra1, ra2, dec1, dec2, n, show, filename):
    """
    This is a function that calculates the area of a non-rectangular area on the sky using a ratio that asserts the
    area of a region will be proportional to the number of points that fall within that region, if the points are
    randomly distributed by area on the sky.

    Args: ra1 = [float]
        ra2 = [float]
        dec1 = [float]
        dec2 = [float]

    Returns: area = [float]

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: calculate area of a non-rectangular area on the sky (in this case, the observational survey from HW3)
    #
    #
    # CALLING SEQUENCE: random_area(151, 170, 31, 39, 1000, False, "hw3_mask")
    #
    #
    # INPUTS: ra1 = degrees of right ascension (lower limit of rectangle)
    #       ra2 = degrees of right ascension (upper limit of rectangle)
    #       dec1 = degrees of declination (lower limit of rectangle)
    #       dec2 = degrees of declination (upper limit of rectangle)
    #       n = number of points used to calculate area
    #       show = whether or not you would like to see a plot and some informative print statements
    #       filename = what you have previously called your .ply file that you wrote
    #
    #      ** These RA and dec bounds must be equal to or larger than the borders of your survey to ensure that the
    #      entire area of the survey is calculated via this method
    #
    # -------------------------------------------------------------
    #
    # MEK initialize mask that describes survey observations that you wrote earlier
    mask = pymangle.Mangle(filename+".ply")

    # MEK randomly populate some rectangle that is larger than or equal to the lat-lon rectangle to ensure you cover
    # the entirety of the survey mask
    big = random_pop(ra1, ra2, dec1, dec2, n)

    # MEK calculate the number of points that fall within the mask
    good = mask.contains(big[0], big[1])
    inmask = np.where(good == True)
    num_in = np.shape(inmask)[1]

    # MEK show the points that fall both in the mask (magenta) and outside of the mask (black) if you'd like
    if show == True:
        print("Mask", mask)
        print("Number in mask", num_in)
        plt.scatter(big[0], big[1], color="black")
        plt.scatter(big[0][inmask], big[1][inmask], color="magenta")
        plt.show()

    # MEK calculate mask area as a ratio of (total area)/(total points) = (mask area)/(points in mask)
    area = num_in*square_area(ra1, ra2, dec1, dec2)[0]/n

    return area

###############################
def pick_points(ra1, ra2, dec1, dec2, show):
    """
    This is a function iterates through a number of points randomly distributed in area over the sky to calculate the
    area of a survey mask of any shape by comparing the ratio of points within the mask to those outside of the mask
    but within a region of known area. The goal of this function is to produce areas that are consistent with one
    another across calculations for 100 iterations of area calculations.

    Args: ra1 = [float]
        ra2 = [float]
        dec1 = [float]
        dec2 = [float]
        show = [Boolean]

    Returns: n = [int]

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: decide a reasonable number of points used to calculate the area of a survey mask without using too many
    # but ensuring that the answer is consistent enough across 100 iterations of area calculations
    #
    #
    # CALLING SEQUENCE: pick_points(151, 170, 31, 39, True)
    #
    #
    # INPUTS: ra1 = degrees of right ascension (lower limit of rectangle)
    #       ra2 = degrees of right ascension (upper limit of rectangle)
    #       dec1 = degrees of declination (lower limit of rectangle)
    #       dec2 = degrees of declination (upper limit of rectangle)
    #       show = whether or not you would like to see a plot
    #
    #       ** These RA and dec bounds must be equal to or larger than the borders of your survey to ensure that the
    #       entire area of the survey is calculated via this method
    #
    # -------------------------------------------------------------
    #
    # MEK start with 100 points to calculate the area of the mask
    n = 100

    # MEK set some large bound for the number of points used to calculate area (we will never reach this many points)
    while n<1e20:

        # MEK calculate the area 100 times for the same number of points, and calculate the standard deviation of these
        # 100 areas
        i = 0
        areas = []
        while i<100:

            # MEK calculate the area of the mask with n number of points
            area = random_area(ra1, ra2, dec1, dec2, n, False, "hw3_mask")
            areas.append(area)
            i += 1
        areas = np.array(areas)
        stddev = np.std(areas)

        # MEK show the trend of standard deviations in area against number of points used to calculate area
        plt.scatter(n, stddev)
        if stddev<0.5:

            # MEK if the standard deviation is below the threshold of 0.5 square degrees, return the number of points
            # that created this area
            plt.plot([0, n], [0.5, 0.5], color="black", label="st dev=0.5")

            # MEK plot the trend in standard deviations if you want!
            if show==True:
                plt.title("Standard Deviations of Areas vs. Number of Points")
                plt.xlabel("Number of Points")
                plt.ylabel("Standard Deviation of Areas")
                plt.legend()
                plt.show()
            return n

        # MEK if the standard deviation is still too large, add 100 points and try again
        else:
            n += 100


###############################
def quasar_mask(area, filename):
    """
    This is a function that reads in quasar locations from SDSS observations, plots them on the sky, and highlights
    which quasars are inside of the mask that describes the survey capabilities of a supposed measurement described
    in this homework. The function creates a plot that labels the area of the mask and the number density of quasars
    that lie within the mask.

    Args: area = [float]
        filename = [str]

    Returns: saves and shows a plot as described above

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: Showing where all the quasars in an SDSS data file are located in the sky and noting which of these
    # quasars lie within the mask of the survey described in this homework
    #
    #
    # CALLING SEQUENCE: quasar_mask(area, "hw3_mask")
    #
    #
    # INPUTS: area = area of the survey mask to use in the legend
    #       filename = name of file where the mask you have written is stored
    #
    # -------------------------------------------------------------
    #
    # MEK read in coordinates in hours-mintues-seconds (right ascension) and degrees-arcminutes-arcseconds (declination)
    # for all quasars in the SDSS data file
    datafile = "../../runnoe/week7/HW3quasarfile.dat"
    data = np.genfromtxt(datafile, unpack=True, dtype=str)
    coordlist = []
    for q in data:
        coord = q[0:2] + " " + q[2:4] + " " + q[4:9] + " " + q[9:12] + " " + q[12:14] + " " + q[14:18]
        coordlist.append(coord)
    coords = np.array(coordlist)

    # MEK create a SkyCoord object of all quasars in the data file, plot them
    quasars = SkyCoord(coords, unit=(u.hourangle, u.deg))
    plt.scatter(quasars.ra, quasars.dec, color="green", label="All Quasars")

    # MEK create a mask to mark which quasars fall inside and outside of the survey area
    mask = pymangle.Mangle(filename+".ply")
    good_q = mask.contains(quasars.ra, quasars.dec)
    inmask_q = np.where(good_q == True)

    # MEK note how many quasars fall inside the mask
    num_in_q = np.shape(inmask_q)[1]

    # MEK plot quasars that are in the mask in a different color, label with mask area and number density of quasars
    # in the mask
    plt.scatter(quasars.ra[inmask_q], quasars.dec[inmask_q], color="yellow", label="mask area="+
                        str(round(area, 3))+" sq deg, num. dens.="+str(round(num_in_q/area, 3))+" quas/sq deg")

    # MEK add in labels to plot
    plt.legend()
    plt.xlabel("RA")
    plt.ylabel("Dec")
    plt.title("SDSS Quasars in Mask")

    # MEK save and show the plot
    plt.savefig("quasars_mask.png")
    plt.show()


# MAIN
###############################
###############################
if __name__ == '__main__':
    # MEK set the timer
    time0 = time.time()

    # MEK create array to describe the lat-lon rectangle as one variable
    rectangle = np.array(latlon(10.25, 11.25, 30, 40))

    # MEK time
    time1 = time.time()
    print(" Time to create lat-lon rectangle [s]: {0:.3f}".format(time1 - time0))

    # MEK create 4 length=4 lists that describe the observation areas as general circular caps on the sky
    gencap1 = gen_cap(155, 34, 2)
    gencap2 = gen_cap(159, 36, 2)
    gencap3 = gen_cap(163, 34, 2)
    gencap4 = gen_cap(167, 36, 2)
    # MEK create array to describe all caps of observations as one variable
    allcaps = np.array([gencap1, gencap2, gencap3, gencap4])

    # MEK time
    time2 = time.time()
    print(" Time to create caps for survey [s]: {0:.3f}".format(time2 - time1))

    # MEK combine the lat-lon rectangle with all 4 observation caps to create the area covered by the survey as one
    # variable
    rectangle = np.concatenate((rectangle, allcaps), axis=0)
    # MEK write this information to a .ply file named after HW3
    write_ply(rectangle, "hw3_mask")

    # MEK time
    time3 = time.time()
    print(" Time to write .ply file [s]: {0:.3f}".format(time3 - time2))

    # MEK pick the number of points that properly calculate the area of the sky mask to within 0.5 degrees
    points = pick_points(151, 170, 31, 39, True)

    # MEK time
    time4 = time.time()
    print(" Time to pick number of points [s]: {0:.3f}".format(time4 - time3))

    # MEK calculate the area that that number of points gives you
    area = random_area(151, 170, 31, 39, points, False, "hw3_mask")
    print(area, "square degrees, calculated using", points, "points")

    # MEK time
    time5 = time.time()
    print(" Time to calculate area for mask [s]: {0:.3f}".format(time5 - time4))

    # MEK print the SDSS quasars that fall within the mask in a specific color, and produce a plot that visualizes
    # where they are on the sky, what the area of the mask is, and what the number density of quasars per square
    # degree are in that mask on the sky
    quasar_mask(area, "hw3_mask")

    # MEK time
    time6 = time.time()
    print(" Time to map quasars [s]: {0:.3f}".format(time6 - time5))

    # MEK time
    time7 = time.time()
    print(" Total time to run HW3 (includes time for plots being open) [s]: {0:.3f}".format(time7 - time0))
