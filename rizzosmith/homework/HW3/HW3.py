# M. Rizzo Smith
# v1 2/23/24
# Homework 3
# ASTR 8080


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
import pymangle
from numpy.random import random
import time

sys.path.insert(0, '../../week6/')
import spherical_caps
from spherical_caps import ra_cap, dec_cap, gen_cap

# FUNCTIONS
###############################
###############################

def lat_lon(rectangle):
    """
    This is a function to generate the 4 caps required to define a lat-lon rectangle. 
    The function takes in the corners of the lat-lon rectangle in degrees as an array, and returns the 4 caps formatted to be written to a .ply file

    Args:
        rectangle (numpy array): RA and Dec bounds of lat long rectangle. RA is passed as 00m00h00s and dec is passed as decimal degrees.

    Returns:
        arr_caps (numpy arary): Array of strings of the spherical caps that define the lat-lon rectangle to be written ta a .ply file
    """
#-----------------------------------------------------------------
#+
# PURPOSE:
#   Produce the 4 caps required to define a lat-lon rectangle
#
# CALLING SEQUENCE:
#   rect_caps = lat_lon(['10h15m00s', '11h15m00s', 30, 40])
#
#
# INPUTS:
#   rectangle - Numpy Array of RA and DEC bounds to lat-lon rectangle
#-
#------------------------------------------------------------------
    # MRS Reformatting RA from HHMMSS into degrees and pulling ra dec bounds from rectangle corners
    ra1 = float(rectangle[0][:2]) + float(rectangle[0][3:5])/60 + float(rectangle[0][6:8])/3600
    ra2 = float(rectangle[1][:2]) + float(rectangle[1][3:5])/60 + float(rectangle[1][6:8])/3600
    dec1 = rectangle[2]
    dec2 = rectangle[3]
    caps = []
    
    # MRS Create spherical caps for each ra dec bound
    ra1_cap = (ra_cap(ra1))
    ra2_cap = (ra_cap(ra2))
    dec1_cap = (dec_cap(dec1))
    dec2_cap = (dec_cap(dec2))

    # MRS Flip last element of RA max and Dec max to define lat-lon rectangle
    ra2_cap[-1] *= -1
    dec2_cap[-1] *= -1
    caps.append(ra1_cap)
    caps.append(ra2_cap)
    caps.append(dec1_cap)
    caps.append(dec2_cap)
    arr_caps = np.array(caps)

    # MRS Reformatting to write nicely to ply file
    caps_str= str(np.array(caps))
    caps_str = caps_str.replace(',', ' ')
    caps_str = caps_str.replace('[', ' ')
    caps_str = caps_str.replace(']', '')
    
    return arr_caps

def write_ply(caps):
    """
    Function to take in an arbitray number of caps. Assumes the first 4 caps correspond to a lat-lon rectangle bounded region, and the follwoing caps are individual spherical caps.

    Args:
        caps (numpy array) - Array of strings that correspond to caps that we wish to write to a ply file.

    Returns:
        Function does not return any variables but creates a .ply file called 'aal_caps.ply' in the current directory
    """
#-----------------------------------------------------------------------
#+
# PURPOSE:
#   Take in an array of caps to write to a single .ply file
#
# CALLING SEQUENCE:
#   write_ply(caps)
#
# INPUTS:
#   caps - numpy array of caps written as strings
#-
#-----------------------------------------------------------------------
    # MRS Create the ply file, specify lat lon as first 4 caps
    f = open('all_caps.ply', 'w')
    latlon = caps[:4]
    circles = caps[4:]

    # MRS Write PLY file assumes 1 for area
    f.write(f'{len(circles)} polygons\n')
    for i in range(0, len(circles)):
        f.write(f'polygon {i+1} ( {len(latlon)+1} caps, 1 weight, 0 pixel, 1 str):\n')
        for cap in latlon:
            f.write(f' {cap[0]:.5f} {cap[1]:.5f} {cap[2]:.5f} {cap[3]:.5f}\n')
        f.write(f' {circles[i][0]:.5f} {circles[i][1]:.5f} {circles[i][2]:.5f} {circles[i][3]:.5f}\n')
    f.close()
    
    return

def rand_pop(rectangle, points):
    """
    A funciton to populate a given lat-lon rectangle with random points equally in area, and return the ra and dec of those randomly generated points. 
    Takes in RA in the format '00h00m00s' and decimal degrees for DEC
    
    Args:
        rectangle (numpy array): The RA and DEC bounds for the lat-lon rectangle in which we want to generate points
        points (integer): Number of random points to generate
    
    Returns:
        ra_sky (array): RA values of randomly generated points
        dec_sky (array: DEC Values of randomly generated points)
    """
#----------------------------------------------------------------
#+
# PURPOSE:
#   Randomly populate a lat-lon rectangle in equal area with a specified number of points
#
# CALLING SEQUENCE:
#   ra, dec = rand_pop(rectangle, 1000)
#
# INPUTS:
#   rectangle - array containg ra and dec bounds of lat-lon rectanlge [RA_Min, RA_max, DEC_min, DEC_max]
#   points - integer of number of points to generate
#-
#----------------------------------------------------------------

    # MRS Reformatting RA and Dec from HHMMSS to degree
    rmin = (float(rectangle[0][:2]) + float(rectangle[0][3:5])/60 + float(rectangle[0][6:8])/3600) * 15
    rmax = (float(rectangle[1][:2]) + float(rectangle[1][3:5])/60 + float(rectangle[1][6:8])/3600) * 15
    dmin = rectangle[2]
    dmax=rectangle[3]

    # MRS Generate RA and Dec randomly in area for a given lat-lon rectangle
    ra_sky = (random(points)*(rmax - rmin)) + (rmin)
    dec_sky = (180./np.pi) * (np.arcsin(random(points)*(np.sin(dmax *(np.pi/180))-(np.sin(dmin *(np.pi/180))))+np.sin(dmin*(np.pi/180))))
    
    return ra_sky, dec_sky

def get_area(rectangle):
    """
    This is a function to generate the area in square degrees of a specified lat-lon rectangle

    Args: 
        rectangle (array): RA and Dec bounds of lat-lon rectangle RA in format '00h00m00s' and Dec in decimal degrees 

    Returns:
        area (float): Area of rectangle in square degrees.
    """
#----------------------------------------------------------------
#+
# PURPOSE:
#   Calcualte the area of a lat-lon rectangle
#
# CALLING SEQUENCE:
#   area = get_area(rectanlge)
#
# INPUTS:
#   rectangle - array of RA and DEC bounds [RA_Min, RA_max, DEC_min, DEC_max]
#-
#----------------------------------------------------------------

    # MRS Reformatting RA and Dec from HHMMSS to degree
    rmin = (float(rectangle[0][:2]) + float(rectangle[0][3:5])/60 + float(rectangle[0][6:8])/3600) * 15
    rmax = (float(rectangle[1][:2]) + float(rectangle[1][3:5])/60 + float(rectangle[1][6:8])/3600) * 15
    dmin = rectangle[2]
    dmax=rectangle[3]
    
    # MRS Calculate width and height of triangle to get area
    h = np.sin(dmax * (np.pi/180)) - np.sin(dmin*(np.pi/180))
    w = ((rmax) * (np.pi/180)) - ((rmin) * (np.pi/180))
    area = h*w *(180/np.pi)**2

    return area
    

def pop_rect(rectangle, ply_file, quasar_file, bool1, bool2):
    """
    Function to numerically calculate the area of a given mask. If passed bool1=0 it will not do the iterating to acheieve the desired 0.5 square degree precesion. By default this function will plot the quasar file in the directory, determine the number density of quasars in the mask and plot it. 

    Args:
        rectangle (array): ra and dec bounds for the lat-lon rectangle of which we want to generate random points in.
        ply_file (string): String pointing to the polygon file of our observing footprint
        quasar_file (string): String pointing to the file holding the positons of the ra and dec of quasars for HW question 3. 
        bool1 (integer): Default value is 1. Boolean statement for if we want to iterate to determine mask area.
        bool2 (integer): Default value is 1, Boolean statement if we want to utilize the quasar file to plot positions

    Returns:
        Does not return a variable, but displays plot of quasar positions.
    """
#-----------------------------------------------------------------
#+
# PURPOSE: 
#   Determine the area of a mask to a precision of 0.5deg^2
#   Plot the positions of quasars that fall inside of our mask, and determine the number density of quasars in mask.
#
# CALLING SEQUENCE:
#   With area precision and quasar plotting:
#       pop_rect(rectangle, 'all_caps.ply', 'Quasar_file.fits', 1, 1)
#   With just plotting
#       pop_rect(rectangle, 'all_caps.ply', 'Quasar_file.fits', 0, 1)
#
# RETURNS:
#   Doesn't return variables but generate Quasar position plots.
#-
#----------------------------------------------------------------
    # MRS Timestamp for this function starting
    time_int = time.time()

    # MRS Initialize mask from PLY file and calculate area of lat-lon rectangles
    m = pymangle.Mangle(ply_file)
    area_rec = get_area(rectangle)
    
    # MRS Let the user know that were doing something and it may take a sec
    print(r'Iterating to meet desired precision of 0.5 square degrees...')
    # MRS Initializing lists for determing number of points requried for desired precision
    areas_tmp = []
    num_points = 100

    # MRS Loop through and calculate the area 100 times for a given number of random points. If precision is met exit, if not do it again
    # MRS Setting large number of points to ensure I dont time out early, def won't need to loop that far.
    if bool1 == 1:
        while num_points < 10000000:
            for i in range(0, 100):
                ra_sky, dec_sky = rand_pop(rectangle, num_points)
                in_mask = m.contains(ra_sky, dec_sky)
                area_mask1 = (area_rec / num_points) * len(ra_sky[in_mask])
                areas_tmp.append(area_mask1)
                i += 1
            areas_tmp = np.array(areas_tmp)
            if (np.std(areas_tmp) < 0.5):
                print(f'Area Mask: {np.mean(areas_tmp):.2f} square degrees', f'Required Num Points: {num_points}')
                break
            else:
                num_points += 100
                i = 0
                areas_tmp = []    
    time_fin = time.time()
    print('Time to reach precision of 0.5 square degree [s]: {0:.3f}'.format((time_fin-time_int)))


    # MRS Plot Quasar positions, determine number density of quasars in mask, generate png of Quasars
    
    if bool2 == 1:
        print('Plotting Quasar positions and cross-matching with polygon mask...')
        time_plot = time.time()

        all_pos = []
        with open(quasar_file, 'r') as f:
            for line in f:
                pos_str = line
                # MRS Reformat the string from the data file into a nicer string to pass to SkyCoord
                pos_str = f'{pos_str[:2]} {pos_str[2:4]} {pos_str[4:9]} {pos_str[9:12]} {pos_str[12:14]} {pos_str[14:]}'

                # MRS Add the quasar position string to the all_pos list
                all_pos.append(pos_str)
        positions = SkyCoord(np.array(all_pos), unit=(u.hourangle, u.deg))

        # MRS Store indices of quasars in the mask
        quasar_match = m.contains(positions.ra.degree, positions.dec.degree)
        num_dens = len(positions.ra[quasar_match]) / area_mask1
        
        # MRS Timestamp for plotting
        time_plot_fin = time.time()
        print('Time to crossmatch quasar positions and plot [s]: {0:.3f}'.format((time_plot_fin-time_plot)))
        
        # MRS Plot all of the points and overplot the points which fall within the mask
        # MRS Couple extra lines in here to place the mask area and number density text onto the plot
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot(111)
        props = dict(facecolor='white', alpha=0.7, edgecolor='none')
        plt.text(0.05, 0.95, f'Mask Area: {area_mask1:.2f} deg squared', transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)
        plt.text(0.05, 0.9, f'Number Density: {num_dens:.2f} per deg squared', transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)
        plt.plot(positions.ra.degree, positions.dec.degree, 'o', color='lightblue')
        plt.plot(positions.ra.degree[quasar_match], positions.dec.degree[quasar_match], 'o', color = 'tomato', label='In mask')
        plt.xlabel('RA')
        plt.ylabel('Dec')
        plt.title('Cross-matching Quasars with Survey Field')
        plt.savefig('quasar_pos.png', dpi=300)
        plt.show()

    return  



###############################
###############################
if __name__ == '__main__':
    
    # MRS Keep track of time
    time0 = time.time()
    
    # MRS Generate lat-lon caps, plate caps, and write to Poly gon file
    rect = ['10h15m00s', '11h15m00s', 30, 40]
    rect_caps = lat_lon(rect)
    spher_caps = np.array([gen_cap(155, 34, 2), gen_cap(159, 36, 2), gen_cap(163, 34, 2), gen_cap(167,36, 2)])
    
    # MRS Format all the caps into one long array
    all_caps = np.concatenate((rect_caps, spher_caps), axis=0)
    write_ply(all_caps)
    
    # MRS First step timestamp
    time1 = time.time()
    print("Time to generate all caps and write to polygon file [ms]: {0:.3f}".format(1e3*(time1-time0)))
    
    # MRS Default values
    prec = 1
    qua_plot = 1

    # MRS Check if arguments are provided
    if len(sys.argv) > 1:
        if sys.argv[1] == 'mcm_off':
            prec = 0
        if sys.argv[1] == 'plot_off':
            qua_plot = 0
        if len(sys.argv) > 2:
            if sys.argv[2] == 'plot_off':
                qua_plot = 0
    
    # MRS Run the mask area calculation and quasar cross-match plotting 
    pop_rect(rect, 'all_caps.ply', '../../../runnoe/week7/HW3quasarfile.dat', prec, qua_plot)
