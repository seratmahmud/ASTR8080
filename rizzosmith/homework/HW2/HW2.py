# M. Rizzo Smith
# v1 2/9/24
# ASTR 8080
# HW 2 - Areas on a sphere


# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
import astropy 
from matplotlib import rc
from astropy import units as u
from astropy.coordinates import SkyCoord
from numpy.random import random
from matplotlib.patches import Rectangle
import healpy
from astropy.io import fits
from collections import Counter

# FUNCTIONS
###############################
###############################

def get_area(rmin, rmax, dmin, dmax, boo):
    """
    This is a function to calculate the area on the sky of a given lat-lon rectangle. 
    This function takes in the corners of a lat-lon rectangle in degrees and calulates its area in square degrees
    If passes the last argument it will plot 4 rectangles of increasing maximum declination and compare their areas (This is to satisfy task 1)

    Args:
        rmin (float): minimum RA value (degrees)
        rmax (float): maxmium RA value (degrees)
        dmin (float): minimum Dec value (degrees)
        dmax (float): maxmium Dec value (degrees)
        boo (integer): If zero plot 4 rectangles

    Returns:
        area (float): The area of the rectangle in square degrees
    """
#--------------------------------------------------------------------
#+
# PURPOSE:
#   Calcualte the area of a given lat-lon rectangle on the sky
#   Provide a plot of 4 increasingly large rectangles and compare areas
#
# CALLING SEQUENCE:
#   With Plot:
#       get_area(rmin, rmax, dmin, dmax, 1)
#   Without Plot:
#       get_area(rmin, rmax, dmin, dmax, 0)
#
# INPUTS:
#   rmin, rmax, dmin, dmax - Rectangle Corner Bounds
#   boo - 1 or 0
#-
#-------------------------------------------------------------------
    # MRS Calculate the height of the triangle 
    h = np.sin(dmax * (np.pi/180)) - np.sin(dmin*(np.pi/180))
    # MRS Calculate the width of the triangle
    w = (rmax * (np.pi/180)) - (rmin * (np.pi/180))
    # MRS Combine both for the area and convert to square degrees
    area = h*w * (180/np.pi)**2

    # MRS If flagged to plot, we make rectangles
    if boo == 1:
        # MRS Defining a corner array to pass to matplotlib.patches.Rectangle
        corners = [rmin, rmax, dmin, dmax]
        # MRS Initializing a variable to increase dmax for larger rectangles 
        j = 0

        # MRS Setup plot
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111, projection='aitoff')
        colors = plt.cm.tab10(range(len(corners)))
        for i in range(0, 4):
            # MRS Bottom corner to pass to Rectangle
            xy = (corners[0]*(np.pi/180.), (corners[2]+j)*(np.pi/180.))
            # MRS Width of the rectangle (RA diff)
            width = (corners[1]-corners[0])*(np.pi/180.)
            # MRS Height of rectangle (Dec Diff)
            height = ((corners[3]+j)-(corners[2]+j))*(np.pi/180.)
            # MRS Create the rectangle object and calculat the area
            rectangle = Rectangle(xy, width, height, edgecolor='black', facecolor=colors[i], label=f'Area: {get_area(corners[0], corners[1], corners[2]+j, corners[3]+j,0):.2f} $deg^{2}$')
            # MRS Add the patch to the plot
            ax.add_patch(rectangle)
            # MRS Make the next rectangle 10degrees taller in Dec
            j+=20
        # MRS Make plot look prettier          
        xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
        ax.set_xticklabels(xlab, weight=800)
        ax.grid(color='black', linestyle='solid', linewidth=1.5)
        ax.set_xlabel('Right Ascension')
        ax.set_ylabel('Declination')
        plt.title('Area Comparison of Increasing \n Lat-Lon Rectangles')
        plt.legend(loc = 'lower right')
        plt.show()
    # MRS If not passed 1 dont plot, and just return the area
    else:
        pass

    return area  

def rand_pop(rmin, rmax, dmin, dmax):
    """
    This function will randomly populate a lat-lon rectangle on the sky in equal area space. 
    The funciton will also provide the ratio of the rectangle area to sky area, and ratio of points in the rectangle proportional to number of points on the whole sky.

    Args:
        rmin (float): minimum RA value (degrees)
        rmax (float): maxmium RA value (degrees)
        dmin (float): minimum Dec value (degrees)
        dmax (float): maxmium Dec value (degrees)
    Returns:
        ra (numpy array): The RA values of all points in rectangle
        dec (numpy array): The Dec values of all points in rectangle
    """
#------------------------------------------------------------
#+
# PURPOSE:
#   Populate a given lat-lon rectangle randomly in equal area
#   Calculate and return the ratio of areas and ratio of points with respect to the whole sky
#
# CALLING SEQUENCE:
#   ra, dec = rand_pop(rmin, rmax, dmin, dmax)
#
# INPUTS:
#   rmin, rmax, dmin, dmax - Rectangle Corner Bounds
#-
#-------------------------------------------------------------
    # MRS Generate points on the whole sky to draw from
    ra_sky = (random(100000)*360)
    dec_sky = (180./np.pi) * (np.arcsin(1-random(100000)*2))
    
    # MRS Pull the RA and Dec that lie within the lat-lon rectangle
    inds = np.where((ra_sky >= rmin) & (ra_sky <= rmax) & (dec_sky >= dmin) & (dec_sky <= dmax))
   
    # MRS Save those RA and Dec values
    ra_rec = ra_sky[inds]
    dec_rec = dec_sky[inds]
   
    # MRS Calculate the ratio of rectangle area to the total sky area
    area = get_area(rmin, rmax, dmin, dmax, 0) 
    a_ratio = area / (4 * np.pi * (180/np.pi)**2)

    # MRS Calculate the ratio of points in the rectangle to total points in sky
    point_ratio = len(ra_rec) / 100000
    
    # MRS Print those ratios
    print(f'Area Ratio: {a_ratio:.4f}', f'Num Points Ratio: {point_ratio:.4f}')
    
    return ra_rec, dec_rec

def get_pix(quasar_file):
    """
    This function will read in a data file of ra and dec for multiple objects, then store the information in the form of a FITS file.
    The data will be formated into [ra, dec, pixnum] where pixnum is a 3-array saving the pixel numbers of each object in the Nside = 4, 8, 16 HEALpix hirerachy.
    Args:
        quasar_file (data file): Data file assumes the RA and Dec of objects are already formatted to HHMMSS.SS+DDMMSS.SS strings.
    Returns:
        N/A Creates a FITS file of the data
    """
#----------------------------------------------------
#+
# PURPOSE:
#   Create a FITS file from the data file passed to the function. Determine pixel positions of objects in the Nside = 4, 8, 16 HEALpix hierarchies. 
#
# CALLING SEQUENCE:
#   get_pix(data_file)
#
# INPUTS:
#   data_file - a file of ra and dec positons of a given set of objects.
#
#
#-
#----------------------------------------------------
    # MRS Initialize list to store each quasar
    all_pos = []
    with open(quasar_file, 'r') as f:
        for line in f:
            pos_str = line
            # MRS Reformat the string from the data file into a nicer string to pass to SkyCoord
            pos_str = f'{pos_str[:2]} {pos_str[2:4]} {pos_str[4:9]} {pos_str[9:12]} {pos_str[12:14]} {pos_str[14:]}'
            
            # MRS Add the quasar position string to the all_pos list
            all_pos.append(pos_str)
    # MRS Create a SkyCoord object for all of the Quasar positions
    positions = SkyCoord(np.array(all_pos), unit=(u.hourangle, u.deg))
    # MRS Define Phi to be passed to the HEALpix functions
    phi = positions.ra.degree * np.pi/180.
    # MRS Define Theta to be passed to the HEALpix functions
    theta = (np.pi / 2) - (positions.dec.degree * (np.pi/180.))
    
    # MRS Return the pixel indices for the Nside = 4, 8, 16 HEALpix hierarchies
    pix4 = healpy.ang2pix(4, theta, phi)
    pix8 = healpy.ang2pix(8, theta, phi)
    pix16 = healpy.ang2pix(16, theta, phi)
    
    # MRS Format a rec array in [ra, de, pixnum] with pixnum being a 3-array [4pix, 8pix, 16pix]
    tp = np.dtype([('ra', 'f8'), ('dec', 'f8'), ('pixnum', 'i8', (3))])
    # MRS Create the empty array of appropriate size to populate later
    all_objs = np.zeros(len(positions), dtype=tp)
    
    # MRS Populate each quasar into the rec array
    all_objs['ra'] = phi * 180./np.pi
    all_objs['dec'] = positions.dec.degree
    all_objs['pixnum'] = [px for px in zip(pix4, pix8, pix16)]
    
    # MRS convert the all_objs rec array into a FITS file format
    hdu = fits.BinTableHDU(all_objs)
    # MRS Write and save the FITS file
    hdu.writeto('Quasars.fits', overwrite=True)
    
    return

def plot_fits(fits_file):
    """
    This function will open and read in a FITS file of Quasar data.
    The function will plot the positions of all the quasars in the data file, and then overplot the quasars which occur in the top 5 most densly populated pixels in the Nside = 4 HEALpix hierarchy. 
    Args:
        fits_file (FITS file): The file containing quasar data.
    Returns:
        N/A Just plots the data in the FITS file
    """
#------------------------------------------------------
#+
# PURPOSE:
#   Open a given FITS file, read in the data, and plot the positions of the objects in the FITS file. Will also overplot the objects in the top 5 densest pixels in the Nside = 4 HEALpix hierarchy
#
# CALLING SEQUENCE:
#   plot_fits(fits_file)
#
# INPUTS:
#   fits_fle - FITS file of data
#-
#------------------------------------------------------
    # MRS Open up the FITS file
    hdul = fits.open(fits_file)
    # MRS Pull out the data
    data = hdul[1].data
    # MRS Take out the Nside=4 HEALpix hierarchy values from the pixnum 3-array
    pix4 = [px[0] for px in data['pixnum']]
    # MRS Initialize a counter object to identify top 5 populated pixels
    counter = Counter(pix4)
    # MRS Store the indeces of the 5 most populated pixels
    five_common = counter.most_common(5)
    colors = ['red', 'yellow', 'green', 'blue', 'purple']
    labels= ['1st Densest', '2nd', '3rd', '4th', '5th']
    # MRS Plot the Quasar positions
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, projection='aitoff')
    ax.grid(color='black', linestyle='solid', linewidth=1.5)
    ax.scatter(data['ra']*(np.pi/180.) - (np.pi), data['dec']*(np.pi/180.), s=10, color='black', alpha=0.9, label='All Quasars')
    # MRS Plot a dummy line to add the label to the legend outside of the loop
    # MRS Loop through and plot the top 5 densest pixels ontop of all the data
    for i in range (0, len(five_common)):
        ind = np.where(pix4 == five_common[i][0])
        ax.scatter((data['ra'][ind]) *(np.pi/180.)-(np.pi), (data['dec'][ind])*(np.pi/180.), marker='v', s=50, color=colors[i], label=labels[i])
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
    ax.set_xticklabels(xlab, weight=800)
    plt.legend()
    plt.show()

    #for i in range(0, len(five_common[0])):


    return

# MAIN
###############################
###############################
if __name__ == '__main__':
    # MRS Calling fucntion to complete Problem 1
    # MRS Create a rectangle boudned by given corners, calcualtes areas, and plots 3 more rectangles with increasing declination to compare area
    get_area(0, 40, 0, 10, 1)

    # MRS Calling Function to complete Problem 2
    # MRS Randomly populate rectanlge field, calculate its area.
    # MRS Will also print the ratio of the area to the whole sky and ratio of number of points as compared to the whole sky.
    ra, dec = rand_pop(0, 20, -10, 10)

    # MRS Calling function to complete Problem 3
    # MRS Writes out a FITS file of the input quasar data.
    get_pix('HW1quasarfile.dat')
    
    # MRS Calling function to complete Problem 4
    # MRS Reads in and plots the data in the FITS file we made.
    plot_fits('Quasars.fits')


