# S. Saad
# v1 4/9/2024
# ASTR 8080 HW4


#Import
############
##########
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from time import time
from astropy.io import fits
from sdss_sweep_data_index import sdss_sweep_data_index
import matplotlib.pyplot as plt
from create_test_data import fetch_sweep_objects
import warnings
warnings.filterwarnings("ignore")

#Functions
############
##########

def plot_get_colorcuts():
    """
    NAME: plot_get_colorcuts

    PURPOSE: Plots colors to get the color cuts for quasars

    REVISION HISTORY:

    v1.0: version Serat Saad, April 9, 2024 to submit as a homework for
    ASTR 8080 class at Vanderbilt University
    """

    # SS Coordinate and radius of the trageted region in the sky
    ra = 180
    dec = 30
    radius = 3

    # SS Either call fetch_sweep_objects to get the objects, or 
    # call a file that already contains the sweep object
     
    objs = fetch_sweep_objects(ra, dec, radius)
    #hdul_objs = fits.open('sources-ra180-dec30-rad3.fits')
    #bjs = hdul_objs[1].data
    #hdul_objs.close()

    # SS getting the coordinates of the sources
    ra_obj = objs["RA"]
    dec_obj = objs["DEC"]
   
    # SS Reading data from the file that contains the quasar sources
    hdul = fits.open("/home/saadsm/ASTR8080/runnoe/week11/qsos-ra180-dec30-rad3.fits")
    data = hdul[1].data
    ra_cat = data["RA"]
    dec_cat = data["DEC"]

    # SS Matching the qusar sources with the other sources to get the quasar objects
    c = SkyCoord(ra=ra_obj*u.degree, dec=dec_obj*u.degree)
    catalog = SkyCoord(ra=ra_cat*u.degree, dec=dec_cat*u.degree)
    idx, d2d, d3d = catalog.match_to_catalog_sky(c)
    quas = objs[idx]
    
    # SS Calculating the flux and mag 
    flux = objs['PSFFLUX']
    mag = np.array(22.5 - 2.5 * np.log10(flux))
    
    # SS Extinction Correction
    extnc = np.array(objs['EXTINCTION'])
    correct_mag = mag - extnc 

    # SS Calculating color cuts
    u_g = [row[0] - row[1] for row in correct_mag]
    g_r = [row[1] - row[2] for row in correct_mag]
    r_i = [row[2] - row[3] for row in correct_mag]
    i_z = [row[3] - row[4] for row in correct_mag]
    
    # SS Calculating the flux and mag for quasars
    quas_flux = quas['PSFFLUX']
    quas_mag = np.array(22.5 - 2.5 * np.log10(quas_flux))

    # SS Extinction Correction for quasars
    quas_extnc = np.array(quas['EXTINCTION'])
    quas_correct_mag = quas_mag - quas_extnc

    # SS Calculating color cuts for quasars
    quas_u_g = [row[0] - row[1] for row in quas_correct_mag]
    quas_g_r = [row[1] - row[2] for row in quas_correct_mag]
    quas_r_i = [row[2] - row[3] for row in quas_correct_mag]
    quas_i_z = [row[3] - row[4] for row in quas_correct_mag]   
        
    # SS PLotting r-i vs i-z color plot
    plt.figure(figsize=(10,6))
    plt.scatter(r_i, i_z, c='black')
    plt.scatter(quas_r_i, quas_i_z, c='red')
    #plt.axhline(y=-0.5, color = 'blue')
    #plt.axhline(y=0.5, color = 'blue')
    #plt.axvline(x=-1, color='blue')
    #plt.axvline(x=0.6, color='blue')
    plt.xlabel('r-i')
    plt.ylabel('i-z')
    plt.show()
    plt.savefig('color_plots/r-i_vs_i-z.png')
    
    # SS PLotting g-r vs u-g color plot
    plt.figure(figsize=(10,6))
    plt.scatter(g_r, u_g, c='black')
    plt.scatter(quas_g_r, quas_u_g, c='red')
    #plt.axhline(y=-0.2, color = 'blue')
    #plt.axhline(y=0.6, color = 'blue')
    #plt.axvline(x=-0.1, color='blue')
    #plt.axvline(x=0.4, color='blue')
    plt.xlabel('g-r')
    plt.ylabel('u-g')
    plt.show()
    plt.savefig('color_plots/u-g_vs_g-r.png')
    
    # SS PLotting w1/w1 (Wise mag) vs u-g color plot
    plt.scatter(objs["W1_NANOMAGGIES"] / objs["W2_NANOMAGGIES"], u_g, color = 'black')
    plt.scatter(quas["W1_NANOMAGGIES"] / quas["W2_NANOMAGGIES"], quas_u_g, color = 'red')
    #plt.axhline(y=-0.2, color = 'blue')
    #plt.axhline(y=0.6, color = 'blue')
    #plt.axvline(x=0.1, color='blue')
    #plt.axvline(x=0.6, color='blue')
    plt.xlabel('W1_W2')
    plt.ylabel('u_g')
    plt.xlim(-10, 10)
    plt.show()
    plt.savefig('color_plots/w1-w2_vs_u-g.png')

#Main
############
##########
if __name__ == '__main__':
    plot_get_colorcuts()





