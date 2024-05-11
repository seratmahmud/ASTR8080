# S. Saad
# v1 4/9/2024
# ASTR 8080 HW4


#Import
############
##########
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import warnings
from time import time
import matplotlib.pyplot as plt
from sdss_sweep_data_index import sdss_sweep_data_index
warnings.filterwarnings("ignore")


#Functions
############
##########


# SS Function to convert flux to magnitude
def flux_to_mag(flux):
    mag = 22.5 - 2.5 * np.log10(flux)
    return mag


# SS Function to apply color and flag cuts for quasar classification
def awesome_function(objs):
    """
    NAME: awesome_function

    PURPOSE: Returning the sources that are qsos

    CALLING SEQUENCE: from the Python prompt
   
      from awesome.py import awesome_function
      qso_array = awesome_function(objs)

    INPUTS:

      objs - a rec array that has columns RA and DEC, and the SDSS sweeps
      columns PSFFLUX, EXTINCTION, RESOLVE_STATUS, OBJC_TYPE, OBJC_FLAGS, 
      OBJC_FLAGS2, FLAGS, and FLAGS2 columns and the WISE W1_NANOMAGGIES 
      and W2_NANOMAGGIES columns.

    OUTPUTS:

      array - an array that has a total length of rows of the objs input
      and the positions which are possile quasars have value 1 in the array
      while other positions are indicated as 0 

    COMMENTS: You'll need to create the objs file using the create_test_data.py
    file in this repo

    REVISION HISTORY:

    v1.0: version Serat Saad, April 9, 2024 to submit as a homework for
    ASTR 8080 class at Vanderbilt University
    """
    
    # SS Recording time at the start
    start_time = time()
    
    # SS Correct for galactic extinction
    corrected_mag = flux_to_mag(objs['PSFFLUX']) - objs['EXTINCTION']
    
    # SS Calculate color indices
    g = corrected_mag[:,1]
    g_cut = g < 20
    
    u_g = corrected_mag[:,0] - corrected_mag[:,1]
    g_r = corrected_mag[:,1] - corrected_mag[:,2]
    r_i = corrected_mag[:,2] - corrected_mag[:,3]
    i_z = corrected_mag[:,3] - corrected_mag[:,4]
    
    # SS Call for wise magnitudes
    w1 = objs['W1_NANOMAGGIES']
    w2 = objs['W2_NANOMAGGIES']
    
    # SS Color Cut according to the find_color_cuts.py file and plots
    color_cut = (u_g > -0.2) & (u_g < 0.6) & (g_r > -0.1) & (g_r < 0.4) & \
      (i_z > -0.5) & (i_z < 0.5) & (r_i > -1) & (r_i < 0.6) &\
      (w1/w2 > 0.1) & (w1/w2 < 0.6) & g_cut
        
        
    
    # SS All necessary Flag cuts. These include whther the source is a point source,
    # whether it's a primary source and saturated + blend flags 
    point_source = (objs["OBJC_TYPE"] == 6)
    satflag = 2**18
    sat_data = (objs["OBJC_FLAGS"] & satflag) == 0
    blendflag = 2 ** 3
    blend_data = (objs["OBJC_FLAGS"] & blendflag) == 0
    
    # SS Final selection includes color cuts and flag filters
    final_selection = color_cut & sat_data & blend_data & point_source #& primary_data
    
    # SS Recording and printing time at the end
    end_time = time()
    print(f"Execution Time: {end_time - start_time:.2f} seconds")
    
    # SS Returning an array indicating quasar (1) or not (0)
    return np.where(final_selection, 1, 0)

