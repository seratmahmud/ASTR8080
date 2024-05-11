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
from awesome_homework import awesome_function
warnings.filterwarnings("ignore")


#Functions
############
##########


def fetch_sweep_objects(ra, dec, radius, objtype='star', sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final'):
    """
    NAME: fetch_sweep_objects

    PURPOSE: Get the data from the SDSS Sweep files as object

    CALLING SEQUENCE: from the Python prompt
   
      from creat_test_data.py import fetch_sweep_objects
      test_objs = fetch_sweep_objects(ra, dec, radius, objtype, sweepdir)

    INPUTS:

      ra - Right Ascension
      dec - Declination
      radius - Radius of the targeted region
      objtype - Object Type (Should be stra in this case)
      sweepdir - Directory where the sweepfiles are at

    OUTPUTS:

      objs - a rec array that has columns RA and DEC, and the SDSS sweeps
      columns PSFFLUX, EXTINCTION, RESOLVE_STATUS, OBJC_TYPE, OBJC_FLAGS, 
      OBJC_FLAGS2, FLAGS, and FLAGS2 columns and the WISE W1_NANOMAGGIES 
      and W2_NANOMAGGIES columns.

    REVISION HISTORY:

    v1.0: version Serat Saad, April 9, 2024 to submit as a homework for
    ASTR 8080 class at Vanderbilt University
    """
    
    # SS Calling the sweep files, then doing string operation to get the sdss sweep files and wise sweep files 
    swfiles = sdss_sweep_data_index(ra, dec, radius, objtype=objtype, sweepdir=sweepdir, all=True, verbose=False)
    star_swfiles = [file.replace("star", "stargal-primary") for file in swfiles]
    wise_swfiles = [file.replace("star", "wise-stargal-primary") for file in swfiles]
    
    # SS Getting the SDSS star objects
    star_objs = []
    for file in star_swfiles:
        try:
            with fits.open(file) as hdul:
                star_objs.append(hdul[1].data)
        except:
            continue
    star_objs = np.hstack(star_objs)

    # SS Getting the WISE Objects
    wise_objs = []
    for file in wise_swfiles:
        try:
            with fits.open(file) as hdul:  
                wise_objs.append(hdul[1].data)
        except:
            continue
    wise_objs = np.hstack(wise_objs)
    
    # SS Matching the SDSS Star Objects and WISE Objects to get only the objects that have data in both of them
    star_target = SkyCoord(star_objs["RA"], star_objs["DEC"], unit=(u.degree, u.degree))
    wise_target = SkyCoord(wise_objs["RA"], wise_objs["DEC"], unit=(u.degree, u.degree))
    match = np.where(star_target.separation(wise_target) < radius*u.degree) [0]
    wise_objs, star_objs = wise_objs[match], star_objs[match]

    # SS Declaring the data name and type to save the data from SDSS Sweep and WISE in the same Object
    type = np.dtype([("RA", 'f8'), ("DEC", 'f8'), ('PSFFLUX', 'f8', (5,)),\
        ("EXTINCTION", 'f8', (5,)), ("RESOLVE_STATUS", 'i8'), ("OBJC_TYPE", 'i8'), \
        ("OBJC_FLAGS", 'i8'), ("OBJC_FLAGS2", 'i8'), ("FLAGS", 'i8', (5,)),\
        ("FLAGS2", 'i8', (5,)), ("W1_NANOMAGGIES", 'f8'), ("W2_NANOMAGGIES", 'f8')])
    
    # SS Inserting data from both SDSS sweep file and WISE file in the object
    test_objs = np.zeros(len(star_objs), dtype=type)
    inputs = ["RA", "DEC", "PSFFLUX", "EXTINCTION", "RESOLVE_STATUS", "OBJC_TYPE", "OBJC_FLAGS", "OBJC_FLAGS2", "FLAGS", "FLAGS2"]
    for str in inputs:
        test_objs[str] = star_objs[str]
    test_objs["W1_NANOMAGGIES"] = wise_objs["W1_NANOMAGGIES"]
    test_objs["W2_NANOMAGGIES"] = wise_objs["W2_NANOMAGGIES"]
    
    # SS Only getitng the sources that are withing the area
    cin = SkyCoord(ra=ra*u.degree, dec=dec*u.degree) 
    csweeps = SkyCoord(ra=test_objs["RA"]*u.degree, dec=test_objs["DEC"]*u.degree)
    sep = cin.separation(csweeps)
    idx = np.where(sep < radius*u.degree)[0]
    test_objs = test_objs[idx]
    
    # SS Saving the object in a fits file to reuse later
    fits_file = fits.BinTableHDU(test_objs)
    fits_file.writeto(f'fits_reuse_file/sources-ra{ra}-dec{dec}-rad{radius}.fits', overwrite=True)
 
    # SS returning test_objs
    return test_objs
