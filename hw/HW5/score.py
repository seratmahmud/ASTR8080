# S. Saad
# v1 4/9/2024
# ASTR 8080 HW4

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import warnings
from time import time
import matplotlib.pyplot as plt
from sdss_sweep_data_index import sdss_sweep_data_index
from awesome_homework import awesome_function
from create_test_data import fetch_sweep_objects
warnings.filterwarnings("ignore")

#Functions
############
##########
def calculate_score(test_ra, test_dec, test_radius, test_file):
    """
    NAME: calculate_score

    PURPOSE: Returning the score that I should get in the awesome_homework

    INPUTS:

      test_ra - Central Right Ascension of the test quasar file area
      test_dec - Central Declination of the test quasar file area
      test_radius - Radius of the test quasar file area from the center
      test_file - Locaiton of the test_file

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

    # SS Read the new QSOS file to get the coordinates
    # SS Reading files vary depending on the format of the given file
    hdul_test = fits.open(test_file)
    test_qsos_data = hdul_test[1].data
    #test_qsos_data = test_qsos_data[0] # SS This line only apply for the qsos-ra-130-dec30-rad3.fits file
    hdul_test.close()
    test_qsos_coords = SkyCoord(ra=test_qsos_data['RA']*u.degree, dec=test_qsos_data['DEC']*u.degree)

    # SS Fetch corresponding objects from the SDSS sweep files
    # SS Or call from local fits file if used before
    #test_objs = fetch_sweep_objects(test_ra, test_dec, test_radius)
    hdul_objs = fits.open('/astro/astr8020/wise-stargal-sweeps-ra138-dec46-rad3.fits')
    test_objs = hdul_objs[1].data
    hdul_objs.close()

    # SS Get time before calling classification funciton
    start_time = time()    

    # SS Apply the classification function to get the array
    classification_result = awesome_function(test_objs)

    # SS Get time after calling classification funciton
    end_time = time()    

    # SS Get the coordinates of the quasars
    classified_qsos = test_objs[np.where(classification_result == 1)]
    classified_qsos_coords = SkyCoord(ra=classified_qsos['RA']*u.degree, dec=classified_qsos['DEC']*u.degree)
    
    # SS Matching to return the matches for classified qsos and test qsos
    match, _, _, _ = test_qsos_coords.search_around_sky(classified_qsos_coords, 1/60/60*u.degree)
    

    # SS Calculating q, t, and f
    area_sq_deg = np.pi * (test_radius ** 2)
    q = len(match) / area_sq_deg 
    t = len(classified_qsos_coords) / area_sq_deg
    print(q, t)
    
    
    f = min(q / t if t != 0 else 0, 0.85)
    q = min(q, 10)

    # SS Calculating accuracy score
    accuracy_score = 30 * (f / 0.85) * (q / 10)

    # SS Printing Number of Quasar, targets and accuracy score
    print(f"q: {q}")
    print(f"t: {t}")
    print(f"Accuracy Score: {accuracy_score:.2f} out of 30.00 points")

    # SS Calculate and printthe speed score
    total_time = end_time - start_time
    t = max(2, min(total_time, 12))
    speed_score = 12 - t
    print(f"Speed Score: {speed_score:.2f} out of 10.00 points")

    # SS Returning accuracy and speed score
    return accuracy_score, speed_score


#Main
############
##########
if __name__ == '__main__':
    # SS Writing the coordinate and radius in a varibale
    test_ra = 138
    test_dec = 46
    test_radius = 3

    # SS Calling the calcualte_score function to return scores
    accuracy_score, speed_score = calculate_score(
    test_ra=test_ra,
    test_dec=test_dec,
    test_radius=test_radius,
    test_file="/astro/astr8020/qsos-ra138-dec46-rad3.fits") 
    #test_file=f"/home/saadsm/ASTR8080/runnoe/week13/moreqsos-ra{test_ra}-dec{test_dec}-rad{test_radius}.fits") 

