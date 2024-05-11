# Final Exam
# ASTR8080
# Serat Saad
# v1 5/1/2024


#Import
############
##########
import numpy as np
from astropy.coordinates import SkyCoord, Angle
from numpy.random import random
from astropy import units as u
from astropy.io import fits
import warnings
import time
import matplotlib.pyplot as plt
import sys
import pymangle

sys.path.append('/home/saadsm/ASTR8080/serat/hw/HW3/')
from hw3 import write_to_mangle, cap_with_radius

warnings.filterwarnings("ignore")


#Functions
############
##########


def sweep_data_rectangular(ra_min, ra_max, dec_min, dec_max, mag_lim, extinction,\
    objtype, sweepdir='.', all=False, verbose=False):
    """
    NAME:
        sweep_data_rectangular

    PURPOSE:
        To fetch and process astronomical data within a specified rectangular area from SDSS sweep files.

    CALLING SEQUENCE:
        objs = sweep_data_rectangular(ra_min, ra_max, dec_min, dec_max, objtype, sweepdir='.', all=False, verbose=False)

    INPUTS:
        ra_min - The minimum right ascension (degrees)
        ra_max - The maximum right ascension (degrees)
        dec_min - The minimum declination (degrees)
        dec_max - The maximum declination (degrees)
        objtype - Type of objects to fetch ('star' or 'gal')
        sweepdir - Directory containing sweep files
        all - Boolean flag to consider all data points, default is False
        verbose - Boolean flag for verbose output, default is False

    OUTPUTS:
        A structured NumPy array containing processed object data.

    COMMENTS:
        Utilizes astropy to handle FITS files and data structures.
        Requires sdss_sweep_data_index for indexing sweeps.

    REVISION HISTORY:
        v1.0: Serat Saad, May 1, 2024 - Initial version for ASTR 8080 Final
    """
    
    indexfile = f'{sweepdir}/datasweep-index-{objtype}.fits'
    
    with fits.open(indexfile) as index_struct:
        index = index_struct[1].data

        if not all:
            index = index[index["NPRIMARY"] > 0]

        # SS Defining the limits of the rectangular section
        idx = (index["RA"] >= ra_min) & (index["RA"] <= ra_max) & \
              (index["DEC"] >= dec_min) & (index["DEC"] <= dec_max)

        if not idx.any():
            raise IOError('no objects in specified rectangular area of SDSS')

        swinfo = index[idx]

    # SS Generating the list of FITS file names to open
    swfiles = [f'{sweepdir}/{sw["RERUN"]}/calibObj-{sw["RUN"]:06}-{sw["CAMCOL"]}-{objtype}.fits.gz'
               for sw in swinfo if sw['RERUN'] == '301']
    
    swfiles = np.array(swfiles)
    swfiles = np.unique(swfiles)
    # SS Opening each FITS file, read the data, and concatenate
    objs = []
    for file in swfiles:
      try:
        with fits.open(file) as obj_struct:
            objs.append(obj_struct[1].data)
      except:
        continue

    # SS Concatenating all fetched data objects into a single array
    if objs:
        objs = np.hstack(objs)
    else:
        objs = np.array([])

    if verbose:
        print('Done...t = {}s'.format(time()-t0))

    # SS Declaring the data name and type to save the data from SDSS Sweep and WISE in the same Object
    type = np.dtype([("RA", 'f8'), ("DEC", 'f8'), ('PSFFLUX', 'f8', (5,)),\
        ("EXTINCTION", 'f8', (5,))])
    
    # SS Inserting data from both SDSS sweep file and WISE file in the object
    test_objs = np.zeros(len(objs), dtype=type)
    inputs = ["RA", "DEC", "PSFFLUX", "EXTINCTION"]
    for str in inputs:
        test_objs[str] = objs[str]
    
    if extinction==True:
        psfmag = 22.5 - 2.5 * np.log10(test_objs['PSFFLUX']) - test_objs['EXTINCTION']
    else:
        psfmag = 22.5 - 2.5 * np.log10(test_objs['PSFFLUX'])
        
    filter_idx = psfmag[:, 3] < mag_lim
    filtered_objs = test_objs[filter_idx]
    
    # SS Saving the object in a fits file to reuse later
    fits_file = fits.BinTableHDU(filtered_objs)
    fits_file.writeto(f'sources-{objtype}.fits', overwrite=True)
    return filtered_objs


def create_caps_for_mangle(star_data):
    """
    NAME:
        create_caps_for_mangle

    PURPOSE:
        Generate caps from star data to use with a mangle mask.

    CALLING SEQUENCE:
        caps = create_caps_for_mangle(star_data)

    INPUTS:
        star_data - Array of star data including positions and photometry.

    OUTPUTS:
        caps - List of caps, each representing a region on the sky.

    COMMENTS:
        The function processes star brightness to filter and create spherical caps.
        
    REVISION HISTORY:
        v1.0: Serat Saad, May 1, 2024 - Initial version for ASTR 8080 Final
    """
    
    # SS Radius in degrees
    radius_deg = 5 / 3600
    # SS Generating caps using filtered star data
    caps = [cap_with_radius(star['RA'], star['DEC'], radius_deg) for star in star_data]
    return caps


def load_and_apply_mangle_mask(filename, galaxy_data):
    """
    NAME:
        load_and_apply_mangle_mask

    PURPOSE:
        Apply a mangle mask to galaxy data to filter galaxies within specified regions.

    CALLING SEQUENCE:
        galaxies_in_mask, filtered_galaxies = load_and_apply_mangle_mask(filename, galaxy_data)

    INPUTS:
        filename - Path to the mangle mask file.
        galaxy_data - Array of galaxy data to filter.

    OUTPUTS:
        galaxies_in_mask - Galaxies that fall within the mask.
        filtered_galaxies - Galaxies that were considered for masking.

    COMMENTS:
        Filters galaxies based on a provided mangle mask.

    REVISION HISTORY:
        v1.0: Serat Saad, May 1, 2024 - Initial version for ASTR 8080 Final
    """
    
    # SS Creating a mask
    mangle_mask = pymangle.Mangle(filename)
    
    #SS Finding the sources that are in the mask
    mask_containment = mangle_mask.contains(galaxy_data['RA'], galaxy_data['DEC'])
    galaxies_in_mask = galaxy_data[mask_containment]
    return galaxies_in_mask


def calculate_mask_area(filename, region_area, ra_max, ra_min, dec_max, dec_min):
    """
    NAME:
        calculate_mask_area

    PURPOSE:
        Calculate the total area covered by spherical caps.

    CALLING SEQUENCE:
        total_area = calculate_mask_area(caps)

    INPUTS:
        filename - filename of the .ply file.
        region_area - the area of the rectangular region
        ra_min - Minimum right ascension.
        ra_max - Maximum right ascension.
        dec_min - Minimum declination.
        dec_max - Maximum declination.
        
    OUTPUTS:
        total_area - Total area covered by the caps in square degrees.

    COMMENTS:
        Summarizes the area of caps to provide total coverage area.

    REVISION HISTORY:
        v1.0: Serat Saad, May 2, 2024 - Initial version for ASTR 8080 Final
    """
    # SS Calculating the mask area
    n_points = 1000000
    ra_random = (random(n_points)*(ra_max - ra_min)) + (ra_min)
    dec_random = (180./np.pi) * (np.arcsin(random(n_points)*(np.sin(dec_max *(np.pi/180))-(np.sin(dec_min *(np.pi/180))))+np.sin(dec_min*(np.pi/180))))
    
    mask = pymangle.Mangle(filename)
    in_mask = mask.contains(ra_random, dec_random)
    area_ratio = np.sum(in_mask) / n_points
    mask_area = region_area * area_ratio
    return mask_area


def calculate_sky_area(ra_min, ra_max, dec_min, dec_max):
    """
    NAME:
        calculate_sky_area

    PURPOSE:
        Calculate the sky area in square degrees within specified RA and Dec bounds.

    CALLING SEQUENCE:
        area = calculate_sky_area(ra_min, ra_max, dec_min, dec_max)

    INPUTS:
        ra_min - Minimum right ascension.
        ra_max - Maximum right ascension.
        dec_min - Minimum declination.
        dec_max - Maximum declination.

    OUTPUTS:
        area - Total sky area in square degrees.

    COMMENTS:
        Computes the area using spherical trigonometry principles.

    REVISION HISTORY:
        v1.0: Serat Saad, May 1, 2024 - Initial version for ASTR 8080 Final
    """
    
    dec_min_rad = np.radians(dec_min)
    dec_max_rad = np.radians(dec_max)
    area = 2 * np.pi * (np.sin(dec_max_rad) - np.sin(dec_min_rad)) * (ra_max - ra_min) / 360
    return area

def estimate_mask_area(ra_min, ra_max, dec_min, dec_max, filename, region_area, n_simulations=10):
    """
    NAME:
        estimate_mask_area

    PURPOSE:
        Calculate the sky area with uncertainty in square radians within specified RA and Dec bounds.

    CALLING SEQUENCE:
        area = estimate_mask_area(ra_min, ra_max, dec_min, dec_max, filename, region_area, n_simulations=10)

    INPUTS:
        ra_min - Minimum right ascension.
        ra_max - Maximum right ascension.
        dec_min - Minimum declination.
        dec_max - Maximum declination.
        filename - filename of the .ply file.
        region_area - the area of the rectangular region
        n_simulations - total number of simulations being used to calculate the uncertainty

    OUTPUTS:
        mean_area - Total sky area in square degrees.
        st_dev_area - The uncertainty in the area calculation.

    COMMENTS:
        Computes the area using spherical trigonometry principles.

    REVISION HISTORY:
        v1.0: Serat Saad, May 1, 2024 - Initial version for ASTR 8080 Final
    """
    n_points = 100000 # SS You can change this and try a large number for more accurate predictions
    area_estimates = []

    for i in range(n_simulations):
        # SS Generating random points within the given RA and Dec ranges
        ra_random = np.random.uniform(ra_min, ra_max, n_points)
        dec_random = (180. / np.pi) * np.arcsin(
            np.random.uniform(np.sin(dec_min * (np.pi / 180)), np.sin(dec_max * (np.pi / 180)), n_points) +
            np.sin(dec_min * (np.pi / 180))
        )

        # SS Loading mask from Mangle and check which points are within the mask
        mask = pymangle.Mangle(filename)
        in_mask = mask.contains(ra_random, dec_random)

        # SS Calculate area based on ratio of points within the mask
        area_ratio = np.sum(in_mask) / n_points
        mask_area = region_area * area_ratio
        area_estimates.append(mask_area)

    # SS Calculating mean and standard deviation of the estimated areas
    mean_area = np.mean(area_estimates)
    std_dev_area = np.std(area_estimates)

    return mean_area, std_dev_area


def plot_aitoff(filtered, in_mask):
    """
    NAME:
        plot_aitoff

    PURPOSE:
        Plot data points in an Aitoff projection to visualize celestial data distribution.

    CALLING SEQUENCE:
        plot_aitoff(filtered, in_mask)

    INPUTS:
        filtered - Dataset of filtered galaxy positions.
        in_mask - Dataset of galaxy positions within a mask.

    OUTPUTS:
        None - Generates a plot.

    COMMENTS:
        Visualizes the distribution of galaxies, highlighting differences between filtered datasets.

    REVISION HISTORY:
        v1.0: Serat Saad, May 1, 2024 - Initial version for ASTR 8080 Final
    """
    # SS Calling the ra and decs
    ra1 = filtered['RA']
    ra2 = in_mask['RA']
    dec1 = filtered['DEC']
    dec2 = in_mask['DEC']
    
    # SS Modifying the ra's for the aitoff plot
    ra1 = ((ra1 + 180) % 360) - 180
    ra2 = ((ra2 + 180) % 360) - 180

    # SS Converting ra's and dec's to radians
    ra1 = np.radians(ra1)
    ra2 = np.radians(ra2)
    dec1 = np.radians(dec1)
    dec2 = np.radians(dec2)
    
    # SS Plotting in aitoff projection
    plt.figure(figsize=(10, 6))
    ax = plt.subplot(111, projection='aitoff')
    ax.scatter(ra1, dec1, color='Black')
    ax.scatter(ra2, dec2, color='Red')
    x_labels_hours = ['14h', '16h', '18h', '20h', '22h', '0h', '2h', '4h', '6h', '8h', '10h']
    ax.set_xticklabels(x_labels_hours, weight=800)
    plt.title('Tergetted Galaxies inside the mask (red)')
    plt.show()
    plt.savefig('aitoff.png')


def main():
    start = time.time()

    # SS Calling the ra and dec max and min values
    ra_max = 210
    ra_min = 150
    dec_max = 60
    dec_min = 30
    
    # SS Getting the area of the whole region
    region_area = calculate_sky_area(ra_min, ra_max, dec_min, dec_max)
    
    # SS Should use the following code if the fits files don't exist
    
    #star_objs = sweep_data_rectangular(ra_min, ra_max, dec_min, dec_max, mag_lim = 10, extinction=False, objtype='star',\
    #    sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final', all=False, verbose=False)
    #gal_objs = sweep_data_rectangular(ra_min, ra_max, dec_min, dec_max, mag_lim=19, extinction=True, objtype='gal',\
    #    sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final', all=False, verbose=False)
    
    # SS Callling the fits files to get the star and galaxy data
    hdul = fits.open('sources-gal.fits')
    gal_data = hdul[1].data
    hdul = fits.open('sources-star.fits')
    star_data = hdul[1].data
    
    # SS Getting the caps for the bright stars
    bright_star_caps = create_caps_for_mangle(star_data)
    
    # SS Writing the bright star caps to a mangle file
    write_to_mangle('bright_star_mask.ply', bright_star_caps, single_polygon=False)

    # SS Calculating the mask area & uncertainty (if you want)
    mask_area = calculate_mask_area('bright_star_mask.ply', region_area, ra_max, ra_min, dec_max, dec_min)
    #mask_area, mask_area_std = estimate_mask_area(ra_min, ra_max, dec_min, dec_max, 'bright_star_mask.ply', region_area, n_simulations=3)
    
    # SS Printing the total area and the mask area and uncertainty (if you want)
    print(f'Total area: {region_area} square radians')
    print(f'Total mask area: {mask_area} square radians')
    #print(f'Uncertainty in mask area calculation: {mask_area} square radians')
    
    # SS Getting the galaxies that are less than 19 mag and plotting the galaxies in the mask
    galaxies_in_mask = load_and_apply_mangle_mask('bright_star_mask.ply', gal_data)
    plot_aitoff(gal_data, galaxies_in_mask)

    # SS Calculating the total number density and the number density by removing galaxies in the mask
    number_density_total = len(gal_data) / region_area
    number_density_adjusted = (len(gal_data) - len(galaxies_in_mask)) / (region_area-mask_area)
    
    # SS Printing number desnities and total number of galaxies in the mask
    print(f'Total number of galaxies in the mask: {len(galaxies_in_mask)}')
    print(f'Number density total: {number_density_total} galaxies per square radians')
    print(f'Number density with mask: {number_density_adjusted} galaxies per square radians')
    
    # SS Printing the total time to run the program
    end = time.time()
    print(f'Total time taken to run the file: {end-start} seconds.')


# MAIN
###############################
###############################
if __name__ == '__main__':
    main()
    
    
    
    
