# S. Saad
# v1 2/18/2024
# ASTR 8080 HW2

# IMPORT BLOCK
###############################
###############################
from astropy.coordinates import SkyCoord 
from astropy.io import fits
import healpy as hp

import matplotlib.pyplot as plt
import numpy as np
import warnings
from astropy.utils.exceptions import ErfaWarning
warnings.filterwarnings('ignore', category=ErfaWarning)

# FUNCTIONS
###############################
###############################


# SS This function calculate area of the field
def area_of_field(alpha_min, alpha_max, delta_min, delta_max):
    # SS Convert declination to radians for calculation
    delta_min_rad = np.radians(delta_min)
    delta_max_rad = np.radians(delta_max)
    # SS Calculate the area in steradians and then convert to square degrees
    area_rad = 2 * np.pi * (np.sin(delta_max_rad) - np.sin(delta_min_rad))\
        * (alpha_max - alpha_min) / 360
    area_deg = np.degrees(area_rad) * (360 / (2 * np.pi))
    return area_deg


# SS Populates a square field with random points within given celestial coordinates
def populate_square_field(alpha_min, alpha_max, delta_min, delta_max, num_points):
    total_area_sphere = 4 * np.pi * (180 / np.pi) ** 2  # SS Calculate total sphere area in sq. degrees
    area_square_field = area_of_field(alpha_min, alpha_max, delta_min, delta_max)  # SS Get area of specified field
    proportion = area_square_field / total_area_sphere  # SS Calculate proportion of field area to total sphere
    num_points_square_field = int(num_points * proportion)  # SS Scale number of points by area proportion
    # SS Generate random points within the field
    random_alphas = np.random.uniform(alpha_min, alpha_max, num_points_square_field)
    random_sin_deltas = np.random.uniform(np.sin(np.radians(delta_min)), np.sin(np.radians(delta_max)), num_points_square_field)
    random_deltas = np.degrees(np.arcsin(random_sin_deltas))
    return list(zip(random_alphas, random_deltas))


# SS Reads quasar coordinates from a file and returns them as SkyCoord objects
def get_quasar_file(file_path):
    quasars = []
    with open(file_path, 'r') as file:
        for line in file:
            ra_hms = line[:2] + "h" + line[2:4] + "m" + line[4:8] + "s"
            dec_dms = line[9:12] + "d" + line[12:14] + "m" + line[15:].strip() + "s"
            quasar_coord = SkyCoord(ra_hms, dec_dms, frame='icrs')
            quasars.append(quasar_coord)
    return quasars


# SS Calculate HEALPix pixels for a list of quasar coordinates at different resolutions
def calculate_healpix_pixels(quasars, nsides=[4, 8, 16]):
    pixnums = []
    for quasar in quasars:
        pixnum_for_nsides = [hp.ang2pix(nside, quasar.ra.degree, quasar.dec.degree, lonlat=True) for nside in nsides]
        pixnums.append(pixnum_for_nsides)
    return pixnums


# SS Plotting function for quasar locations and identifying over-dense pixels
def plot_file(file_name):
    with fits.open(file_name) as hdul:
        data = hdul[1].data
        ras = data['ra']
        decs = data['dec']
        pixnums = data['pixnum'][:, 0]

    plt.figure(figsize=(10, 5))
    plt.scatter(ras, decs, s=1, color='blue', label='All Quasars')
    plt.xlabel('Right Ascension')
    plt.ylabel('Declination')
    plt.title('Quasar Locations and Over-Dense Pixels')

    # SS Identify and highlight over-dense pixels
    unique_pixels, counts = np.unique(pixnums, return_counts=True)
    over_dense_indices = np.argsort(counts)[-5:]
    over_dense_pixels = unique_pixels[over_dense_indices]

    color_list = ['black', 'green', 'red', 'yellow', 'purple']
    for i in range(len(over_dense_pixels)):
        over_dense_mask = pixnums == over_dense_pixels[i]
        plt.scatter(ras[over_dense_mask], decs[over_dense_mask],\
                    s=10, color=color_list[i], label=f'Over-Dense Pixel {over_dense_pixels[i]}')

    plt.legend()
    plt.show()


def main():
    # SS Setting initial parameters for area calculation
    alpha_min = 0
    alpha_max = 30

    # SS Defining different declination regions for area calculation
    regions = [(0, 30),
        (30, 60),
        (60, 75),
        (75, 90)]

    # SS Calculating areas for different regions
    areas = []
    for delta_min, delta_max in regions:
        area = area_of_field(alpha_min, alpha_max, delta_min, delta_max)
        areas.append(area)

    # SS Plotting areas for different declination regions
    plt.figure(figsize=(10, 6))
    for i, ((delta_min, delta_max), area) in enumerate(zip(regions, areas)):
        plt.plot([alpha_min, alpha_max], [i, i], marker='o', label=f"δ: {delta_min}° to {delta_max}°, Area: {area:.2f} sq. deg")

    plt.yticks(range(len(areas)), [f"Region {i+1}" for i in range(len(areas))])
    plt.xlabel("Right Ascension (degrees)")
    plt.ylabel("Regions")
    plt.title("Areas of Fields at Different Declinations")
    plt.legend()
    plt.grid(True)
    plt.show()

    # SS Calculating area of a spherical cap
    spherical_cap_area = area_of_field(0, 360, 0, 90)
    print(spherical_cap_area)
    
    # SS Populating a square field with random points
    alpha_min, alpha_max, delta_min, delta_max = 0, 360, 0, 90  # A spherical cap
    num_points = 100000  # Total number of points for the entire sphere
    points = populate_square_field(alpha_min, alpha_max, delta_min, delta_max, num_points)
    print(len(points), len(points) / num_points)
    
    # SS Reading quasar data from file and calculating HEALPix pixels
    quasars = get_quasar_file('HW1quasarfile.dat')
    pixnums = calculate_healpix_pixels(quasars)

    # SS Creating a FITS file with quasar data and their HEALPix pixel numbers
    dtype = [('ra', 'float'), ('dec', 'float'), ('pixnum', '3int')]
    records = np.recarray((len(quasars),), dtype=dtype)
    for i, quasar in enumerate(quasars):
        records[i] = (quasar.ra.degree, quasar.dec.degree, pixnums[i])
    fits.writeto('quasars_with_pixnums.fits', records, overwrite=True)
    
    # SS Plotting quasar locations and identifying over-dense pixels
    plot_file('quasars_with_pixnums.fits')
    
    
# MAIN
###############################
###############################
if __name__ == '__main__':
    main()