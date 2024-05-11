# Serat
# Homework 3 :)
# Data: 3/3/2024


# IMPORT BLOCK
###############################
###############################
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import pymangle
import matplotlib.pyplot as plt
import time


# FUNCTIONS
###############################
###############################


def sph_to_cart(ra_deg, dec_deg):
    ####
    # SS This function takes the sphecrical coordinates and converts it to  
    # cartesian coordinates, and returns x, y, and z :)
    ####
    c = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree, frame='icrs')
    x = np.cos(c.ra.radian) * np.cos(c.dec.radian)
    y = np.sin(c.ra.radian) * np.cos(c.dec.radian)
    z = np.sin(c.dec.radian)
    return x, y, z

def cap_with_radius(ra, dec, radius):
    ####
    # SS This function takes the ra, dec, and radius and returns
    # cartesian x, y, and z along with 1 - h value for that region :)
    ####
    x, y, z = sph_to_cart(ra, dec)
    one_minus_h = 1 - np.cos(np.radians(radius))
    return [x, y, z, one_minus_h]

def write_to_mangle(filename, caps, single_polygon=True):
    ####
    # SS This function takes the filename and cap values (x, y, z, 1-h)
    # and writes a .ply file in the mangle format :)
    ####
    with open(filename, 'w') as f:
        if single_polygon:
            f.write('1 polygons\n')
            for cap in caps:
                f.write(f'polygon 1 ( {len(caps)} caps, 1 weight, 0 pixel, 0 str):\n')
                f.write(' '.join(map(str, cap)) + '\n')
        else:
            f.write(f'{len(caps)} polygons\n')
            for i, cap in enumerate(caps, start=1):
                f.write(f'polygon {i} ( 1 caps, 1 weight, 0 pixel, 0 str):\n')
                f.write(' '.join(map(str, cap)) + '\n')

def random_points_on_sphere(n_points):
    ####
    # SS This function takes a number and generate that number of 
    # amounts of ra and dec values uniformly spread over the sphere :)
    ####
    ra = np.random.uniform(0, 360, n_points)
    dec = np.arcsin(np.random.uniform(-1, 1, n_points)) * (180 / np.pi)
    return ra, dec

def get_quasar_file(file_path):
    ####
    # SS This function takes the file path of quasar file and returns 
    # the ra and dec files reading it from that file :)
    ####
    with open(file_path, 'r') as file:
        ra = []
        dec = []
        for line in file:
            ra_hms = line[:2] + "h" + line[2:4] + "m" + line[4:8] + "s"
            dec_dms = line[9:12] + "d" + line[12:14] + "m" +\
                line[15:].strip() + "s"
            ra.append(ra_hms)
            dec.append(dec_dms)
        quasar_coord = SkyCoord(ra, dec, frame='icrs')
    return quasar_coord

def plot_quasars(quasars, in_mask, mask_area):
    ####
    # SS This function takes the qusar coordinates, mask, and mask area and returns 
    # a plot of qusars inside and outside of the mask along with showing number density
    # on the plot :)
    ####
    plt.figure(figsize=(10, 5))
    outside_mask = [quasar for quasar, mask in zip(quasars, in_mask) if not mask]
    inside_mask = [quasar for quasar, mask in zip(quasars, in_mask) if mask]
    ra_outside, dec_outside = zip(*[(q.ra.degree, q.dec.degree) for q in outside_mask])
    ra_inside, dec_inside = zip(*[(q.ra.degree, q.dec.degree) for q in inside_mask])
    # SS Plotting
    plt.scatter(ra_outside, dec_outside, color='red', s=1, label='Outside Mask')
    plt.scatter(ra_inside, dec_inside, color='blue', s=1, label='Inside Mask')
    plt.xlabel('Right Ascension (degrees)')
    plt.ylabel('Declination (degrees)')
    plt.title(f'Quasar Distribution (Mask Area: {mask_area:.2f} sq. deg)')
    plt.legend()
    number_density = len(inside_mask) / mask_area
    plt.figtext(0.5, 0.01, f'Number Density: {number_density:.2f} quasars/sq. deg', ha='center')
    plt.savefig('quasar_distribution.png')

def main():
    ##### Problem 1
    #####
    epoch_0 = time.time()
    # SS Writing an array for the required RA, DEC, and radius for the plates
    plate_centers = [(155, 34), (159, 36), (163, 34), (167, 36)]
    plate_radius = 2
    # SS Saving the x, y, z, and 1-h values for all the four possibilities
    caps_manual = []
    for ra, dec in plate_centers:
        caps_manual.append(cap_with_radius(ra, dec, plate_radius))
    # SS writing a new file in .ply format    
    write_to_mangle('mask_manual.ply', caps_manual, single_polygon=False)

    epoch_1 = time.time()
    print("Time taken by Problem 1 in seconds:", epoch_1-epoch_0)

    ##### Problem 2
    #####
    # SS generating the coordinates for 100,000,000 random points uniformly
    # distributed on the sphere
    n_points = 100_000_000
    ra_random, dec_random = random_points_on_sphere(n_points)
    mask = pymangle.Mangle("mask_manual.ply")
    in_mask = mask.contains(ra_random, dec_random)
    # SS Getting the area of the cap
    area_ratio = np.sum(in_mask) / n_points
    total_sphere_area = 4 * np.pi * (180/np.pi)**2 # SS area in square degrees
    estimated_mask_area = total_sphere_area * area_ratio
    print(f'Estimated mask area: {estimated_mask_area} square degrees')

    epoch_2 = time.time()
    print("Time taken by Problem 2 in seconds:", epoch_2-epoch_1)

    ##### Problem 3
    #####
    # SS Getting the coordinates of Quasar
    quasar_coords = get_quasar_file('HW3quasarfile.dat')
    in_mask = mask.contains(quasar_coords.ra.degree, quasar_coords.dec.degree)
    mask_area = estimated_mask_area 
    plot_quasars(quasar_coords, in_mask, mask_area) # SS Plotting all the quasars

    epoch_3 = time.time()
    print("Time taken by Problem 3 in seconds:", epoch_3-epoch_2)
    print("Total time taken in seconds:", epoch_3-epoch_0)

# MAIN
###############################
###############################
if __name__ == '__main__':
    main()