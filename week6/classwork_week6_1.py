# S. Saad
# ASTR 8080 week6, classwork2

# IMPORT BLOCK
###############################
###############################


from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import pymangle
import matplotlib.pyplot as plt





# FUNCTIONS
###############################
###############################


def sph_to_cart(ra_deg, dec_deg):
    c = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree, frame='icrs')
    x = np.cos(c.ra.radian) * np.cos(c.dec.radian)
    y = np.sin(c.ra.radian) * np.cos(c.dec.radian)
    z = np.sin(c.dec.radian)
    return x, y, z

def cap_with_radius(ra, dec, radius):
    x, y, z = sph_to_cart(ra, dec)
    # Calculate 1 - h from the radius
    one_minus_h = 1 - np.cos(np.radians(radius))
    return [x, y, z, one_minus_h]

def cap_with_flipped_radius(ra, dec, radius):
    x, y, z = sph_to_cart(ra, dec)
    # Calculate 1 - h from the radius
    one_minus_h = 1 - np.cos(np.radians(radius))
    return [x, y, z, -one_minus_h]

def write_to_mangle(filename, caps, single_polygon=True):
    with open(filename, 'w') as f:
        if single_polygon:
            f.write('1 polygons\n')
            f.write(f'polygon 1 ( {len(caps)} caps, 1 weight, 0 pixel, 0 str):\n')
            for cap in caps:
                f.write(' '.join(map(str, cap)) + '\n')
        else:
            f.write(f'{len(caps)} polygons\n')
            for i, cap in enumerate(caps, start=1):
                f.write(f'polygon {i} ( 1 caps, 1 weight, 0 pixel, 0 str):\n')
                f.write(' '.join(map(str, cap)) + '\n')

def generate_and_plot_random_points():
    minter = pymangle.Mangle('intersection.ply')
    mbothcaps = pymangle.Mangle('bothcaps.ply')

    ra_minter, dec_minter = minter.genrand(1000)
    ra_mbothcaps, dec_mbothcaps = mbothcaps.genrand(1000)

    plt.figure(figsize=(10, 6))
    plt.scatter(ra_minter, dec_minter, color='blue', label='Intersection')
    plt.scatter(ra_mbothcaps, dec_mbothcaps, color='red', label='Both Caps')
    
    plt.xlabel('Right Ascension (degrees)')
    plt.ylabel('Declination (degrees)')
    plt.title('Random Points from Intersection and Both Caps')
    plt.legend()
    plt.show()

def flipped():
    flipped_1 = pymangle.Mangle('flipped1.ply')
    flipped_2 = pymangle.Mangle('flipped2.ply')

    ra_flipped_1, dec_flipped_1 = flipped_1.genrand(1000)
    ra_flipped_2, dec_flipped_2 = flipped_2.genrand(1000)

    plt.figure(figsize=(10, 6))
    plt.scatter(ra_flipped_1, dec_flipped_1, color='blue')
    #plt.scatter(ra_flipped_2, dec_flipped_2, color='red')
    
    plt.xlabel('Right Ascension (degrees)')
    plt.ylabel('Declination (degrees)')
    plt.title('Random Points from Flipped Caps')
    plt.legend()
    plt.show()

def plot_combined():
    minter = pymangle.Mangle('intersection.ply')
    mflip1 = pymangle.Mangle('flipped1.ply') 
    mflip2 = pymangle.Mangle('mflip2.ply')
    ra_minter, dec_minter = minter.genrand(1000)
    ra_mflip1, dec_mflip1 = mflip1.genrand(1000)
    ra_mflip2, dec_mflip2 = mflip2.genrand(1000)

    plt.figure(figsize=(10, 6))
    plt.scatter(ra_minter, dec_minter, color='blue')
    plt.scatter(ra_mflip1, dec_mflip1, color='green')
    plt.scatter(ra_mflip2, dec_mflip2, color='red')
    plt.xlabel('Right Ascension (degrees)')
    plt.ylabel('Declination (degrees)')
    plt.legend()
    plt.show()

def plot_neg_constraints():
    mneg = pymangle.Mangle('intersection_neg.ply')
    ra_mneg, dec_mneg = mneg.genrand(1000000)
    plt.figure(figsize=(10, 6))
    plt.scatter(ra_mneg, dec_mneg, s=1, alpha=0.1)
    plt.xlabel('Right Ascension (degrees)')
    plt.ylabel('Declination (degrees)')
    plt.show()





# MAIN
###############################
###############################
if __name__ == '__main__':
  cap1 = cap_with_radius(76, 36, 5)
  cap2 = cap_with_radius(75, 35, 5)
  cap3 = cap_with_flipped_radius(76, 36, -5)

  write_to_mangle('intersection.ply', [cap1, cap2], single_polygon=True)
  write_to_mangle('bothcaps.ply', [cap1, cap2], single_polygon=False)


  write_to_mangle('flipped1.ply', [cap3, cap2], single_polygon=True)
  write_to_mangle('flipped2.ply', [cap3, cap2], single_polygon=False)


  generate_and_plot_random_points()
  flipped()


  # Part 5 and 6
  cap2_flipped = cap_with_flipped_radius(75, 35, 5) 
  write_to_mangle('mflip2.ply', [cap1, cap2_flipped], single_polygon=True)
  plot_combined()

  cap1_neg = cap_with_flipped_radius(76, 36, 5)  
  cap2_neg = cap_with_flipped_radius(75, 35, 5)  
  write_to_mangle('intersection_neg.ply', [cap1_neg, cap2_neg], single_polygon=True)

  plot_neg_constraints()

