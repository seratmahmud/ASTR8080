# S. Saad
# ASTR 8080 week3, classwork0

# IMPORT BLOCK
###############################
###############################
from numpy.random import random
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord, search_around_sky, Angle


# FUNCTIONS
###############################
###############################


def main():
    
    # Part 1
    coord1 = SkyCoord(ra=263.75*u.degree, dec=-12.9*u.degree, frame='icrs')
    coord2 = SkyCoord(ra="20h24m59.9s", dec="10d06m00s", frame='icrs')

    cart1 = coord1.cartesian
    cart2 = coord2.cartesian

    
    dot_product = np.dot(cart1.xyz, cart2.xyz)
    norm_cart1 = np.linalg.norm(cart1.xyz)
    norm_cart2 = np.linalg.norm(cart2.xyz)
    cos_angle = dot_product / (norm_cart1 * norm_cart2)
    angle_rad = np.arccos(cos_angle)
    angle_deg = Angle(np.degrees(angle_rad))

    angle_separation = Angle(coord1.separation(coord2))
    
    print("Both way of calculating separation angle worked:",\
          angle_deg, angle_separation)
    
    
    # Part 2
    alpha_min = 2*15
    alpha_max = 3*15
    delta_min = -2
    delta_max = 2
    
    alpha_1 = np.random.uniform(alpha_min, alpha_max, 100)
    delta_1 = np.random.uniform(delta_min, delta_max, 100)
    alpha_2 = np.random.uniform(alpha_min, alpha_max, 100)
    delta_2 = np.random.uniform(delta_min, delta_max, 100)
    
    # Part 3
    set_1 = SkyCoord(ra=alpha_1*u.deg, dec=delta_1*u.deg)
    set_2 = SkyCoord(ra=alpha_2*u.deg, dec=delta_2*u.deg)
    
    id1, id2, d2d, d3d = set_2.search_around_sky(set_1, seplimit=(10/60)*u.deg)
    #print(set_1[id1], set_2[id2])
    
    plt.figure(figsize=(12, 8))
    plt.scatter(alpha_1, delta_1)
    plt.scatter(alpha_2, delta_2)
    #plt.scatter(alpha_1[id1], delta_1[id1], marker="+")
    #plt.scatter(alpha_2[id2], delta_2[id2], marker="+")
    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")
    plt.grid(True)
    plt.show()
    
    
    # Part 4
    new_ra = np.append(alpha_1[id1], alpha_2[id2])
    new_dec = np.append(delta_1[id1], delta_2[id2])
    
    plt.figure(figsize=(12, 8))
    plt.scatter(new_ra, new_dec)
    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")
    plt.grid(True)
    plt.show()
    
    # Part 5
    coord_star = SkyCoord(ra="02h20m05s", dec="00d06m12s", frame='icrs')
    seps_1 = coord_star.separation(set_1)
    seps_2 = coord_star.separation(set_2)
    mask_1 = np.where(seps_1 < 1 * u.deg)
    mask_2 = np.where(seps_2 < 1 * u.deg)
    
    plt.figure(figsize=(10, 5))
    plt.scatter(set_1.ra, set_1.dec, color='blue')
    plt.scatter(set_2.ra, set_2.dec, color='blue')
    plt.scatter(set_1[mask_1].ra, set_1[mask_1].dec, color='red')
    plt.scatter(set_2[mask_2].ra, set_2[mask_2].dec, color='red')
    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")
    plt.grid(True)
    plt.show()
    #print("Set 1 stars:", set_1[mask_1])
    #print("Set 2 stars:", set_2[mask_2])
    
    
    
    

    
# MAIN
###############################
###############################
if __name__ == '__main__':
    main()

