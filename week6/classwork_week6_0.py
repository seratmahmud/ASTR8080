# S. Saad
# ASTR 8080 week6, classwork1

# IMPORT BLOCK
###############################
###############################


from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np





# FUNCTIONS
###############################
###############################



def sph_to_cart(ra_deg, dec_deg):
    c = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree, frame='icrs')
    x = np.cos(c.ra.radian) * np.cos(c.dec.radian)
    y = np.sin(c.ra.radian) * np.cos(c.dec.radian)
    z = np.sin(c.dec.radian)
    return x, y, z


def cap_ra(ra):
    x, y, z = sph_to_cart(ra + 90, 0)  
    h = 1.0  
    return [x, y, z, 1 - h]


def cap_dec(dec):
    x, y, z = sph_to_cart(0, 90)
    h = np.sin(np.radians(dec))
    return [x, y, z, 1 - h]


def cap_with_radius(ra, dec, radius):
    x, y, z = sph_to_cart(ra, dec)
    one_minus_h = 1 - np.cos(np.radians(radius))
    return [x, y, z, one_minus_h]


def write_spherical_caps_to_file(filename, caps):
    hdr = '1 polygons\npolygon 1 ( 3 caps, 1 weight, 0 pixel, 0 str):'
    np.savetxt(filename, caps, fmt='%1.9f', header=hdr, comments='', newline='\n ')




# MAIN
###############################
###############################
if __name__ == '__main__':


    cap_ra_vector = cap_ra(5 * 15) 
    cap_dec_vector = cap_dec(36)  
    cap_with_radius_vector = cap_with_radius(5 * 15, 36, 1) 
    
    print(cap_ra_vector, cap_dec_vector, cap_with_radius_vector)
    
    caps = np.array([cap_ra_vector, cap_dec_vector, cap_with_radius_vector])
    
    filename = 'spherical_caps.txt'
    write_spherical_caps_to_file(filename, caps)






