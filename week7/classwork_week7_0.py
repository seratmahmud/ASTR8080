# S. Saad
# ASTR 8080 week7, classwork0

# IMPORT BLOCK
###############################
###############################


from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
import pymangle





# FUNCTIONS
###############################
###############################

def cap_ra(ra):
    x, y, z = sph_to_cart(ra + 90, 0)  
    h = 1.0  
    return [x, y, z, 1 - h]

def cap_dec(dec):
    x, y, z = sph_to_cart(0, 90)
    h = np.sin(np.radians(dec))
    return [x, y, z, 1 - h]

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

def calculate_area(ra_start, ra_end, dec_start, dec_end):
    ra_start_rad = np.radians(ra_start)
    ra_end_rad = np.radians(ra_end)
    dec_start_rad = np.sin(np.radians(dec_start))
    dec_end_rad = np.sin(np.radians(dec_end))
    return (ra_end_rad - ra_start_rad) * (dec_end_rad - dec_start_rad)

def write_to_mangle(filename, polygons):
    with open(filename, 'w') as f:
        f.write(f'{len(polygons)} polygons\n')
        for i, polygon in enumerate(polygons, start=1):
            caps = polygon['caps']
            weight = polygon['weight']
            area = polygon['area']
            f.write(f'polygon {i} ( {len(caps)} caps, {weight} weight, 0 pixel, {area} str):\n')
            for cap in caps:
                f.write(' '.join(map(str, cap)) + '\n')

def main():
    ra_start_deg = 5 * 15 
    ra_end_deg = 6 * 15 
    dec_start_deg = 30
    dec_end_deg = 40
    area = calculate_area(ra_start_deg, ra_end_deg, dec_start_deg, dec_end_deg)
    weight = 0.9
    cap1 = cap_with_radius(ra_start_deg, dec_start_deg, 5)
    cap2 = cap_with_radius(ra_end_deg, dec_end_deg, 5)

    ra2_start_deg = 10 * 15 
    ra2_end_deg = 12 * 15 
    dec2_start_deg = 60
    dec2_end_deg = 70
    area2 = calculate_area(ra2_start_deg, ra2_end_deg, dec2_start_deg, dec2_end_deg)
    weight2 = 0.2
    cap3 = cap_with_radius(ra2_start_deg, dec2_start_deg, 5)
    cap4 = cap_with_radius(ra2_end_deg, dec2_end_deg, 5)

    polygons = [{'caps': [cap1, cap2], 'area': area, 'weight': weight},
        {'caps': [cap3, cap4], 'area': area2, 'weight': weight2}]

    write_to_mangle('file.ply', polygons)
    
    
    n = 1000000
    ra = np.random.uniform(0, 360, n)
    cos_dec = np.random.uniform(-1, 1, n)
    dec = np.degrees(np.arcsin(cos_dec))
    print(ra, np.max(dec))
    
    mask = pymangle.Mangle("file.ply")
    in_mask = mask.contains(ra, dec)
    ra_in_mask = ra[in_mask]
    dec_in_mask = dec[in_mask]
    
    plt.scatter(ra_in_mask, dec_in_mask, color='red')
    plt.scatter(ra, dec, color='blue')
    plt.show()










# MAIN
###############################
###############################
if __name__ == '__main__':
    main()








