# M. Rizzo Smith
# v1 2/20/24
# Week 7 Lecture 1 Tasks
# ASTR 8080 Masking


# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
import astropy 
from matplotlib import rc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic
from astropy.coordinates import get_body
from astropy.time import Time
import sys

sys.path.insert(0, '../week6/')
import spherical_caps
from spherical_caps import ra_cap, dec_cap


# FUNCTIONS
###############################
###############################
def polygon(ra1, ra2, dec1, dec2, weight, output_name):
    caps = []
    ra1_cap = (ra_cap(ra1))
    ra2_cap = (ra_cap(ra2))
    dec1_cap = (dec_cap(dec1))
    dec2_cap = (dec_cap(dec2))

    ra2_cap[-1] *= -1
    dec2_cap[-1] *= -1
    caps.append(ra1_cap)
    caps.append(ra2_cap)    
    caps.append(dec1_cap)
    caps.append(dec2_cap)
    caps_str= str(np.array(caps))
    caps_str = caps_str.replace(',', ' ')
    caps_str = caps_str.replace('[', ' ') 
    caps_str = caps_str.replace(']', '')
    h = np.sin(dec2 * (np.pi/180)) - np.sin(dec1*(np.pi/180))
    # MRS Calculate the width of the triangle
    w = ((ra2*15) * (np.pi/180)) - ((ra1*15) * (np.pi/180))
    # MRS Combine both for the area and convert to square degrees
    area = h*w
   
    f = open(output_name, 'w')
    f.write('1 polygons\n')
    f.write(f'polygon 1 ( 4 caps, {weight} weight, 0 pixel, {area:.6f} str):\n' )
    f.write(caps_str)
    f.close()

    return area

def add_polygon(main, added):
    with open(added, 'r') as file2:
        lines = file2.readlines()

        lines[1] = lines[1].replace('polygon 1', 'polygon 2')
    with open(main, 'r+') as f1:
        f1.seek(0)
        f1_line = f1.readlines()
        f1_line[0] = f1_line[0].replace('1 polygons', '2 polygons')
        f1.seek(0)
        f1.write(f1_line[0])

    with open(main, 'a') as file1:
        file1.write('\n')
        for line in lines[1:]:
            file1.write(line)

    return

# MAIN
###############################
###############################
if __name__ == '__main__':
    polygon1 = polygon(5, 6, 30, 40, 0.9, 'caps.ply')
    polygon2 = polygon(10, 12, 60, 70, 0.2, 'poly2.ply')
    main_test = polygon(5, 6, 30, 40, 0.9, 'main.ply')
    add_polygon('main.ply', 'poly2.ply')
