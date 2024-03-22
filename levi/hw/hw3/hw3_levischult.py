# L. Schult
# v1 
# ASTR 8080 
# to import this from another directory:
# import sys
# sys.path.insert(0, '../week1/')
# import polycalc
# from polycalc import get_poly_o3

# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import astropy
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
import pdb
import time
import pymangle
import os
import sys
sys.path.insert(0, '../hw2/')
from hw2_levischult import skyarea
import argparse


# FUNCTIONS
###############################
###############################
def testname(x):
    """
    This is a function that 

    Args:
        x (datatype): A...
    Returns:
        y (datatype): A ...
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   You pass butter
#
# CALLING SEQUENCE:
#   variable = testname(x) 
#
# INPUTS:
#   x
#-
#-------------------------------------------------------------  
    y     = x**2. 
    return y 

def latlonrect(ramin, ramax, decmin, decmax):
    """
    This is a function that makes 4 four vectors describing a rectangle in ra/dec
    useful for making mangle polygons
    all arguments should be given in degrees '13d'

    Args:
        ramin, ramax, decmin, decmax (str): 
        minimum/maximum right ascension and declination in format '20d' (for degrees)
    Returns:
        ralim1, ralim2, declim1, declim2 (arrays):
        len=4 vectors specifying the RA caps and DEC caps that specify the box in question
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   make four 4-vectors for mangle polygons that specify a ra/dec bounded rectangle
#
# CALLING SEQUENCE:
#   ra1, ra2, dec1, dec2 = latlonrect('120d', '130d', '20d', '30d')
#
# INPUTS:
#   ramin, ramax, decmin, decmax -- string in format for astropy Angle 
#-
#-------------------------------------------------------------  
    # LSS ra/dec lims in degrees
    ralim1 = racap(ramin)
    ralim2 = racap(ramax)
    # LSS this restricts to box rather than everything but box
    ralim2[-1] = ralim2[-1]*(-1) 
    declim1 = deccap(decmin)
    declim2 = deccap(decmax) 
    declim2[-1] = declim2[-1]*(-1)
    return ralim1, ralim2, declim1, declim2

def racap(rabound):
    """
    This is a function that gets the 4 vector in xyz(1-h)
    rabound should be '5h' string that is in hourangle
    center is found by using RA + 6 hours (90 deg)
    Args:
        rabound (str): a right ascencion bound in format '5h' or something that 
        astropy Angle can understand
    Returns:
        centervector (array): A four vector specifying the center of a cap and 
        its radius
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   make four vector for an RA cap - useful for mangle polygons
#
# CALLING SEQUENCE:
#   racapvector = racap('5h')
#
# INPUTS:
#   rabound - str in units that astropy angle can understand easily - '5h' or '75d'
#-
#-------------------------------------------------------------  
    racenter = Angle(rabound) + 6*u.hourangle # LSS get center of cap
    racenter = SkyCoord(racenter, 0, unit=u.hourangle, frame='icrs') # LSS make into skycoord
    centervector = [racenter.cartesian.x.value, racenter.cartesian.y.value, 
                    racenter.cartesian.z.value, 1]
    return centervector

def deccap(decbound):
    """
    This is a function that gets the 4 vector in xyz(1-h) for a cap bounded in declination
    decbound should be in '30d' in degrees. Note that center is set to +90deg, 
    so must alter function for south pole caps
    Args:
        rabound (str): A declination in '20d' format (something that astropy Angle
        can understand).
    Returns:
        centervector (array): An array of xyz(1-h) values corresponding to a DEC 
        northpole cap
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   make declination cap four vectors for mangle polygons
#
# CALLING SEQUENCE:
#   capvector = deccap('35d') 
#
# INPUTS:
#   decbound
#-
#-------------------------------------------------------------  
    decbound_ang = Angle(decbound) # LSS convert to a degree value
    deccenter = SkyCoord(0, 90, unit=[u.hourangle, u.degree], frame='icrs')
    centervector = [deccenter.cartesian.x.value, deccenter.cartesian.y.value, 
                     deccenter.cartesian.z.value, 1-np.sin(decbound_ang.radian)]
    return centervector

def anyply_writer(numpolygon, capvectors, polygonvects, filename):
    """
    This is a function that writes mangle polygon files (ending in .ply). can do
    any number of polygons with any number of vectors. the only limit is your 
    ability to index. Remember that polygons are AND for vectors within it and
    OR for multiple polygons in a ply file.

    Args:
        numpolygon (int): number of polygons
        capvectors (ndarray): should be n x 4 ndarray with cap 4-vectors
        see racap, deccap, radeccap
        polygonvects (list): len(polygonvects) = numpolygon. each element of
        the list is itself a list that has the index of the cap each polygon will take
        filename (str): polygon filename no ply needed
    Returns: NOTHING! YOU GET NOTHING!
        
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   write Mangle polygon files ending in .ply with arbitrary number of polygons 
#   and vectors. Does not specify weights for polygons
#
# CALLING SEQUENCE:
#   anyply_writer(2, [vect0, vect1, vect2], [[0, 1], [2]]) 
#
# INPUTS:
#   numpolygon, capvectors, polygonvects
#-
#-------------------------------------------------------------  
    with open(filename+'.ply', 'w') as outf:
        first_line = str(numpolygon) + ' polygons\n'
        outf.write(first_line)
        for pidx in range(numpolygon): # LSS taking care of 0 indexing
            polystr = f'polygon {pidx+1} ( {len(polygonvects[pidx])} caps, 1 weight, 0 pixel, 0 str):\n'
            outf.write(polystr) # LSS print first line with number of caps based
            # LSS on length of list of vector indices from polygonvects
            for cap_idx in polygonvects[pidx]:
                vectstr = str(capvectors[cap_idx]).replace('[', '')
                vectstr = vectstr.replace(']', '')
                vectstr = vectstr.replace(',', '')
                outf.write(vectstr+'\n')
        
        outf.close()
    return

def radeccap(ra, dec, theta):
    """
    This is a function that gets the 4 vector in xyz(1-h) for an arbitrary center
    format = '02h35m15s', '25d51m20s'
    theta in degrees
    Args:
        ra (str): A right ascension in '02h35m15s' format
        dec (str): a declination in '24d42m20s' format for skycoord object creation
        theta (int): a cap size in degrees
    Returns:
        centervector (array): A four vector specifying the center of the cap and cap radius
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   make four vector for Mangle polygons using ra/dec center and given cap size
#
# CALLING SEQUENCE:
#   capvector = radeccap('15h36m10s', '24d15m10s', 10) 
#
# INPUTS:
#   ra, dec, theta
#-
#-------------------------------------------------------------  
    capcenter = SkyCoord(ra, dec, unit=[u.hourangle, u.degree], frame='icrs')
    centervector = [capcenter.cartesian.x.value, capcenter.cartesian.y.value, 
                      capcenter.cartesian.z.value, 1-np.cos(np.deg2rad(theta))]
    return centervector


# MAIN
###############################
###############################
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--noplates', action='store_true', help='a way to turn off plate plots')
    parser.add_argument('--nomc', action='store_true', help='a way to turn off area estimation via mc integration')
    args = parser.parse_args()

    # LSS making rectangle on sky
    ralim1, ralim2, declim1, declim2 = latlonrect('153.75d', '168.75d', '30d', '40d')

    # LSS timing 
    starttime = time.time()

    # LSS plate definitions
    p1c = SkyCoord(155, 34, unit=u.degree, frame='icrs')
    p2c = SkyCoord(159, 36, unit=u.degree, frame='icrs')
    p3c = SkyCoord(163, 34, unit=u.degree, frame='icrs')
    p4c = SkyCoord(167, 36, unit=u.degree, frame='icrs')

    plate1 = radeccap(p1c.ra.to_string(u.hour), p1c.dec.to_string(u.degree), 2)
    plate2 = radeccap(p2c.ra.to_string(u.hour), p2c.dec.to_string(u.degree), 2)
    plate3 = radeccap(p3c.ra.to_string(u.hour), p3c.dec.to_string(u.degree), 2)
    plate4 = radeccap(p4c.ra.to_string(u.hour), p4c.dec.to_string(u.degree), 2)

    # LSS making array of 4vectors
    latlon_plates_vectors = [ralim1, ralim2, declim1, declim2, plate1, plate2,
                             plate3, plate4]

    # LSS writing polygons that are restricted to the plates within the rectangle
    anyply_writer(4, latlon_plates_vectors, [[0,1,2,3,4], [0,1,2,3,5], [0,1,2,3,6], [0,1,2,3,7]], 'platesonly')
    # LSS writing polygon that is rectangle for point generation/area estimate
    anyply_writer(1, latlon_plates_vectors, [[0,1,2,3]], 'rectonly')

    rectonly = pymangle.Mangle('rectonly.ply')
    platesonly = pymangle.Mangle('platesonly.ply')

    # LSS generating points for plotting the acceptable plates and area estimation
    Npts = 10000
    rect_ra_rand, rect_dec_rand = rectonly.genrand(Npts)
    plates_ra_rand, plates_dec_rand = platesonly.genrand(Npts)
    if not args.noplates:
        plt.scatter(plates_ra_rand, plates_dec_rand, color='C0', marker=',', s=1, label='plates')
        plt.xlabel('RA (º)')
        plt.ylabel('DEC (º)')
        plt.legend()
        plt.show()
    q1time = time.time()
    print('time for plate plotting + polygon writing (q1):', q1time-starttime)    

    # LSS plate area estimation
    if not args.nomc:
        rectarea = skyarea(153.75, 168.75, 30, 40)
        areastd = 1 # LSS initial values to get while loop running
        Npts=8000
        while areastd > 0.5:
            plateareaeval = [] # LSS only holds area evals for given Npts
            if Npts % 100 == 0: print(f'{Npts=}')
            for i in range(50):
            # LSS generating new points + seeing how many lie within plates
            # LSS doing this 50 times for every Npts to counteract random variations
                rect_ra_rand, rect_dec_rand = rectonly.genrand(Npts)
                # LSS find how many of the random points are within the plates
                maskedpoints = platesonly.contains(rect_ra_rand, rect_dec_rand)
                maskedpoints = np.sum([int(b) for b in maskedpoints])
                # LSS add to array for calculating std dev to see if accuracy threshold met
                plateareaeval.append(rectarea*(maskedpoints/Npts))
            areastd = np.std(plateareaeval)
            Npts += 2000

        platesarea = np.median(plateareaeval)
        print(f'area calculated with a std dev of {areastd} - plate area = {platesarea} sq deg')

        plt.hist(plateareaeval)
        plt.ylabel('counts')
        plt.xlabel('plate area (sqdeg)')
        plt.show()
    
    q2time = time.time()
    print('time for area estimation via random points (q2):', q2time-q1time)    

    # LSS Question 3
    # LSS loading in 17000 quasars
    with open('../../../runnoe/week7/HW3quasarfile.dat', 'r') as infile:
        rawtxt = np.loadtxt(infile, dtype=str)
    coordstrings = [] 
    for line in rawtxt:
        coordstrings.append(f'{line[:2]}h{line[2:4]}m{line[4:9]} {line[9:12]}d{line[12:14]}m{line[14:]}s')
    quasarcoords = SkyCoord(coordstrings)

    # LSS finding quasars within plates
    q_inplates = quasarcoords[platesonly.contains(quasarcoords.ra.deg, quasarcoords.dec.deg)]
    q_noplates = quasarcoords[~platesonly.contains(quasarcoords.ra.deg, quasarcoords.dec.deg)]

    # LSS calculate number density
    qin_n = q_inplates.shape[0] / platesarea
    print('quasar number density in plates =', qin_n)

    # LSS plotting the plates, quasars in plates and those outside of them
    # LSS plot is bounded by imaging survey region
    plt.scatter(plates_ra_rand, plates_dec_rand, color='C0', marker=',', s=1, label='plates')
    plt.scatter(q_inplates.ra.deg, q_inplates.dec.deg, color='k', marker='*', label='quasars in plates')
    plt.scatter(q_noplates.ra.deg, q_noplates.dec.deg, color='g', marker=',', label='quasars not in plates')
    plt.xlim(153.75, 168.75)
    plt.ylim(30, 40)
    plt.xlabel('RA (º)')
    plt.ylabel('DEC (º)')
    plt.title(f'Spectroscopic Survey area = {platesarea} sqdeg,\n quasars/sqdeg = {qin_n}')
    plt.legend()
    plt.savefig('spectsurvey_wquasars.png', format='png')
    plt.show()
    q3time = time.time()
    print('time for plotting quasars in/out of plates (q3):', q3time-q2time)    


    print('hello world XD')

