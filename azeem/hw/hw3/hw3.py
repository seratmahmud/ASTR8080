# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
import time
from astropy import units as u
from astropy.coordinates import SkyCoord
import pymangle
import sys
sys.path.insert(0,'../../week6/')
from sphericalcaps import ra_bound, dec_bound, circle_cap
sys.path.insert(0,'../hw2/')
from HW2 import random_pop

def four_poly(ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4):
    """
    This is a function that makes a "lat-lon" rectangle given 2 ra bounded and
    2 dec bounded caps.

    Args:
        ra1(float): A number that represents the RA of the first RA cap.
        dec1 (float): The dec that correponds to the first RA cap.
        ra2 (float): A number that represents the RA of the second RA cap.
        dec2 (float): The dec that correponds to the first RA cap.
        ra3 (float): A number that represents the RA of the first dec cap.
        dec3 (float): The dec that correponds to the first dec cap.
        ra4 (float): A number that represents the RA of the second dec cap.
        dec4 (float): The dec that correponds to the second dec cap.
    Returns:
        racap1 (array): An output numpy array of all the random RA's in the box.
        racap2 (array): An output numpy array of all the random decs in the box.
        deccap1 (array): An output numpy array of all the random decs in the box.
        deccap2 (array): An output numpy array of all the random decs in the box.
    """

    #-------------------------------------------------------------
    #+
    # PURPOSE:
    #   Gives area for a "lat-long" rectangle given 4 caps that are RA and dec
    #   bounded.
    #
    #
    #
    # CALLING SEQUENCE:
    #   four_poly(ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4)
    #
    # INPUTS:
    #   ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4 - RAs and decs for each of the 4 caps
    #-
    #-------------------------------------------------------------

    #AB call ra and dec cap functions from week 6 Tasks
    #These make ra and dec bounded caps given an RA and dec
    racap1 = ra_bound(ra1, dec1)
    racap2 = ra_bound(ra2, dec2)
    deccap1 = dec_bound(ra3, dec3)
    deccap2 = dec_bound(ra4,dec4)
    #AB flip last term in the second ra and dec caps to make sure we fill in the
    #shape
    racap2[-1] = -1*racap2[-1]
    deccap2[-1] = -1*deccap2[-1]

    return(racap1,racap2,deccap1,deccap2)

def write_ply(caps):
    """
    This is a function that writes a polygon file corresponding to a lat-lon
    rectangle.

    Args:
        caps (array): An array that correponds to spherical and RA and dec bounded
            caps.
    Returns:
        Creates a .ply file
    """

    #-------------------------------------------------------------
    #+
    # PURPOSE:
    #   Take array of caps and make a .ply file.
    #
    #
    #
    # CALLING SEQUENCE:
    #   write_ply(caps)
    #
    # INPUTS:
    #   caps - array of spherical caps and caps in the "lat-lon" rectangle.
    #-
    #-------------------------------------------------------------

    #AB create polygon file and specify caps
    f = open("hw3_caps.ply", "w")
    area = 1 #area is 1
    latlon = caps[0:4] #first 4 caps from lat-lon rectangle
    circles = caps[4:]

    #AB write .ply file
    f.write(f"{len(circles)} polygons\n")
    for i in range(0, len(circles)):
        f.write(f" polygon {i+1} ( 5 caps, 1.0 weight, 0 pixel, {area} str):\n")
        for cap in latlon:
            f.write("  " + str(cap[0]) + " " + str(cap[1]) + " " + str(cap[2]) + " " + str(cap[3]) + "\n")
        f.write("  " + str(circles[i][0]) + " " + str(circles[i][1]) + " " + str(circles[i][2]) + " " + str(circles[i][3]) + "\n")
    f.close()


def rand_area(ra1,dec1,ra2,dec2):
    """
    This is a function that finds the area of a square field in sq. degrees, while
    also finding the random RA and dec postions that populate the area.

    Args:
        ra1 (f1loat): A number that represents the lowest RA of your box.
        dec1 (float): A number that represents the lowest dec of your box.
        ra2 (float): A number that represents the highest RA of your box.
        dec2 (float): A number that represents the highest dec of your box.
    Returns:
        area (float): The area of your square as a function of the inputs
        ra_rand (array): An output numpy array of all the random RA's in the box.
        dec_rand (array): An output numpy array of all the random decs in the box.
    """

    #-------------------------------------------------------------
    #+
    # PURPOSE:
    #   Find area of "lat-lon" rectangle in sq. degrees
    #
    #
    #
    # CALLING SEQUENCE:
    #   rand_area(ra1,dec1,ra2,dec2)
    #
    # INPUTS:
    #   ra1,dec1,ra2,dec2 - RAs and decs of the bounded region.
    #-
    #-------------------------------------------------------------

    #AB first we call the random_pop function from HW2. This randomly populates
    #a square field given the bounds and returns an array of the RAs and decs
    ra_rand, dec_rand = random_pop(ra1,ra2,dec1,dec2)
    #AB formula for area of lat-lon rectangle bounded by ra and dec, in square degrees
    area = (ra2*(np.pi/180)-(np.pi/180)*ra1)*(np.sin(dec2*(np.pi/180))-np.sin((np.pi/180)*dec1))*(180/np.pi)**2

    return area, ra_rand, dec_rand

def pop_area(ra1, ra2, dec1, dec2, ply, quasars):
    """
    This is a function that finds the area of the mask from the .ply file,
    and plots quasar positions.

    Args:
        ra1 (float): A number that represents the lowest RA of your box.
        ra2 (float): A number that represents the highest RA of your box.
        dec1 (float): A number that represents the lowest dec of your box.
        dec2 (float): A number that represents the highest dec of your box.
        ply (string): File path of polygon file.
        quasars (float): File path of quasar file.
    Returns:
        N/A just plots
    """

    #-------------------------------------------------------------
    #+
    # PURPOSE:
    #   Gives area of mask and plots quasar positions
    #
    #
    #
    # CALLING SEQUENCE:
    #   pop_area(ra1, ra2, dec1, dec2, ply, quasars)
    #
    # INPUTS:
    #   ra1, ra2, dec1 ,dec2 - RAs and decs of the bounded rectangle.
    #   ply, quasars - file paths for the .ply file and quasar file
    #-
    #-------------------------------------------------------------
    #AB set the timer
    time0  = time.time()

    #AB use Mangle to make mask from .ply file
    m = pymangle.Mangle("hw3_caps.ply")
    #AB get the random RAs and decs, as well as the area of the rectangle.
    area, ra_rand, dec_rand = rand_area(ra1, dec1, ra2, dec2)

    #AB make empty list and determine number of points
    temp = []
    num_points = len(ra_rand)
    #AB loop through the number of points and determine the precision of calculation
    #If precesion is not met, repeat.
    while num_points < 10000000:
        for i in range(0, num_points):
            #AB determine which points are within mask
            good = m.contains(ra_rand, dec_rand)
            #AB to find area of the mask, use a ratio of areas to points
            area_mask = (area / num_points) * len(ra_rand[good])
            temp.append(area_mask)
            i += 1
        temp = np.array(temp)
        std = np.std(temp)
        #AB check for precision of 0.5 sq degrees
        if std < 0.5:
            print('The area of the mask in sq. degrees: ', np.mean(temp))
            break
        else:
            num_points += 10000000
            i = 0
            temp = []
    time1  = time.time()
    print("    Finding best area of mask: {0:.3f}".format(time1-time0))

    #AB Now for the quasars
    #AB makes empty array to store quasar positions

    time2  = time.time()

    all_pos = []

    #AB To open file and make an array:
    with open(quasars, 'r') as file:
        for line in file:
            qpos = line
            qpos = f'{qpos[:2]} {qpos[2:4]} {qpos[4:9]} {qpos[9:12]} {qpos[12:14]} {qpos[14:]}'
            cquasar = SkyCoord(qpos, unit=(u.hourangle, u.deg)) #converts to SkyCoord object
            all_pos.append(cquasar) #make array
    all_pos = np.array(all_pos)

    #AB make RAs and decs of quasars into SkyCoord object
    c = SkyCoord(all_pos, unit=(u.hourangle, u.deg))

    #AB determine which quasars are in the mask
    q_match = m.contains(c.ra.degree, c.dec.degree)
    #AB number density of quasars and plot
    n_dens = len(c.ra[q_match]) / area_mask
    print("    matching the quasar postions and calculating number density: {0:.3f}".format(1e6*(time2-time1)))
    plt.scatter(c.ra.degree, c.dec.degree, label=f'Mask Area: {area_mask:.2f} sq degs')
    plt.scatter(c.ra.degree[q_match],c.dec.degree[q_match], label=f'Number Density: {n_dens:.2f} per sq deg')
    plt.xlabel('RA')
    plt.ylabel('dec')
    plt.legend()
    time3  = time.time()
    print("    plot [s]: {0:.3f}".format((time3-time2)))
    plt.savefig('quasars.png', dpi=300)


if __name__ == '__main__':
    #AB set the timer
    time0  = time.time()
    #1
    #AB input RAs and decs, need to convert from HMS to degrees
    ra1 = 15*(10+(15/60))
    ra2 = 15*(11+(15/60))
    dec1 = 30
    dec2 = 40
    latlon = np.array(four_poly(ra1,0,ra2,0,0,dec1,0,dec2)) #For RA cap, dec is 0 and vice versa

    #AB Use circle_cap function from week 6 tasks to find "spherical circle" on sky
    #need RA, dec, and search radius, all in degrees
    cap1 = circle_cap(155,34,2)
    cap2 = circle_cap(159,36,2)
    cap3 = circle_cap(163,34,2)
    cap4 = circle_cap(167,36,2)

    #AB make into a single array
    allcaps = np.array([cap1, cap2, cap3, cap4])
    rectangle = np.concatenate((latlon, allcaps), axis=0)
    write_ply(rectangle) #write the .ply file
    time1  = time.time()
    print("    Write polygon file: {0:.3f}".format(time1-time0))

    #2 ans 3
    #AB get area of mask and plot quasars
    pop_area(ra1, ra2, dec1, dec2, 'hw3_caps.ply', '../../../runnoe/week7/HW3quasarfile.dat')
