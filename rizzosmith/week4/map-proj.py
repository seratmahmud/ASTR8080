# M. Rizzo Smith
# v1 1/25/24
# ASTR 8080 HW1


# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic
from astropy.coordinates import get_body
from astropy.time import Time
from numpy.random import random
# FUNCTIONS
###############################
###############################

    

# MAIN
###############################
###############################
if __name__ == '__main__':

    ra = 2 * np.pi * (random(10000)-0.5)
    dec = np.arcsin(1.-random(10000)*2.)

   # plt.scatter(np.degrees(ra), np.degrees(dec), s=1)
  #  plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='aitoff')
    ax.scatter(ra, dec, s=1, color='blue')
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'] 
    ax.set_xticklabels(xlab, weight=800)
    ax.grid(color='black', linestyle='solid', linewidth=1.5)
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='lambert')
    ax.scatter(ra, dec, s=1, color='blue')
    xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'] 
    ax.set_xticklabels(xlab, weight=800)
    ax.grid(color='black', linestyle='solid', linewidth=1.5)
    plt.show()

    big_ra = np.append(ra1, ra2)
    big_dec = np.append(dec1, dec2)
