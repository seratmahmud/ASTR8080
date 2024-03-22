# M. Kaldor
# v1 1/25/2024
# ASTR 8080 HW0

# to import this from another directory:
# import sys
# sys.path.insert(0, '../week1/')
# import hw0
# from hw0 import function

# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt
from matplotlib import rc
import time
import pdb
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import sfdmap
import extinction
from numpy.random import random
import healpy as hp


# FUNCTIONS
###############################
###############################
def func1():
    """
    This is a function that

    Args:

    Returns:

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE:
    #
    #
    # CALLING SEQUENCE:
    #
    #
    # INPUTS:
    #
    # -------------------------------------------------------------
    #

###############################
def func2():
    """
    This is a function

    Args:

    Returns:

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE:
    #
    #
    # CALLING SEQUENCE:
    #
    #
    # INPUTS:
    #
    # -------------------------------------------------------------
    #


# MAIN
###############################
###############################
if __name__ == '__main__':
    # MEK set the timer
    time0 = time.time()

    # phi = ra in rad
    # theta = 90-dec in rad

    ra = 360. * (random(1000000))
    dec = (180 / np.pi) * np.arcsin(1. - random(1000000) * 2.)

    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="aitoff")
    ax.scatter(ra, dec, marker="o", color="red", s=0.01, alpha=0.15)
    ax.grid(color="blue", linestyle="dashed", linewidth=2.5)
    xlab = ["14h", "16h", "18h", "20h", "22h", "0h", "2h", "4h", "6h", "8h", "10h"]
    ax.set_xticklabels(xlab, weight=800)
    plt.title("Aitoff Projection")
    # plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="lambert")
    ax.scatter(ra, dec, marker="o", color="red", s=0.01, alpha=0.15)
    ax.grid(color="blue", linestyle="dashed", linewidth=2.5)
    xlab = ["14h", "16h", "18h", "20h", "22h", "0h", "2h", "4h", "6h", "8h", "10h"]
    ax.set_xticklabels(xlab, weight=800)
    plt.title("Lambert Projection")
    # plt.show()
    '''

    phi = ra*np.pi/180
    theta = (np.pi/2) - dec*np.pi/180
    nside = 1
    pix = hp.ang2pix(nside, theta, phi)
    npix = 12*nside**1

    ASr = 4*np.pi
    acalcrad = ASr/npix
    Adeg = 4*np.pi*(180/np.pi)**2
    acalcdeg = Adeg/npix
    print("I calculated", acalcrad, "Sr")
    print("I calculated", acalcdeg, "sq deg")

    arearad = hp.nside2pixarea(1)
    print("Healpy calculated", arearad, "Sr")
    areadeg = hp.nside2pixarea(1, degrees=True)
    print("Healpy calculated", areadeg, "sq deg")

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(pix, bins=np.arange(13) - 0.5, color="blue", align="mid", edgecolor='black', linewidth=1.2)
    ax.set_xlabel("Value", fontsize=22)
    ax.set_ylabel("Number", fontsize=22)
    fig.tight_layout()
    # plt.show()
    plt.close(fig)

    w2 = np.where(pix == 2)
    w5 = np.where(pix == 5)
    w8 = np.where(pix == 8)

    plt.scatter(phi/np.pi, np.cos(theta), s=0.01, alpha=0.1, color="blue")
    plt.scatter(phi[w2] / np.pi, np.cos(theta[w2]), s=0.01, alpha=0.15, color="red")
    plt.scatter(phi[w5] / np.pi, np.cos(theta[w5]), s=0.01, alpha=0.15, color="yellow")
    plt.scatter(phi[w8] / np.pi, np.cos(theta[w8]), s=0.01, alpha=0.15, color="green")
    plt.show()

    nside2 = 2
    pix2 = hp.ang2pix(nside2, theta, phi)
    w2_5 = np.where(pix2 == w5)
    print(np.shape(w2_5))

    plt.scatter(phi / np.pi, np.cos(theta), s=0.01, alpha=0.1, color="blue")
    plt.scatter(phi[w5] / np.pi, np.cos(theta[w5]), s=0.01, alpha=0.15, color="yellow")
    plt.scatter(phi[w2_5] / np.pi, np.cos(theta[w2_5]), s=0.01, alpha=0.15, color="red")
    plt.show()


    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))
