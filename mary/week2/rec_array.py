# M. Kaldor
# v1 1/16/2024
# ASTR 8080 Rec Array

# to import this from another directory:
# import sys
# sys.path.insert(0, '../week1/')
# import rec_array
# from rec_array import function

# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt
from matplotlib import rc
import pdb
import time


# FUNCTIONS
###############################
###############################
def function(x, params):
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
    #
    # -------------------------------------------------------------
    A, B, C = params
    y = A + B * x + C * x ** 2.
    return y


# MAIN
###############################
###############################
if __name__ == '__main__':
    # MEK set the timer
    time0 = time.time()

    # MEK read in file
    from astropy.io import fits
    fx = fits.open("/Users/marykaldor/accdisk/astro8080/struc.fits")
    print(fx.info())

    # MEK pull out data and headers
    objs = fx[1].data
    hdr = fx[1].header

    # MEK check values
    # print(np.shape(objs))
    # print(np.shape(hdr))
    # print(objs["RA"])
    # print(objs["DEC"])

    # MEK access first column of extinction header
    ext = objs["EXTINCTION"][:,0]
    print("ext", ext)

    # MEK pull indices for extinction>0.22
    ext22 = np.where(ext>0.22)
    print(ext22)

    # MEK plot them!
    # MEK plot RA against DEC
    plt.plot(objs["RA"], objs["DEC"], "bx")
    ra22 = [objs["RA"][idx] for idx in ext22]
    dec22 = [objs["DEC"][idx] for idx in ext22]
    plt.plot(ra22, dec22, "rx")
    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.title("DEC vs. RA")
    plt.show()

    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))