# M. Kaldor
# v1 1/21/2024
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


# FUNCTIONS
###############################
###############################
def gauss_line(m, b, u, l):
    """
    This is a function that selects 10 random points between 1-10, creates a line with the form y=mx+b from these x values,
    adds Gaussian noise with a standard deviation of 0.5 to these points, and returns the resulting y values.

    Args: m = [float]
        b = [float]
        u = [float]
        l = [integer]

    Returns: three arrays: x values, y values with noise, uncertainty on the y values

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: creating a linear relationship with Gaussian noise
    #
    #
    # CALLING SEQUENCE: gauss_line(1, 2, 3, 4)
    #
    #
    # INPUTS:  m = slope of line
    #         b = y intercept of line
    #        u = uncertainty on y
    #       l = length of list
    #
    # -------------------------------------------------------------
    # MEK initiate x and y values
    xlist = 10*np.random.random(l)
    ylist = [m*x+b for x in xlist]
    # MEK add Gaussian noise to y values and save their uncertainty values
    randylist = [np.random.normal(y, u) for y in ylist]
    uncert_y = [u]*l
    # MEK plot the original line in blue and the data points and error bars in purple
    plt.plot((0, 10), (b, m*10), color="blue", label="Original line y="+str(m)+"x+"+str(b), linewidth=4)
    plt.scatter(xlist, randylist, color="purple", s=70, label="Simulated Data")
    plt.errorbar(xlist, randylist, yerr=uncert_y, ecolor="purple", linestyle="", capsize=4)
    return xlist, randylist, uncert_y

###############################
def polyfit(x, y, deg):
    """
    This is a function that fits a line to the data points that are fed in.

    Args: x = [list of floats]
        y = [list of floats]
        deg = [integer]

    Returns: 2d array - polynomial coefficients for an equation of a line (0th item=m=slope, 1st item=b=y intercept)

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: fitting a polynomial to any data produced
    #
    #
    # CALLING SEQUENCE: polyfit(1, 2, 3)
    #
    #
    # INPUTS:  x = x data
    #         y = y data
    #        deg = degree of polynomial for best fit
    #
    # -------------------------------------------------------------
    # MEK fit the polynomial
    fit = np.polyfit(x, y, deg)
    # MEK plot the line of best fit in red
    plt.plot((0, 10), (fit[1], fit[0] * 10), color="red",
             label="Line of Best Fit y="+str(round(fit[0], 3))+"x+"+str(round(fit[1], 3)), linewidth=4)
    return fit


# MAIN
###############################
###############################
if __name__ == '__main__':
    # MEK set the timer
    time0 = time.time()

    # MEK create figure
    plt.figure(figsize=(10, 8))

    # MEK creating Gaussian noise data
    data = gauss_line(1, 2, 0.5, 10)
    print("Simulated data: x =", data[0], "y =", data[1], "y uncertainty =", data[2], "\n")

    # MEK fit simulated data
    linefit = polyfit(data[0], data[1], 1)
    print("Line of best fit: m2 =", linefit[0], "b2 =", linefit[1], "\n")

    # MEK make plot pretty, save it, and display it
    plt.rcParams["font.family"] = "Times New Roman"
    plt.xlabel("X value", font="Times New Roman", fontsize=25)
    plt.ylabel("Y value", font="Times New Roman", fontsize=25)
    plt.title("Y value vs. X value", font="Times New Roman", fontsize=25)
    plt.tick_params(axis='both', labelsize=20)
    plt.legend(fontsize=20)
    print("Plotting figure...\n")
    savepath = "/Users/marykaldor/accdisk/astro8080/hw0.png"
    plt.savefig(savepath)
    print("Figure saved in", str(savepath), "\n")

    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))
    plt.show()