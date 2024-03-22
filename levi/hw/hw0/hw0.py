# L. Schult
# v1 24 jan 2024
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
#import astropy
import pdb
import time
import argparse


# LSS argument parsing from command line:
parser = argparse.ArgumentParser(description= "parsing slope and intercept vals\
                                 from command line for scattering + fitting")

parser.add_argument('--slopeval', required=True, type=float, help = 'slope for line to fit \
                    for.')
parser.add_argument('--interceptval', required=True, type=float, help = 'y-intercept for \
                    line that will be fit')
parser.add_argument('--numpoints', required=False, default=10, type=int, help= 'number of points \
                    generated from given line + scattered for a new linear fit.\
                    default is 10.')

args = parser.parse_args()



# FUNCTIONS
###############################
###############################
def genpointsfromline(m, b, npts=10):
    """
    This is a function that generates 10 random points along a given y=mx+b line
    the function scatters the points around the exact y values for randomly 
    generated x values. The scatter is a random draw from a gaussian centered at
    the y point and has a std dev of 0.5.

    It returns the x values, the scattered y values, and the uncertainty on each
    y value (0.5).

    Args:
        m (float): slope of the line the points will be generated from
        b (float): y intercept of the line the points will be generated from
        npts (int): [OPTIONAL] number of points generated. Default is 10.
    Returns:
        x (ndarray): An array of the points' x vals
        y (ndarray): an array of y vals for given x (NO SCATTER)
        yscat (ndarray): an array of the points' y vals (with scatter)
        y_unc (ndarray): an array of the points' uncertainty in y (0.5) 
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   generate points according to a given line (from slope + intercept values)
#   with scatter added to the y values according to a gaussian with stddev (0.5)
#   returns the x values, scattered y values, and y uncertainty
# CALLING SEQUENCE:
#   xvals, yvals, yscat, yscatuncert = genpointsfromline(m, b) 
#
# INPUTS:
#   m (float) -- slope to generate points from
#   b (float) -- y intercept of line points are generated from.
#-
#-------------------------------------------------------------  
    starttime = time.time()
    x = np.random.uniform(0, 10, size=npts) # LSS generating random points
    y = m*x + b 
    distG_pts = np.random.normal(0, 0.5, size=npts) # LSS gaussian to make noise
    y_scat = y + distG_pts # LSS adding noise to generated values
    y_scat_unc = np.zeros_like(y_scat) # LSS filling array with same shape with 0.5
    y_scat_unc.fill(0.5)
    print(f'points generated from input line: y = {m}x+{b}')
    endtime = time.time()
    totalruntime = endtime - starttime
    print(f'total time to generate points: {totalruntime} s')
    return x, y, y_scat, y_scat_unc


def linfit_and_plot(x, y, y_unc, m, b):
    """
    This is a function that makes a linear fit to input points x and y and plots\
    the points with given uncertainties, the linear fit to the points, and the\
    line used to generate the points.

    Args:
        x (ndarray): A series of x values
        y (ndarray): y values for each x value, with a hopefully linearish \
        relationship for fitting using numpy.polyfit
        y_unc (ndarray): uncertainties for each of the y values
        m (float): the slope of the line used to generate the data. 
        b (float): the y-intercept of the line used to generate the data.
    Returns:
        m_fit (float): The slope of the line fitted to the input data.
        b_fit (float): The y-intercept of the line fitted to the input data.
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   fits a line to input data, plots the data, line fit to it, and the line that\
#   the data was generated from. The original line has to be supplied to the 
#   function through the m and b arguments.
#
# CALLING SEQUENCE:
#   slope_fit, yint_fit = linfit_and_plot(x, y, y_unc, m, b) 
#
# INPUTS:
#   x, y, y_unc, m, b
#-
#-------------------------------------------------------------  
    starttime = time.time()
    # LSS fitting the yvals using linear polyfit
    m_fit, b_fit = np.polyfit(x, y, 1)
    print('line fit to given points')
    # LSS generating exact points from given original line
    yexact = m*x+b
    fitendtime = time.time()
    print(f'time to fit line to data: {fitendtime-starttime} s')

    # LSS plotting
    print('plotting beginning')
    plt.plot(x, yexact, color='C2', label='original line')
    plt.errorbar(x, y, yerr=y_unc, label='data', fmt='.') # LSS data
    plt.plot(x, m_fit*x+b_fit, color='C1', label='polyfit line'\
             ) # LSS best fit line from data
    plt.xlabel('x', weight='bold')
    plt.ylabel('y', weight='bold')
    plt.title('Line Fit To Given Points', weight='bold')
    plt.legend()
    plotendtime = time.time() # LSS function timing + feedback.
    print(f'time to plot: {plotendtime-fitendtime} s')
    # LSS saving figure + displaying
    plt.savefig('./fitlinetopoints.png', format='png', dpi=400, \
                bbox_inches='tight')
    print('file [fitlinetopoints.png] saved in script\'s current directory')
    print('to finish the script, close the plot that has popped up.')
    plt.show()
    print('plotting ending')
    
    return m_fit, b_fit


# MAIN
###############################
###############################
if __name__ == '__main__':
    fullstarttime = time.time()
    xvals, yvals, yscat, yscatunc = genpointsfromline(args.slopeval, \
                                                      args.interceptval, npts=\
                                                        args.numpoints)
    fitslope, fitint = linfit_and_plot(xvals, yscat, yscatunc, args.slopeval,\
                                       args.interceptval)
    fullendtime = time.time()
    print(f'total time for whole program: {fullendtime-fullstarttime} s')
