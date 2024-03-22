# J. Runnoe
# v1 8/26/2020
# v2 1/24/2022
# ASTR 8080 example for HW0

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
import pdb
import time

# FUNCTIONS
###############################
###############################
def get_poly_o2(x,params):
    """
    This is a function that calculates the polynomial A+Bx+Cx^2.

    Args:
        x (1darray): An input array of x values, or a single x value.
        params (4array): Input array of A,B,C constants.
    Returns:
        y (1darray): Output polynomial array with the same dimensions as x.
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   Return a 2nd order polynomial, which has the form:
#
#   y = A + B*x + C*x^2
#
# CALLING SEQUENCE:
#   y = get_poly_o2(x,params) 
#
# INPUTS:
#   params - [A,B,C]
#-
#-------------------------------------------------------------  
    A,B,C = params
    y     = A+B*x+C*x**2. 
    return y 


def get_poly_o3(x,params):
    """
    This is a function that calculates the polynomial A+Bx+Cx^2+Dx^3.

    Args:
        x (1darray): An input array of x values, or a single x value.
        params (4array): Input array of A,B,C,D constants.
    Returns:
        y (1darray): Output polynomial array with the same dimensions as x.
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   Return a 3rd order polynomial, which has the form:
#
#   y = A + B*x + C*x^2 + D*x^3
#
# CALLING SEQUENCE:
#   y = get_poly_o3(x,params) 
#
# INPUTS:
#   params - [A,B,C,D]
#-
#-------------------------------------------------------------  
    A,B,C,D = params
    y       = A+B*x+C*x**2.+D*x**3. 
    return y 

# MAIN
###############################
###############################
if __name__ == '__main__':
    # JCR set the timer
    time0  = time.time()

    # JCR print a test function
    print('polycalc(10,[0,1,2,3]) ==',get_poly_o3(10,[0.,1.,2.,3.]))

    # JCR define two polynomials over [0,100] 
    x  = np.arange(100)
    y  = get_poly_o2(x,[0,1,4])
    y2 = get_poly_o3(x,[0,1,2,3])
    
    # JCR time
    time1  = time.time()
    print("    polynomials calculated [mus]: {0:.3f}".format(1e6*(time1-time0)))

    # JCR plot them
    rc('text', usetex=True) # use latex
    fig, ax = plt.subplots(figsize=(4,3))
    ax.plot(x,y,'.',color='black',markersize=2,label=r'$x+2x^2+3x^3$')
    ax.plot(x,y2,'.',color='red',markersize=2,label=r'$x+4x^2$')
    ax.set_xlabel(r'x',fontsize=12)
    ax.set_ylabel(r'y',fontsize=12)
    ax.set_xlim([0,100])
    ax.set_ylim([0,np.max(y)])
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend(loc='lower right',fontsize=12,markerscale=4)
    fig.tight_layout()
    fig.savefig('poly.png', format='png')

    # JCR time before interactive plot
    #     so it doesn't time plot being open 
    time2  = time.time()
    print("    plot [s]: {0:.3f}".format((time2-time1)))

    plt.show()
    plt.close(fig)

    print("c to continue")
    pdb.set_trace()
