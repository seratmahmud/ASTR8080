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
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import pdb
import time
import sys
sys.path.insert(0, '../../runnoe/week8/')
from sdss_sweep_data_index import sdss_sweep_data_index


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


# MAIN
###############################
###############################
if __name__ == '__main__':
    # checking magnitude to ugriz conversion
    PG1633A_sdssnav_ugriz = [17.28, 15.7, 15.18, 14.71, 14.55]
    PG1633A_WIYN_V_BV_UB_VR_RI = [15.256, 0.873, 0.320, 0.505, 0.511]
    PG1633A_g_conv = PG1633A_WIYN_V_BV_UB_VR_RI[0] + (0.60*PG1633A_WIYN_V_BV_UB_VR_RI[1]) - 0.12
    print(PG1633A_sdssnav_ugriz[1])
    print(PG1633A_g_conv)

    # finding star in SDSS sweeps
    ra_PG = '16h35m26s'
    dec_PG = '09d47m53s'
    PG_skycoord = SkyCoord(ra_PG, dec_PG, unit=[u.hourangle, u.degree])
    ras_PG = []
    decs_PG = []
    for i in range(2):
        ras_PG.append(PG_skycoord.ra.deg)
        decs_PG.append(PG_skycoord.dec.deg)
    sweeps = sdss_sweep_data_index(float(PG_skycoord.ra.value), float(PG_skycoord.dec.value), 0.5/60./60., sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final')
    sweeps = np.unique(sweeps)
    fitsweeps = [fits.open(fi) for fi in sweeps]
    objsweeps = [g[1].data for g in fitsweeps]
    objsweeps = np.hstack(objsweeps)
    objskycoord = SkyCoord(objsweeps['RA'], objsweeps['DEC'], unit=u.degree)
    objectseps = PG_skycoord.separation(objskycoord)
    objectind = np.where(objectseps.arcsecond < 2 )
    
    print(objskycoord[objectind]) 
    # finding ugriz for faint object
    faintobj_sc = SkyCoord(179.34594, -0.33714, unit=u.degree)
    faintobj_ugriz = [24.6, 22.12, 21.41, 21.14, 21.38]
    sweeps = sdss_sweep_data_index(float(faintobj_sc.ra.value), float(faintobj_sc.dec.value), 0.5/60./60., sweepdir='/astro/astr8020/dr15/eboss/sweeps/dr13_final')
    sweeps = np.unique(sweeps)
    fitsweeps = [fits.open(fi) for fi in sweeps]
    objsweeps = [g[1].data for g in fitsweeps]
    objsweeps = np.hstack(objsweeps)
    objskycoord = SkyCoord(objsweeps['RA'], objsweeps['DEC'], unit=u.degree)
    objectseps = faintobj_sc.separation(objskycoord)
    objectind = np.where(objectseps.arcsecond < 2 )
    fo_sweep_ugrizflux = objsweeps[objectind]['PSFFLUX']
    fo_sweep_ugrizmags = [(22.5-(2.5*np.log10(flux))) for flux in fo_sweep_ugrizflux] # convert nanomaggies to magnitudes
    print('navigator vals', faintobj_ugriz)
    print('sweep     vals', fo_sweep_ugrizmags)
    
    print('hello world XD')

