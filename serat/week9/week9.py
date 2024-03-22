# Serat
# Classwork, wek9
# Date: 3/5/2024



#Import
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from time import time
from astropy.io import fits
from sdss_sweep_data_index import sdss_sweep_data_index


#Functions
############
##########

def main():
    ################
    # Part 1

    # UBVRI data
    u_b = 0.320
    b_v = 0.873
    r_i = 0.511
    v_r = 0.505
    v = 15.256

    u_g    =    1.28*(u_b)   + 1.13     
    g_r    =    1.02*(b_v)   - 0.22  
    r_i    =    0.91*(r_i) - 0.20   
    r_z    =    1.72*(r_i) - 0.41    
    g      =    v + 0.60*(b_v) - 0.12    
    r      =    v - 0.42*(b_v) + 0.11

    print('g magnitude from SDSS DR16 navigation tool:', 15.70)
    print('g magnitude through calculation:', g) 

    ######################
    # Part 2

    ra = 248.8583333 # in deg
    dec = 9.798055556 # in deg
    radius = 2/3600


    sweepdir="/astro/astr8020/dr15/eboss/sweeps/dr13_final"
    swfiles = sdss_sweep_data_index(ra, dec, radius, objtype='star',
                                    sweepdir=sweepdir, all=all, verbose=False)
    objs_struct = [ fits.open(file) for file in swfiles ]
    objs        = [obj[1].data for obj in objs_struct ] 
    objs = np.hstack(objs)

    cin = SkyCoord(ra=ra*u.degree, dec=dec*u.degree) 
    csweeps = SkyCoord(ra=objs["RA"]*u.degree, dec=objs["DEC"]*u.degree)
    
    sep = cin.separation(csweeps)
    m2 = np.where(sep < radius*u.degree)[0]

    s_obj = objs[m2]
    g_flux = s_obj['PSFFLUX'][0][1]
    g_mag =  22.5 - 2.5 * np.log10(g_flux)
    print('g magnitude from SDSS Sweep Data Index:', g_mag)
    

    ##################
    # Part 3

    ra_faint = 248.88855
    dec_faint = 9.80169

    # Navigate tool values
    u_nav = 24.30
    g_nav = 22.93
    r_nav = 21.79
    i_nav = 21.52
    z_nav = 21.13
    mag_nav = [u_nav, g_nav, r_nav, i_nav, z_nav]

    # Calculating the sweeps value for this faint obj
    swfiles_faint = sdss_sweep_data_index(ra_faint, dec_faint, radius, objtype='star',
                                    sweepdir=sweepdir, all=all, verbose=False)
    objs_struct_faint = [ fits.open(file) for file in swfiles_faint]
    objs_faint        = [obj[1].data for obj in objs_struct_faint] 
    objs_faint = np.hstack(objs_faint)

    cin_faint = SkyCoord(ra=ra_faint*u.degree, dec=dec_faint*u.degree) 
    csweeps_faint = SkyCoord(ra=objs["RA"]*u.degree, dec=objs["DEC"]*u.degree)
    
    sep_faint = cin_faint.separation(csweeps_faint)
    m2_faint = np.where(sep_faint < radius*u.degree)[0]

    s_obj_faint = objs_faint[m2_faint]
    flux_faint = s_obj_faint['PSFFLUX'][0]
    mag_faint =  22.5 - 2.5 * np.log10(flux_faint)

    print('UGRIZ mag from SDSS Nav Tool:', mag_nav)
    print('UGRIZ mag calculated using sweeps file:', mag_faint)
    



if __name__ == '__main__':
    main()





