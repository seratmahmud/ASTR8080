# Serat
# Classwork, week11
# Date: 3/19/2024



#Import
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from time import time
from astropy.io import fits
from sdss_sweep_data_index import sdss_sweep_data_index
import matplotlib.pyplot as plt

#Functions
############
##########

def main():
    ######################
    # Part 1

    ra = 180 # in deg
    dec = 30 # in deg
    radius = 3


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

    ra_obj = s_obj["RA"]
    dec_obj = s_obj["DEC"]
   
    hdul = fits.open("/home/saadsm/ASTR8080/runnoe/week11/qsos-ra180-dec30-rad3.fits")
    data = hdul[1].data
    ra_cat = data["RA"]
    dec_cat = data["DEC"]

    c = SkyCoord(ra=ra_obj*u.degree, dec=dec_obj*u.degree)
    catalog = SkyCoord(ra=ra_cat*u.degree, dec=dec_cat*u.degree)
    idx, d2d, d3d = catalog.match_to_catalog_sky(c)
    quas = objs[idx]
    print(len(quas))
    
    #Calculating the flux and mag 
    flux = objs['PSFFLUX']
    mag = np.array(22.5 - 2.5 * np.log10(flux))
    
    #Extinction Correction
    extnc = np.array(objs['EXTINCTION'])
    correct_mag = mag + extnc

    #Calculating color cuts
    u_g = [row[0] - row[1] for row in correct_mag]
    g_r = [row[1] - row[2] for row in correct_mag]
    r_i = [row[2] - row[3] for row in correct_mag]
    r_z = [row[3] - row[4] for row in correct_mag]
    
    #Calculating the flux and mag
    quas_flux = quas['PSFFLUX']
    quas_mag = np.array(22.5 - 2.5 * np.log10(quas_flux))

    #Extinction Correction
    quas_extnc = np.array(quas['EXTINCTION'])
    quas_correct_mag = quas_mag + quas_extnc

    #Calculating color cuts
    quas_u_g = [row[0] - row[1] for row in quas_correct_mag]
    quas_g_r = [row[1] - row[2] for row in quas_correct_mag]
    quas_r_i = [row[2] - row[3] for row in quas_correct_mag]
    quas_r_z = [row[3] - row[4] for row in quas_correct_mag]   
    
    
    plt.figure(figsize=(10,6))
    plt.scatter(g_r, u_g, c='black')
    plt.scatter(quas_g_r, quas_u_g, c='red')
    plt.xlabel('g-r')
    plt.ylabel('u-g')
    plt.show()
    plt.savefig('u-g_vs_g-r.png')











     

if __name__ == '__main__':
    main()





