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


#Functions
############
##########
def column(matrix, i):
    return [row[i] for row in matrix]


def main():
    ######################
    # Part 1

    #ra = 336.4388 # in deg
    #dec = -0.8343 # in deg
    #radius = 2/3600

    #sweepdir="/astro/astr8020/dr15/eboss/sweeps/dr13_final"
    #swfiles = sdss_sweep_data_index(ra, dec, radius, objtype='gal',
    #                                sweepdir=sweepdir, all=all, verbose=False)
    #objs_struct = [ fits.open(file) for file in swfiles ]
    #objs        = [obj[1].data for obj in objs_struct ] 
    #objs = np.hstack(objs)

    #cin = SkyCoord(ra=ra*u.degree, dec=dec*u.degree) 
    #csweeps = SkyCoord(ra=objs["RA"]*u.degree, dec=objs["DEC"]*u.degree)
    
    #sep = cin.separation(csweeps)
    #m2 = np.argmin(sep < radius*u.degree)
    #s_obj = objs[m2]
    #print(s_obj['RA'], s_obj['DEC'])
   


    #part 3
    ra = 180 # in deg
    dec = 30 # in deg
    radius = 3


    sweepdir="/astro/astr8020/dr15/eboss/sweeps/dr13_final"
    swfiles = sdss_sweep_data_index(ra, dec, radius, objtype='star', sweepdir=sweepdir, all=all, verbose=False)
    objs_struct = [ fits.open(file) for file in swfiles ]
    objs        = [obj[1].data for obj in objs_struct ]
    objs = np.hstack(objs)

    cin = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    csweeps = SkyCoord(ra=objs["RA"]*u.degree, dec=objs["DEC"]*u.degree)

    sep = cin.separation(csweeps)
    m2 = np.where(sep < radius*u.degree)[0]
    objs = objs[m2]    
    print(len(objs))
    flux = np.array(objs['PSFFLUX'])
    i_flux = column(flux, 3)
    i_mag = 22.5 - 2.5 * np.log10(i_flux)
    a = np.where(i_mag < 20)[0]
    objs = objs[a]
    print(len(objs))    

    ra_obj = objs['RA']
    dec_obj = objs['DEC']
    
    hdul = fits.open("/home/saadsm/ASTR8080/runnoe/week11/qsos-ra180-dec30-rad3.fits")
    data = hdul[1].data
    ra_cat = data["RA"]
    dec_cat = data["DEC"]

    c = SkyCoord(ra=ra_obj*u.degree, dec=dec_obj*u.degree)
    catalog = SkyCoord(ra=ra_cat*u.degree, dec=dec_cat*u.degree)
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)    
    quas = np.unique(objs[idx])
    print(len(quas))

    # part 4
    flag = 2**18
    w = np.where((objs['RESOLVE_STATUS'] & flag) != 0)[0]
    print(len(w))

if __name__ == '__main__':
    main()





