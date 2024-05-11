# Serat Saad

# Classwork 2, week 9
# Date: 3/7/2024


# Serat
# Classwork 2, week9
# Date: 3/7/2024



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
    
    ######################
    # Part 1 & 2
    star_file = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-003615-4-star.fits.gz')
    star_file = star_file[1].data
    stargal_file = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-003615-4-stargal-primary.fits.gz')
    stargal_file = stargal_file[1].data
    wise_stargal_file = fits.open('/astro/astr8020/dr15/eboss/sweeps/dr13_final/301/calibObj-003615-4-wise-stargal-primary.fits.gz')
    wise_stargal_file = wise_stargal_file[1].data
    
    # print(stargal_file)
    # print(wise_stargal_file)
    print('Data Length for Stargal: ', len(stargal_file))
    print('Data Length for Wise-Stargal: ', len(wise_stargal_file))
    

    ra = 143.209 # in deg
    dec = 36.701 # in deg
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
    mag =  22.5 - 2.5 * np.log10(s_obj['PSFFLUX'][0])
    print('PSFFLUX magnitude (star)', mag)
    
    
    
    
    # Part 3
    sweepdir="/astro/astr8020/dr15/eboss/sweeps/dr13_final"
    swfiles = sdss_sweep_data_index(ra, dec, radius, objtype='star',
                                    sweepdir=sweepdir, all=all, verbose=False)
    swfiles_stargal = [file.replace("star", "stargal-primary") for file in swfiles ]
    swfiles_wise_stargal = [file.replace("star", "wise-stargal-primary") for file in swfiles ]
    
    # Stargal
    objs_struct = []
    for file in swfiles_stargal:
        try:
            objs_struct.append(fits.open(file))
        except:
            continue
         
    objs        = [obj[1].data for obj in objs_struct ] 
    
    objs = np.hstack(objs)

    cin = SkyCoord(ra=ra*u.degree, dec=dec*u.degree) 
    csweeps = SkyCoord(ra=objs["RA"]*u.degree, dec=objs["DEC"]*u.degree)
    
    sep = cin.separation(csweeps)
    m2 = np.where(sep < radius*u.degree)[0]

    s_obj = objs[m2]
    mag =  22.5 - 2.5 * np.log10(s_obj['PSFFLUX'][0])
    print('PSFFLUX magnitude (stargal):', mag)
    
    # Wise-Stargal
    objs_struct = []
    for file in swfiles_wise_stargal:
        try:
            objs_struct.append(fits.open(file))
        except:
            continue
         
    objs        = [obj[1].data for obj in objs_struct ] 
    
    objs = np.hstack(objs)

    cin = SkyCoord(ra=ra*u.degree, dec=dec*u.degree) 
    csweeps = SkyCoord(ra=objs["RA"]*u.degree, dec=objs["DEC"]*u.degree)
    
    sep = cin.separation(csweeps)
    m2 = np.where(sep < radius*u.degree)[0]

    s_obj = objs[m2]
    mag_w1 =  22.5 - 2.5 * np.log10(s_obj['W1_NANOMAGGIES'])
    mag_w2 =  22.5 - 2.5 * np.log10(s_obj['W2_NANOMAGGIES'])
    print('W1 Magnitude (wise-stargal):', mag_w1)
    print('W2 Magnitude (wise-stargal):', mag_w2)
    

    
    



if __name__ == '__main__':
    main()






