# S. Saad
# v1 3/26/2024
# ASTR 8080 HW4


##############
#    Import
##############
import numpy as np
from astropy.coordinates import SkyCoord, search_around_sky
from astropy import units as u
import time
from astropy.io import fits
from sdss_sweep_data_index import sdss_sweep_data_index
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")


###############
#  Functions
###############

def main():
    ##############
    #  Problem 1
    ##############
    
    # SS Recording time
    epoch_0 = time.time()
    
    # SS Extracting data of the first sources
    path = "/astro/astr8020/FIRST/first_08jul16.fits"
    hdul = fits.open(path)
    data = hdul[1].data
    hdr = hdul[1].header

    # SS Getting the coordinates of the first sources
    ra_first = data['RA']
    dec_first = data['DEC']
 
    # SS Given coordinates of the designated areas and radius
    ra_one = 163
    ra_two = 167
    dec_one = 50
    dec_two = 50
    rad = 2

    # SS Getting the first sources in the designated area which will be in the survey
    c_one = SkyCoord(ra=ra_one*u.degree, dec=dec_one*u.degree) 
    c_two = SkyCoord(ra=ra_two*u.degree, dec=dec_two*u.degree) 
    c_survey = SkyCoord(ra=ra_first*u.degree, dec=dec_first*u.degree) 
    sep1 = c_one.separation(c_survey)
    m1 = np.where(sep1 < rad*u.degree)
    sep2 = c_two.separation(c_survey)
    m2 = np.where(sep2 < rad*u.degree)
    m = np.unique(np.append(m1, m2))

    # SS Coordinates of the first sources in the survey
    ra_first_in = ra_first[m]
    dec_first_in = dec_first[m]

    # SS Plotting the Right Ascension and Declination of the sources
    plt.figure(figsize=(10,6))
    plt.scatter(ra_first, dec_first, label="All FIRST Sources")
    plt.scatter(ra_first_in, dec_first_in, label="Sources in Survey")
    plt.xlabel("Right Ascension (Degrees)")
    plt.ylabel("Declination (Degrees)")
    plt.title("Coordinates plotting of the FIRST sources in the new survey")
    plt.legend()
    plt.show()
    plt.savefig("first_coord.png")
    
    # SS Recording and printing time
    epoch_1 = time.time()
    print("Total time taken for problem 1:", epoch_1-epoch_0, "s")

    ##############
    #  Problem 2
    ##############
    
    # SS Getting the files names of all the swfiles and wise files
    radius = 1/3600
    sweepdir="/astro/astr8020/dr15/eboss/sweeps/dr13_final"
    swfiles = sdss_sweep_data_index(ra_first_in, dec_first_in,\
        radius, sweepdir=sweepdir, all=all, verbose=False)
    #swfiles_stargal = [file.replace("star", "stargal-primary") for file in swfiles ]
    swfiles_wise_stargal = [file.replace("star", "wise-stargal-primary") for file in swfiles ]

    # SS Reading through SDSS swfiles
    objs = []
    for file in swfiles:
        try:
            hdul_objs = fits.open(file)
            objs.append(hdul_objs[1].data)
        except:
            objs.append(nan)
            continue
    objs = np.hstack(objs)

    # SS Checking for the sources that are within 1" of the sources in the survey
    csweeps = SkyCoord(ra=objs["RA"]*u.degree, dec=objs["DEC"]*u.degree)
    cin = SkyCoord(ra=ra_first_in*u.degree, dec=dec_first_in*u.degree)
    idx1, idx2, _, _ = search_around_sky(csweeps, cin, (radius)*u.deg)
    matched_objs = objs[idx1]
    
    # SS Getting the flux and magnitude of the sources
    psfflux =  np.array(matched_objs['PSFFLUX'])
    mag =  np.array(22.5 - 2.5 * np.log10(psfflux))
    #print("PSFFLUX Magnitude:", mag)

    # SS Getting the RA and DEC of these sources for problem 4
    ra_match = matched_objs["RA"]
    dec_match = matched_objs["DEC"]

    # SS Reading through wise files
    objs = []
    for file in swfiles_wise_stargal:
        try:
            hdul_objs = fits.open(file)
            objs.append(hdul_objs[1].data)
        except:
            objs.append(nan)
            continue

    # SS Getting the data for the wise files
    objs = np.hstack(objs)
    matched_objs = objs[idx1]

    # SS Getting the flux values from the wise data
    w1flux = np.array(matched_objs['W1_NANOMAGGIES'])
    w2flux = np.array(matched_objs['W2_NANOMAGGIES'])
    #mag_w1 =  np.array(22.5 - 2.5 * np.log10(w1flux))
    #mag_w2 =  np.array(22.5 - 2.5 * np.log10(w2flux))
    #print('W1 Magnitude (wise-stargal):', mag_w1)
    #print('W2 Magnitude (wise-stargal):', mag_w2)
    
    # SS Recording and printing time
    epoch_2 = time.time()
    print("Total time taken for problem 2:", epoch_2-epoch_1, "s")

    ##############
    #  Problem 3
    ##############
    
    # SS Getting the index of the brightest source
    ir1 = np.where(w1flux == np.max(w1flux))[0]

    # SS Getting the flux of the brightest sources and putting them in one array 
    w1flux = w1flux * 309.54e-9
    w2flux = w2flux * 171.787e-9
    wflux_ir1 = np.array([w1flux[ir1], w2flux[ir1]])
    psfflux_ir1 = psfflux[ir1] * 3631e-9
    flux_ir1 = np.append(psfflux_ir1, wflux_ir1)
    
    # SS The given wavelength for each of the fluxes
    wavelength_mag = np.array([3543, 4770, 6231, 7625, 9134, 34000, 46000])

    # SS Plotting the wavelength vs flux
    plt.figure(figsize=(10, 6))
    plt.scatter(wavelength_mag, flux_ir1)
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux (Jansky)')
    plt.title('FLux vs wavelngth for the brightest source')
    plt.show()
    plt.savefig("wl_vs_mag.png")
    
    # SS Recording and printing time
    epoch_3 = time.time()
    print("Total time taken for problem 3:", epoch_3-epoch_2, "s")
    
    ##############
    #  Problem 4
    ##############

    # SS Getting the coordinate of the brightest source
    ra_ir1 = ra_match[ir1]
    dec_ir1 = dec_match[ir1]

    # SS Printing the coordinates of the brightest source
    print(f"The coordinates of the brightest object are: RA={ra_ir1} deg and DEC={dec_ir1} deg)")

    # SS Recording and printing time
    epoch_4 = time.time()
    print("Total time taken for problem 4:", epoch_4-epoch_3, "s")
    print("Total time taken to run the whole program:", epoch_4-epoch_0, "s")


if __name__ == '__main__':
    # SS Calling the main function
    main()
    
#######
###   Comment for Problem 4:
###     The spectra from SDSS shows that it's a quasar. It's bright
###     and it has been identified as QSO in SDSS Explore section.
###     In question 3 plot we saw that its w1 and w2 flux are much 
###     greater than the PSFFLUX. The source being a cosmological 
###     (very far) but bright source (like a quasar) also justifies 
###     the fact that it has greater brightness in infrared but not
###     in optical range.