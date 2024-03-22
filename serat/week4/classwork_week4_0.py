# S. Saad
# ASTR 8080 week3, classwork0

# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import sfdmap
from astropy.io import fits
from astropy.wcs import WCS
from extinction import fitzpatrick99


# FUNCTIONS
###############################
###############################


def main():

    # Part 1
    #The following data is achieved from SDSS Nav tool

    r_o1=18.74
    g_o1=18.81
    i_o1=18.81

    r_o2=18.79
    g_o2=19.10
    i_o2=18.72
    
    plt.scatter(g_o1-r_o1, r_o1-i_o1)
    plt.scatter(g_o2-r_o2, r_o2-i_o2)
    plt.show()

    # They don't have the same color, let's correct them and see what happens

    dust_map = fits.open('sfddata-master/SFD_dust_4096_ngp.fits')
    ebv_data = dust_map[0].data

    coords_o1 = SkyCoord('246.933d', '40.795d', frame='icrs')
    wcs_info_o1 = WCS(dust_map[0].header)
    pixel_coords_o1 = wcs_info_o1.world_to_pixel(coords_o1)
    ebv_value_o1 = ebv_data[int(pixel_coords_o1[1]), int(pixel_coords_o1[0])]
    R_V = 3.1
    A_V_o1 = R_V * ebv_value_o1

    coords_o2 = SkyCoord('236.562d', '2.440d', frame='icrs')
    wcs_info_o2 = WCS(dust_map[0].header)
    pixel_coords_o2 = wcs_info_o2.world_to_pixel(coords_o2)
    ebv_value_o2 = ebv_data[int(pixel_coords_o2[1]), int(pixel_coords_o2[0])]
    R_V = 3.1
    A_V_o2 = R_V * ebv_value_o2

    # Getting the extinction for each of the bands u, g, r, i, z
    wave = np.array([3543., 4770., 6231, 7625., 9134.])
    A_o1 = fitzpatrick99(wave, A_V_o1)
    A_o2 = fitzpatrick99(wave, A_V_o2)

    r_o1_c=r_o1-A_o1[2]
    g_o1_c=g_o1-A_o1[1]
    i_o1_c=i_o1-A_o1[3]

    r_o2_c=r_o2-A_o2[2]
    g_o2_c=g_o2-A_o2[1]
    i_o2_c=i_o2-A_o2[3]

    plt.scatter(g_o1-r_o1, r_o1-i_o1, color='blue', marker='^')
    plt.scatter(g_o2-r_o2, r_o2-i_o2, color='blue', marker='^')
    plt.scatter(g_o1_c-r_o1_c, r_o1_c-i_o1_c, color='red')
    plt.scatter(g_o2_c-r_o2_c, r_o2_c-i_o2_c, color='red')
    plt.show()

    # Part 2

    ra_center_1, dec_center_1 = 236.6, 2.4
    ra_bins_1, dec_bins_1 = 1, 1
    size = 100
    ra_range_1 = np.linspace(ra_center_1 - 50 * ra_bins_1, ra_center_1 + 49 * ra_bins_1, size)
    dec_range_1 = np.linspace(dec_center_1 - 50 * dec_bins_1, dec_center_1 + 49 * dec_bins_1, size)
    ra_grid_1, dec_grid_1 = np.meshgrid(ra_range_1, dec_range_1)
    sc_1 = SkyCoord(ra_grid_1, dec_grid_1, unit='deg')
    sc_1.galactic
    dustdir = 'sfddata-master/'
    m = sfdmap.SFDMap(dustdir, scaling=1)
    ebv_o1 = m.ebv(sc_1.ra.value,sc_1.dec.value, frame='galactic')

    ra_center_2, dec_center_2 = 246.9, 40.8
    ra_bins_2, dec_bins_2 = 1.3, 1  
    ra_range_2 = np.linspace(ra_center_2 - 50 * ra_bins_2, ra_center_2 + 49 * ra_bins_2, size)
    dec_range_2 = np.linspace(dec_center_2 - 50 * dec_bins_2, dec_center_2 + 49 * dec_bins_2, size)
    ra_grid_2, dec_grid_2 = np.meshgrid(ra_range_2, dec_range_2)
    sc_2 = SkyCoord(ra_grid_2, dec_grid_2, unit='deg')
    sc_2.galactic
    dustdir = 'sfddata-master/'
    m = sfdmap.SFDMap(dustdir, scaling=1)
    ebv_o2 = m.ebv(sc_2.ra.value,sc_2.dec.value, frame='galactic')
    
    ebv_max_1 = np.max(ebv_o1)
    ebv_max_2 = np.max(ebv_o2)
    
    print(ebv_max_1, ebv_max_2)
    
    clevels = np.array([0.999, 0.99, 0.9])
    clevels_o1 = (1-clevels)*ebv_max_1
    clevels_o2 = (1-clevels)*ebv_max_2
    
    
    plt.contour(sc_1.ra.value, sc_1.dec.value, ebv_o1, levels=clevels_o1)
    plt.xlabel("Galactic Longitude (l)")
    plt.ylabel("Galactic Latitude (b)")
    plt.show()
    
    plt.contour(sc_2.ra.value, sc_2.dec.value, ebv_o2, levels=clevels_o2)
    plt.xlabel("Galactic Longitude (l)")
    plt.ylabel("Galactic Latitude (b)")
    plt.ylim(-20, 20)
    plt.xlim(265, 310)
    plt.show()
    
    
# MAIN
###############################
###############################
if __name__ == '__main__':
    main()