# S. Saad
# ASTR 8080 week5, classwork1

# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp



# FUNCTIONS
###############################
###############################



def main():
    
    # Part 1
    ra = 360. * np.random.random(1000000)
    dec = (180/np.pi) * np.arcsin(1. - np.random.random(1000000) * 2.)
    
    #plt.figure(figsize=(10, 5))
    #plt.scatter(ra, dec, s=1, alpha=0.5)
    #plt.xlabel('RA (deg)')
    #plt.ylabel('Dec (degrees)')
    #plt.show()
    
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    
    theta = np.pi/2 - dec_rad
    phi = ra_rad
    
    nside = 1
    pixels = hp.ang2pix(nside, theta, phi)
    
    #pixel_area = hp.nside2pixarea(nside, degrees=True)

    #print(pixel_area)

    pixel_counts, _ = np.histogram(pixels, bins=np.arange(-0.5, hp.nside2npix(nside)+0.5, 1))
    

    w2 = np.where(pixels == 2)[0]
    w5 = np.where(pixels == 5)[0]
    w8 = np.where(pixels == 8)[0]


    plt.figure(figsize=(10, 5))
    plt.scatter(ra, dec, s=1, alpha=0.5, color='k', marker='.', label='All Points')

    plt.scatter(ra[w2], dec[w2], s=1, color='red', marker='.', label='Pixel 2')
    plt.scatter(ra[w5], dec[w5], s=1, color='blue', marker='.', label='Pixel 5')
    plt.scatter(ra[w8], dec[w8], s=1, color='green', marker='.', label='Pixel 8')

    plt.xlabel('RA (degrees)')
    plt.ylabel('Dec (degrees)')
    plt.legend()
    plt.savefig("Equatorial_pixel_plot.png")
    plt.show()

    #nside1_pixel_number = 5
    #nside2_pixels_start = nside1_pixel_number * 4
    #nside2_pixels_end = nside2_pixels_start + 4 
    #nside2_pixels_within_nside1_pixel5 = list(range(nside2_pixels_start, nside2_pixels_end))

    #print(nside2_pixels_within_nside1_pixel5)

    #nside = 2
    #pixels_nside2 = hp.ang2pix(nside, theta, phi)


    #indices = np.hstack([np.where(pixels_nside2 == x)[0] for x in nside2_pixels_within_nside1_pixel5])


    #plt.scatter(ra[indices], dec[indices], s=1, color='red', marker='.')

    #plt.xlabel('RA')
    #plt.ylabel('Dec')
    #plt.legend()
    #plt.show()

    

    
    
    

    
# MAIN
###############################
###############################
if __name__ == '__main__':
    main()



