# S. Saad
# v1 8/26/2020
# v2 1/24/2022
# ASTR 8080 example for HW0


# IMPORT BLOCK
###############################
###############################
from astropy.io import fits
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import polycalc

# FUNCTIONS
###############################
###############################

def main():
    # Opening the FITS file and extractong the data we have
    fx = fits.open("struc.fits")
    objs = fx[1].data
    ra = objs["RA"]
    dec = objs["DEC"]

    # Generating random numbers in set1, set2, and set3
    set1 = np.random.randint(low=0, high=100, size=100)
    set2 = np.random.randint(low=0, high=100, size=100)
    set3 = np.random.randint(low=0, high=100, size=100)

    # Combine into a single array of 3-element arrays
    randomnum = np.array(list(zip(set1, set2, set3)))

    # Create a recarray
    data = np.core.records.fromarrays([ra, dec, randomnum], 
                                      names='ra,dec,randomnum',
                                      formats=['f8', 'f8', '(3,)i4'])

    # Converting and saving in Astropy Table format
    t = Table(data)
    t.write('output.fits', format='fits', overwrite=True)

    # Final Plotting
    extnc = objs["EXTINCTION"][:,0]
    a = np.where(extnc > 0.22)[0]
    plt.plot(ra[a], dec[a], "bx")
    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.show()

    

# MAIN
###############################
###############################
if __name__ == '__main__':
    main()
    print(help(polycalc))
    print(polycalc.get_poly_o3.__doc__)
