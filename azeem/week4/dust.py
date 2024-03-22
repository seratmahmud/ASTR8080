# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pdb
import time

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates import AltAz
from astropy.coordinates import Galactic, Longitude

if __name__ == '__main__':
    #AB get our coordinates and transform to galactic coords
    ra = '00h42m30s'
    dec = '+41d12m00s'
    c = SkyCoord(ra,dec).galactic

    #AB import SFDMap and obtain reddining
    import sfdmap
    dustdir = 'C:/ASTR8080/ASTR8080/azeem/week4/sfddata/'
    m = sfdmap.SFDMap(dustdir, scaling=1)
    ebv = m.ebv(c.l.value,c.b.value,frame='galactic')
    c = SkyCoord(ra,dec)
    ebv = m.ebv(c.ra.value,c.dec.value)

    #AB import extinction and find extinction for SDSS ugriz filters
    import extinction
    wave = np.array([3543.,4770.,6231.,7625.,9134.])
    A = extinction.fitzpatrick99(wave, 3.1*ebv)
    print(A)

    #tasks
    #1
    #AB g,r,and i magnitudes for quasar 1
    gmag1 = 18.81
    rmag1 = 18.74
    imag1 = 18.81

    #AB g,r,and i magnitudes for quasar 2
    gmag2 = 19.10
    rmag2 = 18.79
    imag2 = 18.72

    #AB find g-r and r-i colors for both quasars and plot
    gr1 = gmag1 - rmag1
    ri1 = rmag1 - imag1

    gr2 = gmag2 - rmag2
    ri2 = rmag2 - imag2

    plt.scatter(ri1, gr1, label='quasar 1')
    plt.scatter(ri2,gr2, label='quasar 2')
    plt.xlabel('r-i')
    plt.ylabel('g-r')
    plt.legend()
    plt.show()
    #AB: No, the quasars do not have the same colors, and they should considering quasars have the same spectral shape relatively

    #AB input coords of the two quasars
    q1 = SkyCoord(246.933, 40.795, frame='icrs', unit = 'deg')
    q2 = SkyCoord(236.562, 2.440, frame='icrs', unit = 'deg')


    #AB find reddening for both given their ra and dec, and ugriz filters
    ebv1 = m.ebv(q1.ra.value,q1.dec.value)
    ebv2 = m.ebv(q2.ra.value,q2.dec.value)


    A1 = extinction.fitzpatrick99(wave, 3.1*ebv1)
    A2 = extinction.fitzpatrick99(wave, 3.1*ebv2)


    #AB find corrected magnitudes and plot
    corrgmag1 = gmag1 - A1[1]
    corrrmag1 = rmag1 - A1[2]
    corrimag1 = imag1 - A1[3]

    corrgmag2 = gmag2 - A2[1]
    corrrmag2 = rmag2 - A2[2]
    corrimag2 = imag2 - A2[3]

    plt.scatter(corrrmag1-corrimag1, corrgmag1-corrrmag1, label='quasar 1')
    plt.scatter(corrrmag2-corrimag2, corrgmag2-corrrmag2, label ='quasar 2')
    plt.xlabel('r-i$_{corr}$')
    plt.ylabel('g-r$_{corr}$')
    plt.legend()
    plt.show()

    #2
    #AB input center coordinates around first and second object
    xcent1 = 236.6
    ycent1 = 2.4

    xcent2 = 246.9
    ycent2 = 40.8

    #AB steps needed for 1 deg and 1.3 degree bins
    step1 = 100
    step2 = 103
    #AB make array centered around coord center for first object and make meshgrid
    xv = np.arange(xcent1 - (step1/2), xcent1 + (step1/2), 1)
    yv = np.arange(ycent1 - (step1/2), ycent1 + (step1/2), 1)
    xmap, ymap = np.meshgrid(xv,yv)

    #AB repeat for second object
    xw = np.arange(xcent2 - (step2/2), xcent2 + (step2/2), 1)
    yw = np.arange(ycent2 - (step2/2), ycent2 + (step2/2), 1)
    xmap2, ymap2 = np.meshgrid(xw,yw)
    #AB plot
    plt.plot(xmap,ymap, marker=',')
    plt.show()
    plt.plot(xmap2,ymap2, marker=',')
    plt.show()
