# M. Rizzo Smith
# v1 1/25/24
# Week 3 Lecture 2 Tasks
# ASTR 8080 coordinate system transformations


# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
import astropy 
from matplotlib import rc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import sfdmap
import extinction

# FUNCTIONS
###############################
###############################

# None for now

# MAIN
###############################
###############################
if __name__ == '__main__':
    # MRS Initial setup of dust at one position
    # MRS Setup SkyCoord object of given RA and Dec
    ra, dec = '00h42m30s', '+41d12m00s'
    c = SkyCoord(ra, dec, unit=(u.hourangle, u.degree)).galactic

    # MRS Import the dust map data from directory holding sfddata
    dust_dir = './sfddata'
    m = sfdmap.SFDMap(dust_dir, scaling=1)

    # MRS Obtain reddening at ra and dec position
    ebv = m.ebv(c.l.value,c.b.value, frame='galactic')
   
    # MRS Obtain reddening with just RA and Dec
    c = SkyCoord(ra, dec)
    ebv = m.ebv(c.ra.value,c.dec.value)

    # MRS Extinction magnatiudes for the ugriz filters
    wave = np.array([3543., 4770., 6231, 7625., 9134.])
    A = extinction.fitzpatrick99(wave, 3.1*ebv)

    # MRS Initialize Sky Coord objects of the two test Quasars
    obj1 = SkyCoord('246.933 40.795', unit=(u.degree, u.degree))
    obj2 = SkyCoord('236.562 2.440', unit=(u.degree, u.degree))

    # MRS Store the ugri magnitudes of quasar 1 and 2 
    mag1 = np.array([18.82, 18.81, 18.74, 18.81])
    mag2 = np.array([19.37, 19.10, 18.79, 18.72])

    # MRS Plot the g-r vs the r-i colors of each quasar
    plt.plot(mag1[1]-mag1[2], mag1[2]-mag1[3], 'o', color = 'red', label='Quasar 1')
    plt.plot(mag2[1]-mag2[2], mag2[2]-mag2[3], 'o',  color = 'orange', label='Quasar 2')
    plt.xlabel('g-r')
    plt.ylabel('r-i')
    plt.legend()
    plt.title('g-r vs r-i For Two Quasars') 
    plt.show()
    
    # MRS Obtain the extinction mags for each quasar
    ebv1 = m.ebv(obj1.ra.value, obj1.dec.value) 
    A1 = extinction.fitzpatrick99(wave, 3.1*ebv1)
    
    # MRS correct for dust for quasar 1
    cor_mag1 = mag1 - A1[0:4]
    
    # MRS Obtain and correct for dust for quasar 2
    ebv2 = m.ebv(obj2.ra.value, obj2.dec.value)
    A2 = extinction.fitzpatrick99(wave, 3.1*ebv2)
    cor_mag2 = mag2 - A2[0:4]

    # MRS Plot the newly corrected colors of both quasars
    plt.plot(mag1[1]-mag1[2], mag1[2]-mag1[3], 'o', color = 'red', label='Quasar 1')
    plt.plot(mag2[1]-mag2[2], mag2[2]-mag2[3], 'o',  color = 'orange', label='Quasar 2')
    plt.plot(cor_mag1[1]-cor_mag1[2], cor_mag1[2]-cor_mag1[3], 's', color = 'red', label = 'Extinction Corrected Quasar 1')
    plt.plot(cor_mag2[1]-cor_mag2[2], cor_mag2[2]-cor_mag2[3], 's',  color = 'orange', label='Extinction Corrected Quasar 2')
    plt.xlabel('g-r')
    plt.ylabel('r-i')
    plt.title('g-r vs r-i For Two Quasars')
    plt.legend()
    plt.show()
    
    # MRS Create a 100x 100 grid spcaed by 1 degree in ra and dec centered at 236.6 2.4
    ra_range1 = np.arange(236.6-49, 236.6 + 50, 1)
    dec_range1 = np.arange(2.4-50, 2.4 + 50, 1)
    ra_mesh1, dec_mesh1 = np.meshgrid(ra_range1, dec_range1)
    
    # MRS Create a 100x100 grid spaced by 1.3 degrees in ra and 1 degree in dec centered at 246.9 40.8
    ra_range2 = np.arange(246.9 -65, 246.9+65, 1.3)
    dec_range2 = np.arange(40.8 - 50, 40.8+50, 1)
    ra_mesh2, dec_mesh2 = np.meshgrid(ra_range2, dec_range2)
    
    # MRS Create SkyCoord objects for each of the points in the two grids
    map1 = SkyCoord(ra_mesh1, dec_mesh1, unit=(u.degree, u.degree)).galactic
    map2 = SkyCoord(ra_mesh2, dec_mesh2, unit=(u.degree, u.degree)).galactic
    
    # MRS Obtain the reddening at each position of each map
    ebv_map1 = m.ebv(map1.l.value, map1.b.value, frame='galactic')
    ebv_map2 = m.ebv(map2.l.value, map2.b.value, frame='galactic')

    # MRS Initializing the coordinates of the Galactic Plane
    GP_l = np.arange(0,360, 1)
    GP_b = np.zeros(360)

    # MRS Transform from Galactic l,b to ra, dec
    Gal_plane = SkyCoord(GP_l, GP_b, frame='galactic', unit=(u.degree, u.degree))
    Gal_plane = Gal_plane.transform_to('icrs')
    
    # MRS Plot the Contour of extinction for each grid point centered on quasar 1, galactic plane, and both quasar positions
    ebv_max1 = np.max(ebv_map1)
    clevels1 = [0.999, 0.99, 0.9]
    clevels1 = (1-np.array(clevels1))*ebv_max1
    plt.figure()
    cs = plt.contour(ra_mesh1, dec_mesh1, ebv_map1, levels=clevels1, colors=['k', 'k', 'k'])
    plt.clabel(cs)
    plt.title('Dust Map Centered on Quasar 1')
    plt.xlim(np.min(ra_range1)-1,np.max(ra_range1)+1 )
    plt.ylim(np.min(dec_range1)-1, np.max(dec_range1)+1)
    plt.plot(Gal_plane.ra.degree, Gal_plane.dec.degree, '-', color='red')
    plt.plot(obj1.ra.degree, obj1.dec.degree, 'o', color = 'red', label='Quasar 1' )
    plt.plot(obj2.ra.degree, obj2.dec.degree, 'o', color='orange', label='Quasar 2')
    plt.legend()
    plt.show()
   
    # MRS Plot the contour of extinction for each grid point centered on quasar 2, galactic plane, and both quasar positions
    plt.figure()
    cs = plt.contour(ra_mesh2, dec_mesh2, ebv_map2, levels=1000)
    plt.clabel(cs)
    plt.title('Dust Map Centered on Quasar 2')
    plt.xlim(np.min(ra_range2)-1,np.max(ra_range2)+1 )
    plt.ylim(np.min(dec_range2)-1, np.max(dec_range2)+1)
    plt.plot(Gal_plane.ra.degree, Gal_plane.dec.degree, '-', color='red')
    plt.plot(obj1.ra.degree, obj1.dec.degree, 'o', color = 'red', label='Quasar 1' )
    plt.plot(obj2.ra.degree, obj2.dec.degree, 'o', color='orange', label='Quasar 2')
    plt.legend()
    plt.show()
   


