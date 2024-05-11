# Serat
# Classwork, week15
# Date: 4/4/2024

from astropy.io import fits
import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from spec_sql_query import sdss_crossid_sql_query, gen_wget


def get_random_objects(filename, num_objects=10):
    with fits.open(filename) as hdul:
        data = hdul[1].data
        hdr = hdul[1].header
        print(hdr)
        indices = np.random.randint(0, len(data), size=num_objects)
        return data[indices]

qsos = get_random_objects('/home/saadsm/ASTR8080/runnoe/week11/qsos-ra180-dec30-rad3.fits')
stars = get_random_objects('/home/saadsm/ASTR8080/runnoe/week11/stars-ra180-dec30-rad3.fits')

ra = qsos['RA']  
dec = qsos['RA']

sql_query = 'SELECT s.ra, s.dec,s.mjd,s.plate,s.fiberid,s.run2d,p.psfMag_g\
    FROM #upload u\
        JOIN #x x ON x.up_id = u.up_id\
        JOIN SpecObjAll s ON s.specObjID = x.specObjID\
        JOIN PhotoTag p ON p.objID = s.bestObjID'


sql_data = sdss_crossid_sql_query(sql_query,ra,dec,FROM_FILE=True)
wget_coms,spec_files = gen_wget(sql_data['mjd'],sql_data['plate'],sql_data['fiberid'],sql_data['run2d'])

for com in wget_coms:
    os.system(com)
        
ra = stars['RA']  
dec = stars['RA']

sql_query = 'SELECT s.ra, s.dec,s.mjd,s.plate,s.fiberid,s.run2d,p.psfMag_g\
    FROM #upload u\
        JOIN #x x ON x.up_id = u.up_id\
        JOIN SpecObjAll s ON s.specObjID = x.specObjID\
        JOIN PhotoTag p ON p.objID = s.bestObjID'


sql_data = sdss_crossid_sql_query(sql_query,ra,dec,FROM_FILE=True)
wget_coms,spec_files = gen_wget(sql_data['mjd'],sql_data['plate'],sql_data['fiberid'],sql_data['run2d'])

for com in wget_coms:
    os.system(com)
    
for i in range(len(qsos)):
    print(f'A{i+1}', qsos['RA'][i], qsos['DEC'][i])
    
#cmd = "wget -i download_url.txt"
#returned_value = os.system(cmd)
#print('returned value:', returned_value)

filenames = glob.glob('*.fits')
for file in filenames:
    with fits.open(file) as hdul:
        data = hdul[1].data
        flux = data['flux']
        wl = 10**data['loglam']
        name = file.replace('.fits', '')
        plt.figure(figsize=(10, 6))
        plt.plot(wl, flux)
        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Flux (10^-17 ergs/s/cm2/Å)')
        plt.show()
        plt.savefig(f'plots/{name}.png')
    os.system(f'rm {file}')
    