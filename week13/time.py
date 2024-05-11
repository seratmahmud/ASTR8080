# Serat
# Classwork, week13
# Date: 4/2/2024


# Import
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

df = pd.read_csv('result.csv')

plt.figure(figsize=(10,6))
plt.scatter(df['ra'], df['dec'])
plt.xlabel('Right Ascension (Degrees)')
plt.ylabel('Decliantion (Degrees)')
plt.show()
plt.savefig('ra_vs_dec.png')

a = np.where((df['run'] != 206) & (df['run'] != 106))[0]
ra_c = df['ra'][a]
dec_c = df['dec'][a]

plt.figure(figsize=(10,6))
plt.scatter(ra_c, dec_c)
plt.xlabel('Right Ascension (Degrees)')
plt.ylabel('Decliantion (Degrees)')
plt.show()
plt.savefig('ra_vs_dec_c.png')


plt.figure(figsize=(10,6))
plt.scatter(ra_c, dec_c)
plt.xlabel('Right Ascension (Degrees)')
plt.ylabel('Decliantion (Degrees)')
plt.xlim(14.654, 14.658)
plt.ylim(0.72, 0.76)
plt.show()
plt.savefig('ra_vs_dec_zoomedin.png')

print(len(a))

ra_target = 14.654 + 0.0064
dec_target = 0.744
radius = 2/3600

cin = SkyCoord(ra_target*u.degree, dec_target*u.degree) 
csweeps = SkyCoord(ra_c*u.degree, dec_c*u.degree, unit = 'degree')
    
sep = cin.separation(csweeps)
m2 = np.where(sep < radius*u.degree)[0]

print(np.shape(m2))

plt.figure(figsize=(10,6))
plt.scatter(ra_c, dec_c, color='b')
plt.scatter(ra_target, dec_target, color='r')
plt.xlabel('Right Ascension (Degrees)')
plt.ylabel('Decliantion (Degrees)')
plt.xlim(14.654, 14.658)
plt.ylim(0.72, 0.76)
plt.show()
plt.savefig('ra_vs_dec_zoomedin.png')