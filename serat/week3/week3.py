# S. Saad
# ASTR 8080 week3, classwork1

# to import this from another directory:
# import sys
# sys.path.insert(0, '../week1/')
# import polycalc
# from polycalc import get_poly_o3

# IMPORT BLOCK
###############################
###############################
from astropy.coordinates import AltAz
from astropy.coordinates import EarthLocation
import astropy.units as u
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, Galactic, AltAz, EarthLocation
from astropy.coordinates import get_body, HeliocentricTrueEcliptic



# FUNCTIONS
###############################
###############################
    

# Function to convert RA and Dec to Galactic coordinates
def ra_dec_to_galactic(ra, dec):
    ra_ngp = 192.859508  # RA of the North Galactic Pole
    dec_ngp = 27.128336  # Dec of the North Galactic Pole
    l_omega = 32.931918  # Galactic longitude of the ascending node of the Galactic plane

    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)
    ra_ngp_rad = np.radians(ra_ngp)
    dec_ngp_rad = np.radians(dec_ngp)

    b = np.arcsin(np.sin(dec_rad) * np.sin(dec_ngp_rad) + 
                  np.cos(dec_rad) * np.cos(dec_ngp_rad) * np.cos(ra_rad - ra_ngp_rad))
    
    l = np.arctan2(np.sin(dec_rad) - np.sin(b) * np.sin(dec_ngp_rad),
                   np.cos(dec_rad) * np.sin(ra_rad - ra_ngp_rad) * np.cos(dec_ngp_rad)) + np.radians(l_omega)
    
    l = np.mod(np.degrees(l), 360)
    b = np.degrees(b)

    return l, b


def main():
    
    # Part 1:
    ra = 10.68458  
    dec = 41.26917 
    distance = 1  

    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, distance=distance*u.dimensionless_unscaled, frame='icrs')
    cartesian_coords = c.represent_as('cartesian')
    x, y, z = cartesian_coords.x, cartesian_coords.y, cartesian_coords.z

    x_eq = np.cos(np.radians(dec)) * np.cos(np.radians(ra))
    y_eq = np.cos(np.radians(dec)) * np.sin(np.radians(ra))
    z_eq = np.sin(np.radians(dec))

    print("Checking if the coordinate transform works:")
    print(x, y, z)
    print(x_eq, y_eq, z_eq)
    
    # Part 2:
    # Coordinates for Galactic Center
    ra = '17h45m40.0409s'
    dec = '-29d00m28.118s'

    galactic_center = SkyCoord(ra, dec, frame='icrs')
    constellation = galactic_center.get_constellation()

    print("The Galactic Center is in the constellation:", constellation)
    
    
    # Part 3:
        
    nashville_lat = 36.0  # Latitude of Nashville in degrees
    ra_monthly = np.linspace(0, 24, 100) * 15  # Convert hours to degrees
    dec = nashville_lat
    galactic_coords = []
    
    for ra in ra_monthly:    
        galactic_coords.append(ra_dec_to_galactic(ra, dec))

    l_values, b_values = zip(*galactic_coords)

    plt.figure(figsize=(10, 6))
    plt.scatter(l_values, b_values, marker='o')
    plt.xlabel('Galactic Longitude (l) in degrees')
    plt.ylabel('Galactic Latitude (b) in degrees')
    plt.show()
    
    # Part 4:
    now = Time.now()
    planets = ['mercury', 'venus', 'mars']
    planet_positions = {planet: get_body(planet, now).transform_to(HeliocentricTrueEcliptic()) for planet in planets}

    longitudes = [planet_positions[planet].lon for planet in planets]
    latitudes = [planet_positions[planet].lat for planet in planets]

    plt.figure(figsize=(8, 6))
    for i, planet in enumerate(planets):
        plt.scatter(longitudes[i].degree, latitudes[i].degree, label=planet.capitalize())
        plt.text(longitudes[i].degree, latitudes[i].degree, f' {planet.capitalize()}')

    plt.title('Current Positions of Mercury, Venus, and Mars in Ecliptic Coordinates')
    plt.xlabel('Ecliptic Longitude (degrees)')
    plt.ylabel('Ecliptic Latitude (degrees)')
    plt.legend()
    plt.show()
     

# MAIN
###############################
###############################
if __name__ == '__main__':
    main()