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
from astropy.coordinates import Galactic, Longitude, get_body

def planets_ecliptic(planets,dates):
    """
    This is a function that plots the positions of any specified object
    given its position and time of observation in the ecliptic coordinate system.

    Args:
        planets (list): An input list that specifies the object name as a string
        dates (list): Input array list of observation times in UTC in the format
        MM-DD-YYYY.
    Returns:
        N/A, just plots
    """

    #-------------------------------------------------------------
    #+
    # PURPOSE:
    #   Plot planet positions in ecliptic coordiantes.
    #
    #
    #
    # CALLING SEQUENCE:
    #   planets_ecliptic(planets,dates)
    #
    # INPUTS:
    #   planets - ['mercury', 'venus', ...]
    #   dates - ['Observation Dates']
    #-
    #-------------------------------------------------------------

    #AB input Earth location at 0,0,0 and intialize Time object
    loc = EarthLocation(0, 0, 0)
    time = Time(dates)

    #AB loop over each planet, in this case we have 5
    for i in np.arange(0,5):
        planet = get_body(planets[i],time,loc) #SkyCoord object of each planet in list
        ecliptic = planet.transform_to('heliocentrictrueecliptic') #coordinate transformation
        plt.scatter(ecliptic.lon.deg,ecliptic.lat.deg, label=planets[i])
    #AB plot
    plt.legend()
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('First 5 planet positions from 01-01-2020 to 01-01-2030 at 7:00 AM and PM MST in ecliptic coordinates')
    plt.show()
    return


def planets_equatorial(planets,dates):
    """
    This is a function that plots the positions of any specified object
    given its position and time of observation in the equatorial coordinate system.
    It also finds the lowest airmass for each object and specifies position, and therefore
    best observing time, on the plot.

    Args:
        planets (list): An input list that specifies the object name as a string
        dates (list): Input array list of observation times in UTC in the format
        MM-DD-YYYY.
    Returns:
        N/A, just plots
    """

    #-------------------------------------------------------------
    #+
    # PURPOSE:
    #   Plot planet positions in equatorial coordiantes and
    #   find the lowest airmass for each.
    #
    #
    # CALLING SEQUENCE:
    #   planets_equatorial(planets,dates)
    #
    # INPUTS:
    #   planets - ['mercury', 'venus', ...]
    #   dates - ['Observation Dates']
    #-
    #-------------------------------------------------------------

    #AB make empty array to store lowest airmass
    best_airmass = np.empty((len(planets),3), dtype=object)
    #AB give location of observing site (APO in this case)
    loc1 = EarthLocation.of_site('APO')
    time = Time(dates, scale='utc') #time object
    #AB loop over each planet, in this case 5
    for i in np.arange(0,5):
        c1 = get_body(planets1[i], time,loc1) #SkyCoord object of planets
        #AB find altitude and azimuth given observing site and location
        altaz = AltAz(location=loc1, obstime=time)
        planetaa = c1.transform_to(altaz)
        #AB find airmasses
        airmass = planetaa.secz
        plt.scatter(c1.ra.deg,c1.dec.deg, label=planets[i]) #initialize plot
        #AB since we want lowest airmass, but don't want negatives, make sure all values are positive
        pos = np.where(airmass>0)
        best_air = airmass[pos]

        # AB find lowest positive value for each airmass in loop
        if best_air.size >0:
            min_index = np.argmin(best_air) #Stores index
            best_airmass[i,0] = best_air[min_index].value #Stores minimum positive airmass
            best_airmass[i,1] = np.where(airmass ==  best_air[min_index])[0][0] #Keep index of planet positions
            best_airmass[i,2] = dates[best_airmass[i,1]] #Keeps record of planet observation date

            # AB Print out planet name, airmass, and date of best airmass
            print(f'At APO, {planets[i]} Lowest Airmass:', best_airmass[i,0], 'On:', best_airmass[i,2], 'UTC')
            #AB Need ra and dec for best airmass of each planet
            low_ra = c1[best_airmass[i, 1]].ra.degree
            low_dec = c1[best_airmass[i, 1]].dec.degree
            plt.plot(low_ra, low_dec, '*', markeredgecolor='black', markersize=10, label='lowest airmass')
    #plot
    plt.legend()
    plt.xlabel('RA')
    plt.ylabel('dec')
    plt.title('First 5 planet positions from 01-01-2020 to 01-01-2030 at 7:00 AM and PM MST in equatorial coordinates at APO')
    plt.show()
    return

def vis_qso(month, year):
    """
    This is a function that finds the best quasar to observe (i.e. lowest airmass) given a specified month
    and year. It will also give the quasar's position in ra and dec, the date of observation,
    and airmass. The quasar is retrieved from a list of input objects.

    Args:
        month (int): An integer for the month of observation ranging from 1-12.
        year (int): An integer for the year of observation in the YYYY format.
    Returns:
        N/A, just prints
    """

    #-------------------------------------------------------------
    #+
    # PURPOSE:
    #   Find lowest airmass quasar given a list.
    #   Prints best quasar's position, date of observation, and airmass.
    #
    #
    # CALLING SEQUENCE:
    #   vis_qso(month, year)
    #
    # INPUTS:
    #   month - 1
    #   year - 2024
    #-
    #-------------------------------------------------------------

    #AB makes empty array to store quasar positions
    all_pos = []
    quasars = 'HW1quasarfile.dat' #quasar file
     #AB To open file and make an array:
    with open(quasars, 'r') as file:
        for line in file:
            qpos = line
            qpos = f'{qpos[:2]} {qpos[2:4]} {qpos[4:9]} {qpos[9:12]} {qpos[12:14]} {qpos[14:]}'
            cquasar = SkyCoord(qpos, unit=(u.hourangle, u.deg)) #converts to SkyCoord object
            all_pos.append(cquasar) #make array
    all_pos = np.array(all_pos)

    #AB Need number of days in each month, accounting for leap year
    days = [31, 29 if year % 4 == 0 else 28, 31, 30,31, 30, 31, 31, 30, 31, 30, 31]

    # AB Observing at 11PM MST, need to convert to UTC
    timeUTC = Time(f"{year}-{month:02d}-01 11:00:00") + 7 * u.hour

    #AB store for later
    best_airmass = np.inf
    best_quasar = None
    best_time = None

    loc1 = EarthLocation.of_site('APO') #site of observing

    #AB loop to get quasar position and best airmass
    for offset in range(days[month-1]):
        time = timeUTC + offset * u.day #needed for obstime in MM-DD-YYYY
        for quasar in all_pos:
            #AB convert to altitude and azimuth
            altaz = AltAz(location=loc1, obstime=time)
            quasaraa = quasar.transform_to(altaz)
            Qairmass = quasaraa.secz #find airmass of quasar
            #AB given that the quasar's airmass is positive (above horizon), and
            #we find the smallest airmass:
            if Qairmass < best_airmass and quasaraa.alt > 0 * u.deg:
                best_airmass = Qairmass
                best_quasar = quasar
                best_time = time
    #print
    if best_quasar:
        print(f"The optimal quasar to observe is located at:" f" {best_quasar.to_string('hmsdms')}" f" ,At time {best_time.iso} UTC" f" ,with an airmass of {best_airmass}")
        return



if __name__ == '__main__':

    dates1 = ["2020-01-01 14:00", "2021-01-01 14:00", "2022-01-01 14:00", "2023-01-01 14:00", "2024-01-01 14:00","2025-01-01 14:00","2026-01-01 14:00","2027-01-01 14:00","2028-01-01 14:00","2029-01-01 14:00","2030-01-01 14:00","2020-01-02 02:00", "2021-01-02 02:00", "2022-01-02 02:00", "2023-01-02 02:00", "2024-01-02 02:00","2025-01-02 02:00","2026-01-02 02:00","2027-01-02 02:00","2028-01-02 02:00","2029-01-02 02:00","2030-01-02 02:00"]
    planets1 = ["mercury", 'venus', 'mars', 'jupiter', 'saturn']
    planets_ecliptic(planets1, dates1)
    planets_equatorial(planets1, dates1)
    vis_qso(1,2024)
