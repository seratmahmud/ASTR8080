# M. Rizzo Smith
# v1 1/25/24
# ASTR 8080 HW1


# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import rc
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Galactic
from astropy.coordinates import get_body
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates import AltAz

# MRS I have this warning supression here to stop getting that pesky dubious year error
import warnings
warnings.filterwarnings('ignore')

# FUNCTIONS
###############################
###############################

def generate_epochs(start_time):
    """
    This is a function to generate a set of observation epochs based on a given start date and time.
    It will take the initial date and create an array of astropy.Time objects containing that date for the next 10 years at the entered time and 12 hours later.
    
    Args:
        start_time (str): Date and time in format YYYY-MM-DD 00:00:00
    
    Returns:
        all_epochs (numpy array): Array of all observation epochs
    """

#------------------------------------------
#+
# Purpose:
#   Return an array of observation epochs beginning at a specified time and date.
# 
# CALLING SEQUENCE:
#   all_epochs = generate_epochs('YYYY-MM-DD 00:00')
#
# INPUTS:
#   start_time - 'string date time'
#
#-
#------------------------------------------


    # MRS Pull the year out of the start date
    year = int(start_time[:4])
    
    # MRS initialize an array to store the dates
    date_array = []
    # MRS Generate the next 10 years of dates
    for i in range(11):
        # MRS Append the next year's date
        date_array.append(f'{year + i}-{start_time[5:]}')

    # MRS creat astropy.Time objects for each 7AM epoch
    morning = Time(date_array, format='iso', scale='utc')

    # MRS Add half a Julian Day to shift 12 hours for 7pm epochs
    night = morning.jd + 0.5
    # MRS Make the 7pm epochs a Time object
    night = Time(night, format='jd')
    night.format = 'iso'
    
    # MRS Loop through and combine all 7am and 7pm epochs into one long array
    all_epochs = []
    for i in range(0, len(morning)):
        all_epochs.append(morning[i])
        all_epochs.append(night[i])
    
    return np.array(all_epochs)


def planet_ecliptic(planets, date):
    """
    This is a function to plot the positions of a specified list of planets over a given range of dates in eccliptic coordinates.

    Args:
        planets (string/array of strings): List of planets
        dates (numpy array): All desired observation dates

    Returns:
        all_pos (numpy array): Array of SkyCoord objects of each planet on each date specified
    """

#-----------------------------------------
#+
#
# PURPOSE:
#   Calculate positions of specified planets on a given range of dates. 
#   Provide a plot of each planets positon in eccliptic coordinates
#
# CALLING SEQUENCE:
#   positions = planet_eccliptic(planets, dates)
#
# INPUTS:
#   planets - ['Mercury', 'Venus', ...]
#   dates - ['Observation Dates']
#
#-
#-------------------------------------------
    dates = generate_epochs(date)
    # MRS Initialize array that will hold positions for all planets
    all_pos = np.empty((len(planets), int(len(dates))), dtype=object)
    
    # MRS Loop through each planet
    for i, planet in enumerate(planets):
        # MRS Loop through each observation epoch
        for j, date in enumerate(dates):
            # MRS Store the position of a given planet for each epoch
            all_pos[i,j] = get_body(planet, date)
            all_pos[i,j] = all_pos[i,j].transform_to('heliocentrictrueecliptic')
    # MRS Now pull the ra and dec in degrees and plot each planet
    # MRS Preset colors for each planet (WONT WORK IF I KEEP LIKE THIS)
    colors = plt.cm.tab10(range(len(planets))) 

    # MRS Need to pull ra and dec out of my big array of positions
    # MRS This is so that I can plot all of the planets without pulling ra and dec from each cell of the all_pos array.
    for i, planet in enumerate(planets):
        lat = [temp.lat.degree for temp in all_pos[i, :]]
        lon = [temp.lon.degree for temp in all_pos[i, :]]

        plt.plot(lon, lat,'o', color = colors[i], label=f'{planet}')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title(f'Yearly Planetary Positions from \n {dates[0]} to {dates[-1]} ')
    plt.legend()
    plt.show()
    return all_pos
    
    
def planet_equatorial(positions, planets, date):
    """
    This is a function to plot the positions in equatorial coordinates of a specified list of planets for a given range of observations.
    This function will print out the date of lowest airmass for each planet and add a visually distinct point on the plot for this date.
    
    Args:
        positions (numpy array): An array of SkyCoord objects for each planet on each epoch of the date range. THIS IS AN OUTPUT FROM planet_eccliptic() ABOVE
        planets (string/array of string): List of planets
        dates (numpy array): All desired observation dates
    
    Returns:
        N/A
        Does print lowest airmass for each planet 
        Produces plot of equatorial planet positions
    """
#--------------------------------------------------
#+
# PURPOSE:
#   Convert positions of specified planets into equatorial coordinates on a given range of dates. 
#   Provide a plot of each planets positon, print out the lowest airmass and date of lowest airmass, add visually distinct point to plot for lowest airmass
#   The positions array can be easily obtained by calling planet_eccliptic from above :D
#
# CALLING SEQUENCE:
#   planet_equatorial(positions, planets, dates)
#
# INPUTS:
#   positions - [Skycoord object of planet]
#   planets - ['Mercury', 'Venus', ...]
#   dates - ['Observation Dates']
#
#-
#---------------------------------------------------
    # MRS Generate range of next 10 years from start date
    dates = generate_epochs(date)
    # MRS Set APO location 
    APO = EarthLocation.of_site('APO')
    # MRS Initialize Alt Az frames from APO for each date
    all_aa_frames = AltAz(location=APO, obstime=dates)

    # MRS Initialize array to hold planet alt az coordinates
    aa_coords = np.empty(positions.shape, dtype=object)
    
    # MRS Initialize array to hold airmass at each observation
    airmass = np.empty(positions.shape, dtype=object)

    # MRS Initialize array to hold lowest airmass epoch for each planet
    lowest_airmass = np.empty((len(planets),3), dtype=object)

    # MRS Loop through the each planet and each observation convert to ICRS for RA and DEC
    # MRS Loop through each planet and store alt az coordinates 
    for i in range(0, positions.shape[0]):
        for j  in range(0, positions.shape[1]):
            positions[i,j] = positions[i,j].transform_to('icrs')
            aa_coords[i,j] = positions[i,j].transform_to(AltAz(location=APO, obstime=dates[j]))
            # MRS Calculate the airmass for each observation
            airmass[i,j] = aa_coords[i,j].secz
    
        # MRS Only consider the lowest airmass above the horizon (positive airmass)
        positive_air = airmass[i, airmass[i, :] > 0]

        # MRS Loop through and find the minimum positive airmass for each planet
        if positive_air.size >0:
            # Store the index from the positive_air array of the min airmass
            min_index = np.argmin(positive_air)
            # Store the minimum positive airmass 
            lowest_airmass[i,0] = positive_air[min_index].value
            # MRS Keep the index of the positions array at which this lowest airmass occurs
            lowest_airmass[i,1] = np.where(airmass[i, :] ==  positive_air[min_index])[0][0]
            # MRS Record the date for each planet on which the lowest airmass occurs
            lowest_airmass[i,2] = dates[lowest_airmass[i,1]]
            
            # MRS Print out each planet, its lowest airmass, and the date
            print(f'At APO, {planets[i]} Lowest Airmass sec(z):', lowest_airmass[i,0], 'On:', lowest_airmass[i,2], 'UTC')
    
    # MRS Select some cool color map to fit all your planets
    colors = plt.cm.tab10(range(len(planets))) 

    # MRS Plot all positions planet by planet
    for i in range(0, positions.shape[0]):
        ra = [temp.ra.degree for temp in positions[i, :]]
        dec = [temp.dec.degree for temp in positions[i, :]]
        ra_low = positions[i, lowest_airmass[i, 1]].ra.degree
        dec_low = positions[i, lowest_airmass[i, 1]].dec.degree
        
        plt.plot(ra, dec,'o', color = colors[i], label=f'{planets[i]}')
        # MRS Overlay a star for the lowest airmass epochs
        plt.plot(ra_low, dec_low, '*', color=colors[i], markeredgecolor='black', markersize=10)
    # MRS Plot a dummy nan coordinate to add a label to my legend
    plt.plot(np.NaN, np.NaN, '*', color='black', label='Lowest Airmass')
    plt.xlabel('Right Asencion')
    plt.ylabel('Declination')
    plt.title(f'Yearly Planetary Positions from \n {dates[0]} to {dates[-1]} ')
    plt.legend()
    plt.show()
    
    return

def vis_quasars(month):
    """
    This is a function that given a month of the year 2024 will return what Quasar, from a list of objects, is the most optimal target
    It will return the RA and Dec as well as date for the lowest airmass observation of ALL targets.

    Args:
        month (float of integer): If entering a leading zero for single digit numbers you must add a decimal to make the entry a float ex. 09. additionally an integer would work just as well ex. 9.

    Returns:
        N/A Does print the optimal target and date of observation for given month
    """
#-------------------------------------------
#+
# PURPOSE:
#   Return the optimal observation date of the lowest airmass target for a specified month, given a list of targets.
#
# CALLING SEQUENCE:
#   vis_quasars(month)
#
# INPUTS:
#   month - 09.
#
#-
#---------------------------------------------------

    # Initialize a list to hold the quasar positions read in from quasar_file
    all_pos = []
    quasar_file = 'HW1quasarfile.dat'

    # MRS Loop through the quasar file and split the string into 
    # MRS the RA and Dec format 'HH MM SS.SS +DD MM SS.SS' 
    with open(quasar_file, 'r') as f:
        for line in f:
            pos_str = line
            pos_str = f'{pos_str[:2]} {pos_str[2:4]} {pos_str[4:9]} {pos_str[9:12]} {pos_str[12:14]} {pos_str[14:]}'
            
            # MRS Turn the ra dec string into a SkyCoord object
            pos = SkyCoord(pos_str, unit=(u.hourangle, u.deg))
            # MRS Add that quasar to the list of all quasars
            all_pos.append(pos)
    all_pos = np.array(all_pos)
    
    # MRS Determine the number of days in the given month, not always 31 or 20
    if 1 <= month <=12:
        if month in {1,3,5,7,8,10,12}:
            num_days = 31
        elif month in {4,6,9,11}:
            num_days = 30
        else:
            num_days = 28

    # MRS With the number of days create an array for each day of the month
    full_month = np.linspace(1, num_days, num_days) 
    # MRS For each day in the array, add that day into a full date string with
    # MRS the format 'YYYY-MM-DD HH:MM:SS' 
    full_month = [f'2024-{str(int(month)).zfill(2)}-{int(day):02} 23:00:00' for day in full_month]
    # MRS Convert each daily observation into a time object in MST
    all_epochs = Time(full_month, format='iso', scale='utc')
    # MRS Convert to UTC by adding 7 hours
    all_epochs_utc = all_epochs + 7*u.hour
    
    # MRS Set APO as the observing location
    APO = EarthLocation.of_site('APO')
    # MRS Initialize an empty array that will hold the airmass of every quasar on each day of the mont
    airmass = np.empty((len(all_pos),len(all_epochs_utc)))
    # MRS Create a large value such that we can determine the largest positive airmass later and save it
    lowest_air = 10000000
    # MRS index 0 of this arry will return the quasar position, and index 1 will return the date of observation
    lowest_air_index = [0, 0]
    
    # MRS Warn user that retrieving best target may take time
    print('Retrieving Optimal Target Quasar, may take a bit...')
    # MRS Loop through all of the quasars, and calcualte the airmas on each day of the month
    for i, quasar in enumerate(all_pos):
        for j, epoch in enumerate(all_epochs_utc):
            # MRS Transform from RA and Dec to Alt Az for each observation epoch
            all_pos_aa = quasar.transform_to(AltAz(location=APO, obstime=epoch))
            # MRS Calcualte the airmass on a given date
            airmass[i,j] = all_pos_aa.secz.value
            # MRS If statement to see if this is the lowest positive airmass we have encountered yet
            if 0 < airmass[i,j] < lowest_air:
                # MRS Save the lowest airmass value
                lowest_air = airmass[i,j]
                # MRS Save which quasar (i), and which date (j) that this occurs
                lowest_air_index = [i, j]
    # MRS Super long print statement to return which quasar is the optimal target and when
    print(f'For the month {str(int(month)).zfill(2)} the optimal quasar observation is Object located at: RA', all_pos[lowest_air_index[0]].ra.degree,'Dec', all_pos[lowest_air_index[0]].dec.degree, 'On:', all_epochs[lowest_air_index[1]], 'MST from APO')
    return 

def combined_tasks(date, planets, month):
    # MRS Task 1 Generate the planetary position plot in ecliptic coords
    pos = planet_ecliptic(planets, date)

    # MRS Task 2 Generate the planetary positions plot in equatorial coords
    # MRS Task 2 also prints the best observation date for each planet based on lowest airmass
    planet_equatorial(pos, planets, date)
    
    # MRS Task 3 Print out the optimal target for a given month from a list of quasars
    vis_quasars(month)

    return
    
# MAIN
###############################
###############################
if __name__ == '__main__':
    # MRS Makes sure the user enters a sufficient number of arguments
    if len(sys.argv) <3:
        print("Usage: python HW1.py 'YYYY-MM-DD HH:MM' 'Planet' 'Month' ")
        print('User must enter atleast 1 planet to run script')
    else:
        # MRS Start date for 10 year range
        date = sys.argv[1]
        # MRS Planets to plot
        planets = sys.argv[2:-1]
        # MRS Month to sample best quasar observation for 2024
        month =float(sys.argv[-1])
        # MRS Do all three tasks with parameters from the command line
        combined_tasks(date, planets, month)

