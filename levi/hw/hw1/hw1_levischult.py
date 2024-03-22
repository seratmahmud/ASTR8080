# L. Schult
# v1 31 Jan 2024
# ASTR 8080 
# to import this from another directory:
# import sys
# sys.path.insert(0, '../week1/')
# import polycalc
# from polycalc import get_poly_o3

# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import astropy
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body
import pdb
import time
import calendar

# FUNCTIONS
###############################
###############################
def testname(x):
    """
    This is a function that 

    Args:
        x (datatype): A...
    Returns:
        y (datatype): A ...
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   You pass butter
#
# CALLING SEQUENCE:
#   variable = testname(x) 
#
# INPUTS:
#   x
#-
#-------------------------------------------------------------  
    y     = x**2. 
    return y 

def generatedates(startdate, nyears=10):
    """
    This is a function that generates dates over a given year range or for every
    day of the month. The time specified in the start date will be used for every 
    date time generated. To do a month, set nyears=0. If multiple times a day for
    a month, the year and month that will be used are chosen from startdate's 
    first entry. NOTE: assumes MST (UTC - 7 hours)

    date format: 'yyyy-mm-ddT20:00:00'

    Args:
        startdate (list): A list of datetime strings
        nyears (int): number of years to generate data for. set to 0 for every day of month
    Returns:
        utc_datetimes (array): An array of generated dates as an array.
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   generate repeated days for number of years at same given time/day or same 
#   time every day for a month
#
# CALLING SEQUENCE:
#   utctimes = generatedates(['2020-01-01T20:00:00'], nyears=4) 
#
# INPUTS:
#   startdate
#   nyears
#-
#-------------------------------------------------------------  
    datetimes = []
    
    # LSS generating N years of data on same date
    if nyears >= 1:
        startyr = int(startdate[0][:4])
        for i in range(startyr, startyr+nyears+1): # LSS taking care of [0, 1) type things
            for date in startdate: # LSS iterating over dates given
                datetimes.append((str(i))+date[4:]) # LSS using string indexing
                # LSS 1 year is not a calendar year!
    
    # LSS generating dates of one month for given year/month
    elif nyears == 0:
        workdate = startdate[0] # LSS assumes dates given in list are the same year and month, just different times
        yr = int(workdate[:4])
        mth = int(workdate[5:7])
        # LSS utilizing calendar functions (0 days are outside given month)
        cal = calendar.Calendar()
        monthdays = np.array([day for day in cal.itermonthdays(year=yr, month=mth)])
        monthdays = monthdays[monthdays > 0]
        for day in monthdays:
            for date in startdate: # LSS iterating over dates given
                if len(str(day)) < 2:
                    stringday = '0' + str(day)
                    datetimes.append(date[:8] + stringday + date[10:])
                else:
                    stringday = str(day)
                    datetimes.append(date[:8] + stringday + date[10:])

    # LSS convert the times to UTC
    #for i in datetimes: print(i)
    mdts = Time(datetimes, format='isot', scale='utc') # LSS mountain time
    udts = mdts + 7 * u.hour # LSS convert to UTC

    return udts

def plotplanets(datetime, nyears=10, coordframe='ecliptic', airmass=True):
    """
    This is a function that plots the locations of Mercury, Venus, Mars, Jupiter
    and Saturn for given dates + times in a given timezone, in given coordinates.
    If airmass=True (default), then the date/time of lowest airmass for each object
    will be printed as well as plotted with a bigger marker in the plot. NOTE: 
    assumes MST (UTC - 7 hours)

    Args:
        datetime (list): list of datetimes to plot the planets at. Should be in
        format: 'yyyy-mm-ddT20:00:00'
        nyears (int): number of years to generate datetimes for. default is 10.
        if nyears=0, then every day of the month from the first datetime entry
        will be generated for the times specified. 
        coordfram (str): 'ecliptic', 'equatorial' are the options here - for plotting
        airmass (bool): whether or not to find datetime of lowest airmass for each
        object.

        
    Returns:
        udts (array): an array of the datetimes the coords are generated from.
        planetcoord (array): array of planet coordinates in planets x datestimes
        planet order is [mercury, venus, mars, jupiter, saturn]
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   plot first 5 non-earth planets for user specified dates, times at a given 
#   timezone + coordinate frame. Also return the datetime with lowest airmass.
#
# CALLING SEQUENCE:
#   datetimes, planetcoords = plotplanets([2022-02-01T23:00:00], coordframe='equatorial', airmass=True) 
#
# INPUTS:
#   datetimes
#   coordframe (optional)
#   airmass (optional)
#-
#-------------------------------------------------------------  
    # LSS generate dates
    udts = generatedates(datetime, nyears=nyears)
    
    # LSS getting the location of planets
    loc = EarthLocation.of_site('Apache Point Observatory') # LSS setting location
    with solar_system_ephemeris.set('builtin'):
        merc = get_body('mercury', udts, loc)
        ven = get_body('venus', udts, loc)
        mars = get_body('mars', udts, loc)
        jup = get_body('jupiter', udts, loc) 
        strn = get_body('saturn', udts, loc)
    
    # LSS setting SkyCoords to ecliptic coords
    if coordframe=='ecliptic':
        coordtotrans = 'heliocentrictrueecliptic' 
    elif coordframe=='equatorial':
        coordtotrans = 'icrs'
    merc = merc.transform_to(coordtotrans)
    ven = ven.transform_to(coordtotrans)
    mars = mars.transform_to(coordtotrans)
    jup = jup.transform_to(coordtotrans)
    strn = strn.transform_to(coordtotrans)

    # LSS plotting
    planet_coords = [merc, ven, mars, jup, strn]
    planet_names = ['mercury', 'venus', 'mars', 'jupiter', 'saturn']
    mrkrs = ['s', 'p', 'x', '+', '*']
    # LSS making size arrays that will make lowest airmass point larger.
    # LSS if airmass is false, all points will be the same size
    if airmass==False:
        # LSS array of marker sizes for plotting. when airmass is false all equal
        sizearr = np.zeros((len(planet_coords), len(planet_coords[0])))
        sizearr.fill(8.0)
        sizearr = list(sizearr)
    else:
        # LSS making array to hold the airmasses for each object at each time
        planet_am = np.zeros((len(planet_coords), len(planet_coords[0])), dtype='object')
        # making isotropic size array initially.
        sizearr = np.zeros((len(planet_coords), len(planet_coords[0])))
        sizearr.fill(8.0)
        sizearr = list(sizearr)

        for idx in range(len(planet_coords)): # LSS transform to altaz + get airmass
            planet_am[idx, :] = planet_coords[idx].transform_to(
                astropy.coordinates.AltAz(obstime=udts, location=loc)).secz
            # LSS restrict airmasses to those above horizon of observatory
            abovehorizon = np.where(planet_am[idx] > 0)[0]
            # LSS find lowest positive airmass
            # LSS indexing by above horizon restricts to positive values
            # LSS argmin finds the minumum and returns the index of abovehorizon
            # LSS that is the minimum, and then re-indexes above horizon to get 
            # LSS the actual index of planet_am that is the lowest airmass.
            lowest_am_index = abovehorizon[planet_am[idx][abovehorizon].argmin()]
            print(f'\n{planet_names[idx]} has lowest airmass on {udts[lowest_am_index]} UTC, airmass={planet_am[idx][lowest_am_index]} \n')
            # LSS added new lines for ease of reading + warning messages clog up 
            # LSS command line

            # LSS changing sizearr to make lowest airmass date larger on plot
            sizearr[idx][lowest_am_index] = 60.0
        

    # LSS handling coordinate frame names/units in plot
    if coordframe=='ecliptic':
        for idx, name in enumerate(planet_names):
            plt.scatter(planet_coords[idx].lon.deg, planet_coords[idx].lat.deg, 
                        label=name, s=sizearr[idx], marker=mrkrs[idx])
        plt.xlabel(coordtotrans+' longitude (º)')
        plt.ylabel(coordtotrans+' latitude (º)')

    elif coordframe=='equatorial':
        for idx, name in enumerate(planet_names):
            plt.scatter(planet_coords[idx].ra.deg, planet_coords[idx].dec.deg, 
                        label=name, s=sizearr[idx], marker=mrkrs[idx])
        plt.xlabel(coordframe+' RA (º)')
        plt.ylabel(coordframe+' Dec (º)') 
    
    firstdt = udts[0]
    lastdt = udts[-1]
    plt.title(f'Locations of first 5 non-Earth planets from \n {firstdt} to {lastdt}', fontsize=10)
    plt.legend()
    plt.show()

    return udts, np.array(planet_coords)

def monthofquasarobs(datetime):
    """
    This is a function that finds the quasar with the lowest airmass on every day
    at a given time over the course of a month In the datetime argument, 
    the year, month, and time matter. the day does not have to be the first.
    NOTE: assumes MST (UTC - 7 hours)


    date format: 'yyyy-mm-ddT20:00:00'

    Args:
        datetime (list): list of strings to specify the month and time of interest, 
        year is also taken into account. NOTE: the year and month are taken from 
        the first entry in this list, any additional entries other than the first
        are to specify other TIMES only.
    Returns:
        quasarskycoords (array of SkyCoord objects): quasars in equatorial 'icrs'
        format.
    """
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   find the lowest airmass quasar at a certain time for every day of a given month
#
# CALLING SEQUENCE:
#   variable = monthofquasarobs(['2020-01-01T20:00:00']) 
#
# INPUTS:
#   datetime
#-
#-------------------------------------------------------------  
    # LSS generate one month of dates for given year + month
    udts = generatedates(datetime, nyears=0)
    print('\n quasar month observation dates generated \n')
    # LSS load in data + parse raw text into skycoords
    with open('./HW1quasarfile.dat', 'r') as infile:
        rawtxt = np.loadtxt(infile, dtype=str)
    coordstrings = [] 
    strtim = time.time()
    for line in rawtxt:
        coordstrings.append(f'{line[:2]}h{line[2:4]}m{line[4:9]} {line[9:12]}d{line[12:14]}m{line[14:]}s')
    quasarcoords = SkyCoord(coordstrings)

    # LSS getting airmasses of quasars
    loc = EarthLocation.of_site('Apache Point Observatory') # LSS setting location
    quasar_am = np.zeros((len(quasarcoords), len(udts))) # creating array to hold airmasses
    for idx, qso in enumerate(quasarcoords):
        quasar_am[idx, :] = qso.transform_to(astropy.coordinates.AltAz(obstime=udts, location=loc)).secz

    # LSS going to iterate over each day of the month and find the quasar with 
    # the lowest airmass. That quasar will be printed.
    for idx, day in enumerate(udts):
        abovehorizon = np.where(quasar_am[:, idx] > 0)[0] # iterate over days, not quasars
        '''the below does something similar to plotplanets. It first selects only
        positive airmasses with np.where, becoming the abovehorizon variable. We
        then find the airmass index that has the lowest value when masked for abovehorizon
        values using argmin. We use this index to return the above horizon value
        at that index which is the index of the true quasar_am array, hence finding
        the quasar with the lowest airmass on that day
        '''
        quasar_lowam_idx = abovehorizon[quasar_am[:, idx][abovehorizon].argmin()]
        print(f'\n At {udts[idx]} UTC the quasar with lowest airmass is {quasarcoords[quasar_lowam_idx]}, airmass={quasar_am[quasar_lowam_idx, idx]} \n')

    return quasarcoords


# MAIN
###############################
###############################
if __name__ == '__main__':
    strtim = time.time()
    timesforq12 = ['2020-01-01T07:00:00', '2020-01-01T19:00:00'] # LSS time in mountain time
    # LSS below is the function call for problem 1
    _ = plotplanets(timesforq12)
    
    # LSS below is the function call for problem 2
    _ = plotplanets(timesforq12, coordframe='equatorial', airmass=True)

    # LSS below is the function call for problem 3
    quasardate = ['2022-02-01T23:00:00']
    _ = monthofquasarobs(quasardate)
    endtimes = time.time() - strtim
    print(f'total runtime: {endtimes}s')
    #print('hello world XD')

