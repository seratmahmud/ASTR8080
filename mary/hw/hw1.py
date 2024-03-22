# M. Kaldor
# v1 2/1/2024
# ASTR 8080 HW1

# to import this from another directory:
# import sys
# sys.path.insert(0, '../week1/')
# import hw0
# from hw0 import function

# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib
matplotlib.use("MacOSX")
import matplotlib.pyplot as plt
from matplotlib import rc
import time
import pdb
import astropy
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import sfdmap
import extinction


# FUNCTIONS
###############################
###############################
def planet_positions_ecliptic(planets, begin, end, times):
    """
    This is a function that takes in a list of planets, years, and times of day and produces a plot of their
    locations at those times in ecliptic coordinates.

    Args: planets = [list of strings]
        begin = int
        end = int
        times = [list of ints]

    Returns: list of dates as SkyCoord objects, planet positions as astropy coordinate objects

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: Visualizing planet positions for any time with respect to the ecliptic
    #
    #
    # CALLING SEQUENCE: planet_positions_ecliptic(["mercury", "venus", "mars", "jupiter", "saturn"], 2020, 2030, [7, 19])
    #
    #
    # INPUTS: planets = planet names
    #   begin = year for starting
    #   end = year for ending
    #   time = time of day in 24-hour reference
    #
    # -------------------------------------------------------------
    #
    # Establish astropy time objects in the UTC frame for all days and all years
    year = begin
    datelist = []
    while year < end+1:
        for time in times:
            if time+7<24:
                stringtime = str(year)+"-01-01 "+str(time+7)+":00:00"
            else:
                stringtime = str(year)+"-01-02 " + str(time+7-24)+":00:00"
            timetime = astropy.time.Time(stringtime)
            datelist.append(timetime)
        year += 1


    # Initialize the planet positions array by making it the length of however many planets you have and the height of
    # however many dates you will have in total
    planet_pos = np.zeros([len(planets), len(times) * (end - begin + 1)], dtype=object)
    i = 0
    for planet in planets:
        # Do all ecliptic conversions for one planet at a time
        planet_list_ec = []
        for date in datelist:
            planloc = astropy.coordinates.get_body(planet, date)
            planet_ec = planloc.transform_to("heliocentrictrueecliptic")
            planet_list_ec.append(planet_ec)
        planet_list_ec = np.array(planet_list_ec)
        # Place the planet information into the planet positions array
        planet_pos[i,:] = planet_list_ec
        i += 1

    # Make color list to iterate through when plotting planets so they're not all the same color
    colorlist = ["blue", "red", "yellow", "green", "purple", "orange", "pink", "black", "gray"]

    # Pull out latitude and longitude from all planets, plot them in a scatter plot.
    for i in np.arange(0, len(planets)):
        lon = [planet_pos[i, ind].lon.deg for ind in np.arange(0, len(times)*(end-begin+1))]
        lat = [planet_pos[i, ind].lat.deg for ind in np.arange(0, len(times)*(end-begin+1))]
        plt.scatter(lon, lat, color=colorlist[i], label=planets[i])
    plt.legend()
    plt.xlabel("Longitude (degrees)")
    plt.ylabel("Latitude (degrees)")
    plt.title("Latitude vs. Longitude")
    plt.show()

    return datelist, planet_pos

###############################
def planet_positions_equatorial(planets, begin, end, times):
    """
    This is a function that plots the planet positions in RA and DEC over a span of years at multiple times of day and returns their
    airmasses and marks the lowest airmass for each planet on a plot.

    Args: planets = [list of strings]
        begin = int
        end = int
        times = [list of ints]

    Returns: list of dates as SkyCoord objects, planet positions as astropy coordinate objects

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: Visualizing planet positions for any time in RA and DEC and finding the best observing targets (those
    # with low airmass)
    #
    #
    # CALLING SEQUENCE: planet_positions_equatorial(["mercury", "venus"], 2020, 2030, [7, 19])
    #
    #
    # INPUTS: planets = planet names
    #     begin = year for starting
    #     end = year for ending
    #     time = time of day in 24-hour reference
    #
    # -------------------------------------------------------------
    #
    # Establish astropy time objects in the UTC frame for all days and all years
    year = begin
    datelist = []
    while year < end + 1:
        for time in times:
            if time + 7 < 24:
                stringtime = str(year) + "-01-01 " + str(time + 7) + ":00:00"
            else:
                stringtime = str(year) + "-01-02 " + str(time + 7 - 24) + ":00:00"
            timetime = astropy.time.Time(stringtime)
            datelist.append(timetime)
        year += 1

    # Set our location to APO
    obs_loc = astropy.coordinates.EarthLocation.of_site("apo")
    # Initialize the planet positions array by making it the length of however many planets you have and the height of
    # however many dates you will have in total
    planet_pos = np.zeros([len(planets), len(times) * (end - begin + 1)], dtype=object)
    i = 0
    # Make lists to track the minimum airmasses that we will find for planets
    airmassmin_list = []
    airmassdate_list = []
    airmassmin_ind = []
    for planet in planets:
        planet_list_eq = []
        airmassmin = 100
        airmassdate = astropy.time.Time.now()
        j = 0
        # Find airmass by taking secant of the zenith angle after putting the planet into altitude-azimuth frame
        for date in datelist:
            planloc = astropy.coordinates.get_body(planet, date)
            planet_eq = planloc.transform_to("icrs")
            planet_list_eq.append(planet_eq)
            altaz = planloc.transform_to(astropy.coordinates.AltAz(obstime=date, location=obs_loc))
            airmass = altaz.secz
            # Check that the airmass is a physically relevant (positive) value, then note down if the value is the
            # minimum overall per planet
            if 0 < airmass < airmassmin:
                airmassmin=airmass
                airmassdate=date
                airmassind=j
            j += 1
        airmassmin_list.append(float(airmassmin))
        airmassdate_list.append(airmassdate)
        airmassmin_ind.append(airmassind)
        planet_list_eq = np.array(planet_list_eq)
        # Place the planet information into the planet positions array
        planet_pos[i, :] = planet_list_eq
        i += 1

    print("Airmass minimum dates: ", airmassdate_list)
    print("Airmass minimum values: ", airmassmin_list)

    # Make color list to iterate through when plotting planets so they're not all the same color
    colorlist = ["blue", "red", "yellow", "green", "purple", "orange", "pink", "black", "gray"]

    # Pull out RA and DEC from all planets, plot them in a scatter plot. Overlay the overall minimum for each planet
    # in yellow and label it.
    fig, ax = plt.subplots()
    for i in np.arange(0, len(planets)):
        ra = [planet_pos[i, ind].ra.deg for ind in np.arange(0, len(times) * (end - begin + 1))]
        dec = [planet_pos[i, ind].dec.deg for ind in np.arange(0, len(times) * (end - begin + 1))]
        plt.scatter(ra, dec, color=colorlist[i], label=planets[i])
        plt.scatter(ra[airmassmin_ind[i]], dec[airmassmin_ind[i]], color="yellow", marker="x", s=30)
        ax.annotate(text="X_min "+str(planets[i]), xy=(ra[airmassmin_ind[i]], dec[airmassmin_ind[i]]))
    plt.legend()
    plt.xlabel("RA (degrees)")
    plt.ylabel("DEC (degrees)")
    plt.title("DEC vs. RA")
    plt.show()

    return datelist, planet_pos

###############################
def quasar_airmass(month, year):
    """
    This is a function that takes in a month and year and tells you which quasar is the best observing target (lowest
    airmass source) is on each day of that month

    Args: month = int
        year = int

    Returns: list of airmasses for targets, RA and DEC of targets, and list indices of targets from quasar list

    """
    # -------------------------------------------------------------
    # +
    # PURPOSE: selecting the best quasars to view ahead of time (like if you were to go observing soon!)
    #
    #
    # CALLING SEQUENCE: quasar_airmass(2, 2024)
    #
    #
    # INPUTS: month = integer for observing month with 12-month reference
    #   year = year for observing
    #
    # -------------------------------------------------------------
    # Read in quasar file and turn it into usable SkyCoord objects
    datafile = "/Users/marykaldor/ASTR8080/runnoe/week3/HW1quasarfile.dat"
    data = np.genfromtxt(datafile, unpack=True, dtype=str)
    coordlist = []
    for q in data:
        coord = q[0:2] + " " + q[2:4] + " " + q[4:9] + " " +q[9:12] + " " + q[12:14] + " " + q[14:18]
        coord = SkyCoord(coord, unit=(u.hourangle, u.deg))
        coordlist.append(coord)
    datelist = []
    # Tell yourself what you think you grabbed (just a checkpoint)
    print("Month is", month)
    print("Year is", year)

    # Set month number of days based on what month (and year) it is - don't forget about leap years! Make astropy time
    # objects for every day of the month in UTC
    if month==1 or month==3 or month==5 or month==7 or month==8 or month==10 or month==12:
        stringtime = ["2024-" + str(month).zfill(2) + "-" + str(num).zfill(2) + " 23:00:00" for num in np.arange(1,32)]
        print("This is a 31 day month")
        timetime = [astropy.time.Time(string) + 7*u.hour for string in stringtime]
        datelist.append(timetime)

    if month==4 or month==6 or month==9 or month==11:
        stringtime = ["2024-" + str(month).zfill(2) + "-" + str(num).zfill(2) + " 23:00:00" for num in np.arange(1, 31)]
        print("This is a 30 day month")
        timetime = [astropy.time.Time(string) + 7*u.hour for string in stringtime]
        datelist.append(timetime)

    if month==2 and year%4 != 0:
        stringtime = ["2024-" + str(month).zfill(2) + "-" + str(num).zfill(2) + " 23:00:00" for num in np.arange(1, 29)]
        print("This is a 28 day month - not a leap year")
        timetime = [astropy.time.Time(string) + 7*u.hour for string in stringtime]
        datelist.append(timetime)

    if month==2 and year%4 == 0:
        stringtime = ["2024-" + str(month).zfill(2) + "-" + str(num).zfill(2) + " 23:00:00" for num in np.arange(1, 30)]
        print("This is a 29 day month - a leap year")
        timetime = [astropy.time.Time(string) + 7*u.hour for string in stringtime]
        datelist.append(timetime)

    # Set our location to APO
    obs_loc = astropy.coordinates.EarthLocation.of_site("apo")

    # I think that this might (?) be helpful, but I couldn't get it to work. This section isn't used but might inspire
    # me in the future when I look back on this!
    dates = np.array(datelist[0])
    coords = np.array(coordlist)
    #altaz = coords.transform_to(astropy.coordinates.AltAz(obstime=dates, location=obs_loc))
    #print(altaz)
    # quasar_dates = np.recarray(dates, coords, dtype=object)
    # print("quasar dates", np.shape(quasar_dates))

    # List for minimum value indices
    mins = []

    # Find airmass by taking secant of the zenith angle after putting the planet into altitude-azimuth frame
    print("Warning! This takes a super long time! Go get a coffee. Sorry!")
    for date in dates:
        altaz = [quasar.transform_to(astropy.coordinates.AltAz(obstime=date, location=obs_loc)) for quasar in coordlist]
        airmass = [q.secz for q in altaz]
        airmass = np.array(airmass)
        # Check that the airmass is a physically relevant (positive) value, then note down if the value is the
        # minimum overall per day across all quasars
        valid_idx = np.where(airmass >= 0)[0]
        out = valid_idx[airmass[valid_idx].argmin()]
        mins.append(out)

    # Return helpful values for picking targets to observe
    print("Indices of targets: ", mins)
    airmasses = [airmass[ind] for ind in mins]
    print("Airmasses of targets: ", airmasses)
    targets = [coordlist[ind] for ind in mins]
    print("RA and DEC of targets: ", targets)

    return airmasses, targets, mins






# MAIN
###############################
###############################
if __name__ == '__main__':
    # MEK set the timer
    time0 = time.time()

    planet_positions_ecliptic(["mercury", "venus", "mars", "jupiter", "saturn"], 2020, 2030, [7, 19])

    planet_positions_equatorial(["mercury", "venus"], 2020, 2030, [7, 19])

    quasar_airmass(2, 2024)

    # MEK time
    time1 = time.time()
    print(" Time [mus]: {0:.3f}".format(1e6 * (time1 - time0)))
