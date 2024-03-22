# S. Saad
# v1 1/19/2024
# ASTR 8080 HW1

# IMPORT BLOCK
###############################
###############################
from astropy.coordinates import get_body, EarthLocation, AltAz
from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import HeliocentricTrueEcliptic, SkyCoord
from astropy.time import Time
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import warnings
from astropy.utils.exceptions import ErfaWarning
warnings.filterwarnings('ignore', category=ErfaWarning)

# FUNCTIONS
###############################
###############################


# SS Function to get plot the planets position in ecliptic coordinates
def plot_planet_positions_minimal_loops():
    # SS Calling all the planets , years, times of day, and other infos
    planets = ['mercury', 'venus', 'mars', 'jupiter', 'saturn']
    years = np.arange(2020, 2031)
    times_of_day = ['07:00', '19:00'] # SS MST Time
    colors = ['grey', 'orange', 'red', 'brown', 'yellow']
    markers = ['o', '^']

    # SS Time strings for all years and times of day
    datetime_strings = [f"{year}-01-01 {time}" for year in years\
                        for time in times_of_day]

    # SS Getting UTC Time
    times_utc = Time(datetime_strings, format='iso') + 7 * u.hour

    plt.figure(figsize=(10, 6))

    with solar_system_ephemeris.set('builtin'):
        for planet, color in zip(planets, colors):
            planet_positions_lon = np.empty((len(years), len(times_of_day)))
            planet_positions_lat = np.empty((len(years), len(times_of_day)))
            for i, time_utc in enumerate(times_utc):
                planet_position = get_body(planet, time_utc).transform_to(
                    HeliocentricTrueEcliptic(obstime=time_utc))
                year_index = i // len(times_of_day)
                time_index = i % len(times_of_day)
                planet_positions_lon[year_index, time_index] =\
                    planet_position.lon.degree
                planet_positions_lat[year_index, time_index] =\
                    planet_position.lat.degree

            # SS Plotting for each time of day
            for time_index, marker in enumerate(markers):
                plt.scatter(planet_positions_lon[:, time_index],\
                    planet_positions_lat[:, time_index], color=color,\
                    label=f"{planet.capitalize()} {times_of_day[time_index]}",\
                    marker=marker)
    
    # SS Other plotting staffs
    plt.xlabel('Ecliptic Longitude (degrees)')
    plt.ylabel('Ecliptic Latitude (degrees)')
    plt.title('Positions of Non-Earth Planets at 7AM and 7PM MST (2020-2030)')
    plt.legend(loc='upper right')
    plt.grid(True) 
    plt.savefig('planet_position_ecliptic.png') # SS Save plot
    plt.show()


# SS Function to get plot the planets position in equatorial coordinates
# and getting the lowest airmass situation among the plotted planets
def plot_positions_and_lowest_airmass():
    # SS Calling all the planets , years, times of day, and other infos
    planets = ['mercury', 'venus', 'mars', 'jupiter', 'saturn']
    years = np.arange(2020, 2031)
    times_of_day = ['07:00', '19:00']  # SS MST Time
    colors = ['grey', 'orange', 'red', 'brown', 'yellow']
    apo_location = EarthLocation.of_site('Apache Point Observatory')
    lowest_airmass_info = {planet: {'time': None, 'year': None,\
                                    'airmass': np.inf, 'ra': None,\
                                    'dec': None} for planet in planets}
    
    plt.figure(figsize=(12, 8))
    
    # SS Labels for each planet
    planet_labels_added = {planet: False for planet in planets}
    
    for year in years:
        for time in times_of_day:
            # SS Convert to UTC
            adjusted_time = Time(f'{year}-01-01 {time}') + 7 * u.hour
            for planet, color in zip(planets, colors):
                with solar_system_ephemeris.set('builtin'):
                    planet_pos = get_body(planet, adjusted_time,\
                                          location=apo_location)
                    planet_altaz = planet_pos.transform_to(
                        AltAz(obstime=adjusted_time, location=apo_location))
                    if planet_altaz.alt.degree > 0:
                        airmass = planet_altaz.secz
                        if airmass < lowest_airmass_info[planet]['airmass']:
                            display_time = adjusted_time -\
                                7 * u.hour  # Convert back to MST for display
                            lowest_airmass_info[planet] =\
                                {'time': display_time, 'year': year,\
                                 'airmass': airmass, 'ra': planet_pos.ra,\
                                     'dec': planet_pos.dec}
                        if not planet_labels_added[planet]:
                            plt.scatter(planet_pos.ra.hour,
                                        planet_pos.dec.degree, color=color,\
                                        label=planet.capitalize(), marker='o')
                            planet_labels_added[planet] = True
                        else:
                            plt.scatter(planet_pos.ra.hour,\
                                planet_pos.dec.degree, color=color, marker='o')

    # SS Adding lowest airmass label
    lowest_airmass_label_added = False
    for planet, info in lowest_airmass_info.items():
        if info['time']:
            if not lowest_airmass_label_added:
                plt.scatter(info['ra'].hour, info['dec'].degree,\
                            edgecolors='black', facecolors='none', s=100,\
                            marker='o', label="Lowest airmass")
                lowest_airmass_label_added = True
            else:
                plt.scatter(info['ra'].hour, info['dec'].degree,\
                            edgecolors='black', facecolors='none',\
                            s=100, marker='o')

    for planet, info in lowest_airmass_info.items():
        if info['time']:
           print(f"{planet.capitalize()}: Lowest airmass observed at "
      f"{info['time'].iso} MST with airmass "
      f"{info['airmass']:.2f}")


    plt.xlabel('Right Ascension (hours)')
    plt.ylabel('Declination (degrees)')
    plt.title('Equatorial Positions and Lowest Airmass Observations '
          'of Non-Earth Planets (2020-2030)')

    plt.legend()
    plt.grid(True)
    plt.savefig('planet_position_equator.png') # SS Save plot
    plt.show()


# SS Function to get read the quasar file
def get_quasar_file(file_path):
    quasars = []
    with open(file_path, 'r') as file:
        for line in file:
            # SS Getting RA and Dec
            ra_hms = line[:2] + "h" + line[2:4] + "m" + line[4:8] + "s"
            dec_dms = line[9:12] + "d" + line[12:14] + "m" +\
                line[15:].strip() + "s"
            # SS Creating SkyCoords for each quasar
            quasar_coord = SkyCoord(ra_hms, dec_dms, frame='icrs')
            quasars.append(quasar_coord)
    return quasars


# SS Function to get the lowest airmass quasar
def find_lowest_airmass_quasars(month, year=2020):
    quasars = get_quasar_file('HW1quasarfile.dat')
    apo_location = EarthLocation.of_site('Apache Point Observatory')
    
    # SS Number of days in the month
    days_in_month = [31, 29 if year % 4 == 0 else 28, 31, 30,\
                     31, 30, 31, 31, 30, 31, 30, 31]
    
    # SS 11 PM MST (+7 hours to convert to UTC) on the 1st of the month
    start_time_utc = Time(f"{year}-{month:02d}-01 11:00:00") + 7 * u.hour

    lowest_airmass = np.inf
    best_quasar = None
    best_time = None
    
    for day_offset in range(days_in_month[month-1]):
        single_time = start_time_utc + day_offset * u.day
        for quasar in quasars:
            quasar_altaz = quasar.transform_to(AltAz(obstime=single_time,\
                                                     location=apo_location))
            if quasar_altaz.secz < lowest_airmass\
                and quasar_altaz.alt > 0 * u.deg:
                lowest_airmass = quasar_altaz.secz
                best_quasar = quasar
                best_time = single_time
    
    # SS Printing the best quasar
    if best_quasar:
        print(f"Best quasar to observe at lowest airmass:"
             f" {best_quasar.to_string('hmsdms')}")
        print(f"Time (UTC): {best_time.iso}")
        print(f"Airmass: {lowest_airmass}")


def main():
    plot_planet_positions_minimal_loops()
    plot_positions_and_lowest_airmass()
    find_lowest_airmass_quasars(1, 2020)
    
    



# MAIN
###############################
###############################
if __name__ == '__main__':
    main()