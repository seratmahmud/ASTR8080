# S. Saad
# v1 8/26/2020
# v2 1/24/2022
# ASTR 8080 example for HW0

# to import this from another directory:
# import sys
# sys.path.insert(0, '../week1/')
# import polycalc
# from polycalc import get_poly_o3

# IMPORT BLOCK
###############################
###############################
from astropy.coordinates import SkyCoord
from astropy.coordinates import AltAz
from astropy.coordinates import EarthLocation
import astropy.units as u
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt


# FUNCTIONS
###############################
###############################


def print_ra_dec():
    dec_dms = (20, 30, 10)
    ra_hms = (10, 20, 30)
    coord = SkyCoord(ra=ra_hms, dec=dec_dms, frame='icrs', unit=(u.hourangle, u.deg))

    print(f"RA in deg: {coord.ra.deg}")
    print(f"RA in hms: {coord.ra.hms}")
    print(f"RA in hours: {coord.ra.hour}")

    print(f"DEC in deg: {coord.dec.deg}")
    print(f"DEC in hms: {coord.dec.hms}")
    print(f"DEC in hours: {coord.dec.hour}")


def check_jd_mjd_difference():
    jd_now = Time.now().jd
    mjd_now = Time.now().mjd
    jd_mjd_diff = jd_now - mjd_now
    print("The difference between jd and mjd is", jd_mjd_diff)


def print_near_mjds():
    mjd_now = Time.now().mjd
    near_mjd = np.arange(-7, 7, 1.5) + mjd_now
    print("Some near mjds are:", near_mjd)


def print_airmass_at_times():
    APO = EarthLocation(lat="32d46m49.30s", lon="-105d49m13.5s", height=2788*u.m)
    t_dec_dms = (28, 56, 0)
    t_ra_hms = (5, 46, 0)
    target_coord = SkyCoord(ra=t_ra_hms, dec=t_dec_dms, frame='icrs', unit=(u.hourangle, u.deg))

    time_one = Time('2024-01-18 1:00', scale='utc') # 11PM APO Time
    time_two = Time('2023-10-18 8:00') # 2AM APO Time 

    aa_one = AltAz(location=APO, obstime=time_one)
    aa_two = AltAz(location=APO, obstime=time_two)

    aa_one = target_coord.transform_to(aa_one)
    aa_two = target_coord.transform_to(aa_two)

    print("Air mass at 2024-01-18 11PM at APO:", aa_one.secz)
    print("Air mass at 2023-10-18 5PM at APO:", aa_two.secz)


def plot_airmass_vs_time():
    APO = EarthLocation(lat="32d46m49.30s", lon="-105d49m13.5s", height=2788*u.m)
    t_dec_dms = (28, 56, 0)
    t_ra_hms = (5, 46, 0)
    target_coord = SkyCoord(ra=t_ra_hms, dec=t_dec_dms, frame='icrs', unit=(u.hourangle, u.deg))

    time_one = Time('2024-01-18 1:00', scale='utc') # 11PM APO Time
    time_two = Time('2023-10-18 8:00') # 2AM APO Time 

    delta_hours = np.linspace(0, 6, 100)*u.hour # For utc

    full_night_frames_one = AltAz(location=APO, obstime = time_one + delta_hours)
    full_night_frames_one = target_coord.transform_to(full_night_frames_one)
    full_night_frames_two = AltAz(location=APO, obstime = time_two + delta_hours)
    full_night_frames_two = target_coord.transform_to(full_night_frames_two)

    plt.plot(delta_hours, full_night_frames_one.secz)
    plt.gca().invert_yaxis()
    plt.xlabel("Time (hours) from utc 11 PM")
    plt.ylabel("Airmass")
    plt.show()

    plt.plot(delta_hours, full_night_frames_two.secz)
    plt.gca().invert_yaxis()
    plt.xlabel("Time (hours) from utc 6 AM")
    plt.ylabel("Airmass")
    plt.show()


def main():
    print_ra_dec()
    check_jd_mjd_difference()
    print_near_mjds()
    print_airmass_at_times()
     

# MAIN
###############################
###############################
if __name__ == '__main__':
    main()