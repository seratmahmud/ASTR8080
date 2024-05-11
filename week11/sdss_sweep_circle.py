import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from time import time
#import fitsio
from astropy.io import fits
from sdss_sweep_data_index import sdss_sweep_data_index

def sdss_sweep_circle(ra, dec, radius, 
                      objtype='star', sweepdir='./', all=False, verbose=False):
    """
    NAME: sdss_sweep_circle

    PURPOSE: Return all of the objects in the SDSS sweeps files within a
      circular (angular) radius of a SINGLE ra, dec position

    CALLING SEQUENCE: from the Python prompt
   
      from sdss_sweep_circle import sdss_sweep_circle
      objs = sdss_sweep_circle(ra, dec, radius, [, objtype=, sweepdir=, all=, verbose=])

    INPUTS:

      ra, dec - central location (J2000 deg)
      radius - search radius (deg)

    OPTIONAL INPUTS:

      objtype - type of object to search for, from 'star', 'gal', 'sky'
          [default 'star']
      sweepdir - the location of the directory that contains (the 301 
          directory that contains) sweeps files from an SDSS Data Release
      all - send as True to keep all objects, not just SURVEY_PRIMARY
      verbose - send as True to print timing messages

    OUTPUTS:

      objs - the list of objects in the SDSS sweep files that are within the passed
          radius of the passed ra, dec position

    COMMENTS: ra, dec and radius must all be floats
              will take ~2 minutes to run for a ~3o radius for objtype='star'
              and about ~6 minutes to run for a ~3o radius for objtype='gal'

    REQUIREMENTS: sdss_sweep_data_index.py must be in a location that can be imported
    
    REVISION HISTORY:
    
    v1.1: version Jessie C. Runnoe, October 16, 2019 updated to use
          astropy.io.fits instead of fitsio for ASTR8020 data mining class
    v1.0: version Adam D. Myers, March 25, 2017 as an example for
          ASTR5160 data mining class
    """
  
    #ADM start timer
    t0 = time()

    #ADM determine the sweeps files of interest
    swfiles = sdss_sweep_data_index(ra, dec, radius, objtype=objtype,
                                    sweepdir=sweepdir, all=all, verbose=False)
    if verbose:
        print('Determined sweeps files of interest...t = {:.1f}s'.format(time()-t0))

    #ADM read in all of the objects of interest
    #JCR updated for astropy.io.fits instead of fitsio
    #objs = [ fitsio.read(file) for file in swfiles ]
    objs_struct = [ fits.open(file) for file in swfiles ]
    objs        = [obj[1].data for obj in objs_struct ] 


    if verbose:
        print('Read in sweeps files...t = {:.1f}s'.format(time()-t0))

    #ADM reform the list of objects into one long array
    objs = np.hstack(objs)

    if verbose:
        print('Formatted sweeps files...t = {:.1f}s'.format(time()-t0))

    #ADM unless we want all objects, restrict to just PRIMARY objects
    primaryflag = 2**8
    if not all:
        w = np.where( (objs["RESOLVE_STATUS"] & primaryflag) != 0)
        #ADM if there are no primary objects, tell the user
        if len(w[0]) == 0:
            raise IOError("There are no primary objects in the SDSS within {} degrees of (RA,Dec)=({},{})".format(radius, ra, dec))
        objs = objs[w]

    #ADM set up skycoord objects
    cin = SkyCoord(ra=ra*u.degree, dec=dec*u.degree) 
    csweeps = SkyCoord(ra=objs["RA"]*u.degree, dec=objs["DEC"]*u.degree)

    #ADM find coordinate matches between objects in the sweeps and
    #ADM the input RA, Dec, using the input matching radius
    sep = cin.separation(csweeps)
    m2 = np.where(sep < radius*u.degree)[0]

    #ADM return the objects that were a coordinate match
    if verbose:
        print('Done...t = {:.1f}s'.format(time()-t0))

    return objs[m2]
