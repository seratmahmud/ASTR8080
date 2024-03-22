import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from time import time
from astropy.io import fits 

def sdss_sweep_data_index(ra, dec, radius, 
                          objtype='star', sweepdir='.', all=False, verbose=False):
    """
    NAME: sdss_sweep_data_index

    PURPOSE: Return the calibObj* (sweeps) objects that a set of 
             ra, dec, radius positions intersect 

    CALLING SEQUENCE: from the Python prompt
   
      from sdss_sweep_data_index import sdss_sweep_data_index
      swfiles = sdss_sweep_data_index(ra, dec, radius, [, objtype=])

    INPUTS:

      ra, dec - central location (J2000 deg)
      radius - search radius (deg)

    OPTIONAL INPUTS:

      objtype - type of object to search for, from 'star', 'gal', 'sky'
          [default 'star']
      sweepdir - the location of the directory that contains (the 301 
          directory that contains) sweeps files from an SDSS Data Release
      all - send as True to keep all objects, not just SURVEY_PRIMARY
      verbose - send as True print timing messages

    OUTPUTS:

      swfiles - the unique list of sweeps files that intersect the passed
          ras, decs and radius

    COMMENTS: ra and dec inputs can be arrays, RADIUS must be a single float

    REVISION HISTORY:
    v2.0: version Jessie Runnoe, October 6, 2019 as an example for
          ASTR8020 data mining class. Edited for astropy.io.fits 
          instead of fitsio

    v1.0: version Adam D. Myers, March 4, 2017 as an example for
          ASTR5160 data mining class
    """
  
    #ADM start timer
    t0 = time()

    #ADM read in index file dependent on passed type of object
    #JCR updated for astropy.io.fits instead of fitsio
    indexfile = sweepdir+'/datasweep-index-'+objtype+'.fits'
    index_struct = fits.open(indexfile)
    index = index_struct[1].data

    #ADM unless we want all objects, don't check sweeps files without primary objs
    if not all:
        w = np.where(index["NPRIMARY"] > 0)
        index = index[w]

    #ADM set up skycoord objects
    cin = SkyCoord(ra=ra*u.degree, dec=dec*u.degree) 
    cindex = SkyCoord(ra=index["RA"]*u.degree, dec=index["DEC"]*u.degree)

    #ADM find matching sweeps files (margin of indexing is 0.36 degress)
    #ADM use search_around_sky for an array, separation for a scalar
    if type(ra) == type(1.) or type(ra) == type(1):
        sep = cin.separation(cindex)
        m2 = np.where(sep < (radius+0.36)*u.deg)[0]
    else:
        m1, m2, d2d, d3d = cindex.search_around_sky(cin, (radius+0.36)*u.deg)

    if len(m2) == 0:
        raise IOError('no objects in SDSS area')

    #ADM find unique runs and camcols
    #JCR updated for python3
    #rc = index[m2]["RUN"]*6L+index[m2]["CAMCOL"]-1
    rc = index[m2]["RUN"]+index[m2]["CAMCOL"]-1
    #ADM uniqindices contains an index to each unique combination
    #ADM of runs and camcols
    uniq, uniqindices = np.unique(rc, return_index=True)
    #ADM the information about the unique sweep files of interest
    swinfo = index[m2][uniqindices]
    
    #JCR remove rerun=157 reduction and return only 301
    latest_red = (swinfo['RERUN'])=='301'

    #ADM now we just need to format the string correctly
    #JCR edited to remove 157 rerun
    swfiles = [ '{}/{}/calibObj-{:06}-{}-{}.fits.gz'.format(sweepdir,sw["RERUN"],sw["RUN"],sw["CAMCOL"],objtype) for sw in swinfo[latest_red] ]

    if verbose:
        print('Done...t = {}s'.format(time()-t0))

    return swfiles
