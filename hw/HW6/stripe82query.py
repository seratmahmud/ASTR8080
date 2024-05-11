class sdssQuery:
    """
    NAME: sdssQuery
 
    PURPOSE: class that can be initialized using SoapPy to 
    send an SQL command to SDSS web services
 
    CALLING SEQUENCE: from the UNIX command line:
      
      python stripe82query.py ra dec

    INPUTS: ra and dec shoud be sent from the command line:

      ra - Right Ascension of object to match to SDSS
      dec - declination of object to match to SDSS

    OPTIONAL INPUTS: 

    OUTPUTS: the result of the SQL command called "query" in the
    code, below, is executed by the SDSS DR7 SQL web services and
    printed out at the command line
    
    COMMENTS: This is hacked from an example provide by the SDSS 
    in an early tutorial on web services

    Note that the SQL command passed as "query" can be changed to
    any valid SDSS string

 
    EXAMPLES: At the Unix command line:

    python stripe82query.py 29.2256832 0.4208970
      
    should return
      
    ra dec psfMag_u psfMag_g psfMag_r psfMag_i psfMag_z psfMagErr_u psfMagErr_g psfMagErr_r psfMagErr_i psfMagErr_z mjd_u mjd_g mjd_r mjd_i mjd_z
    29.22568203 0.42089991 18.536137 17.307686 16.916077 16.779457 16.728699 0.021424 0.019576 0.015824 0.014918 0.016103 5.29094022E4 5.29094039E4 5.29094006E4 5.29094014E4 5.29094031E4
    29.2256856 0.42089434 18.612911 17.393892 16.930103 16.769346 16.72949 0.022699 0.041158 0.016824 0.014189 0.017772 5.40083443E4 5.4008346E4 5.40083427E4 5.40083435E4 5.40083452E4
    29.22568311 0.42089287 18.546591 17.343668 16.944519 16.775738 16.706455 0.020216 0.013211 0.013618 0.014646 0.018711 5.22252604E4 5.22252621E4 5.22252587E4 5.22252596E4 5.22252612E4
    29.2256852 0.42090197 18.535133 17.335161 16.930634 16.757013 16.747959 0.020706 0.017592 0.01382 0.014552 0.022086 5.25772885E4 5.25772902E4 5.25772869E4 5.25772877E4 5.25772893E4
    29.22568055 0.42089266 18.625237 17.44462 16.946463 16.820351 16.760046 0.035377 0.062818 0.021703 0.018584 0.019226 5.37052929E4 5.37052945E4 5.37052912E4 5.3705292E4 5.37052937E4
    ...etc...

    REVISION HISTORY: 
    v3.0: This version Serat Saad, April 20, 2024 updated
        for the HW6 of ASTR 8008 class spring 2024, Vanderbilt Univ.
        To add more columns.

    v2.0: This version Jessie C. Runnoe, November 4, 2019 
          updated for Python 3 

    v1.0: This version Adam D. Myers, April 8, 2015 as an example for
          ASTR5160 data mining class
    """

    url = 'http://cas.sdss.org/stripe82/en/tools/search/x_sql.asp'
    format = 'csv'
    def __init__(self):
        self.query = ''
        self.cleanQuery = ''
    def executeQuery(self):
        #JCR changed for python3
        from urllib.parse import urlencode
        from urllib.request import urlopen
        self.filterQuery()
        params = urlencode({'cmd': self.cleanQuery, 'format':self.format})
        return urlopen(self.url + '?%s' % params)
    def filterQuery(self):
        from os import linesep
        self.cleanQuery = ''
        #JCR changed for python3
        tempQuery = str.lstrip(self.query)
        for line in tempQuery.split('\n'):
            self.cleanQuery += line.split('--')[0] + ' ' + linesep;
if __name__ == '__main__':
    import sys, string
    from time import process_time, sleep

    ra1,dec1 = sys.argv[1],sys.argv[2]

    #JCR changed for python3
    s = process_time()
    out = sys.stdout
    qry = sdssQuery()

    # SS changed query for getting more columns
    query = """SELECT p.ra, p.dec, p.psfMag_u, p.psfMag_g, p.psfMag_r, p.psfMag_i, p.psfMag_z, f.mjd_u, f.mjd_g, f.mjd_r, f.mjd_i, f.mjd_z FROM fGetNearbyObjEq("""+str(ra1)+""","""+str(dec1)+""",0.03) n, PhotoPrimary p, Field f WHERE n.objID = p.objID and f.fieldID = p.fieldID"""
    
    print(query)

    qry.query = query
    for line in qry.executeQuery():
        #JCR changed for python3
        a = line.strip()

        list = a.split(b",")

        print((b" ".join(list)).decode('utf-8'))

    sleep(1)
   
