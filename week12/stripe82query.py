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
      
    ra dec psfmag_g mjd_g 
    29.22568203 0.42089991 17.307686 5.29094039E4 0.01130151
    29.2256856 0.42089434 17.393892 5.4008346E4 0.01289675
    29.22568311 0.42089287 17.343668 5.22252621E4 0.0148656
    29.2256852 0.42090197 17.335161 5.25772902E4 0.01928333
    29.22568055 0.42089266 17.44462 5.37052945E4 0.01831076
    29.22568726 0.42089877 17.305477 5.40613266E4 0.01595195
    29.2256791 0.42089037 17.320415 5.22880964E4 0.02807115
    29.22567693 0.4208946 17.343426 5.25762359E4 0.02417368
    29.22567698 0.42090044 17.362961 5.2910374E4 0.02557164
    29.22568585 0.42090325 17.426786 5.36403587E4 0.02444315
    ...etc...

    REVISION HISTORY: 
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

    query = """SELECT p.ra, p.dec, p.psfmag_g, f.mjd_g, n.distance*60 FROM fGetNearbyObjEq("""+str(ra1)+""","""+str(dec1)+""",0.03) n, PhotoPrimary p, Field f WHERE n.objID = p.objID and f.fieldID = p.fieldID"""
    
    print(query)

    qry.query = query
    for line in qry.executeQuery():
        #JCR changed for python3
        a = line.strip()

        list = a.split(b",")

        print((b" ".join(list)).decode('utf-8'))

    sleep(1)
   
