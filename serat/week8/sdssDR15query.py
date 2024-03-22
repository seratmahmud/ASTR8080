class sdssQuery:
    """
    NAME: sdssQuery
 
    PURPOSE: class that can be initialized using SoapPy to 
    send an SQL command to SDSS web services
 
    CALLING SEQUENCE: from the UNIX command line:
      
      python sdssDR15query.py ra dec

    INPUTS: ra and dec shoud be sent from the command line:

      ra - Right Ascension of position to query around in SDSS DR15
      dec - declination of position to query around in SDSS DR15

    OPTIONAL INPUTS: 

    OUTPUTS: the result of the SQL command called "query" in the
    code, below, is executed by the SDSS DR15 SQL web services and
    printed out at the command line
    
    COMMENTS: This is hacked from an example provide by the SDSS 
    in an early tutorial on web services

    Note that the SQL command passed as "query" can be changed to
    any valid SDSS string
 
    EXAMPLES: At the Unix command line:

      python sdssDR15query.py 145.285854 34.741254
      
    should return:

      45.28585385223 34.7412541787531 21.12497 20.06535 19.61032 19.37788 19.42766 0.000777943878117123
    
    REVISION HISTORY: 

    v2.0: This DR15 python3 version Jessie Runnoe, October 6, 2019 as an example for
          ASTR8020 data mining class

    v1.0: This DR9 version Adam D. Myers, March 4, 2017 as an example for
          ASTR5160 data mining class
    """

    url='http://skyserver.sdss.org/dr15/en/tools/search/x_sql.aspx'
    format = 'csv'
    def __init__(self):
        self.query = ''
        self.cleanQuery = ''
    def executeQuery(self):
        from urllib.parse import urlencode
        from urllib.request import urlopen
        self.filterQuery()
        params = urlencode({'cmd': self.cleanQuery, 'format':self.format})
        return urlopen(self.url + '?%s' % params)
    def filterQuery(self):
        from os import linesep
        self.cleanQuery = ''
        tempQuery = str.lstrip(self.query)
        for line in tempQuery.split('\n'):
            self.cleanQuery += line.split('--')[0] + ' ' + linesep;
if __name__ == '__main__':
    import sys, string
    from time import process_time, sleep

    ra1,dec1 = sys.argv[1],sys.argv[2]

    s = process_time()
    out = sys.stdout
    qry = sdssQuery()

    query = """SELECT top 1 ra,dec,u,g,r,i,z,GNOE.distance*60 FROM PhotoObj as PT JOIN dbo.fGetNearbyObjEq("""+str(ra1)+""","""+str(dec1)+""",0.03) as GNOE on PT.objID = GNOE.objID ORDER BY GNOE.distance"""

    qry.query = query
    for line in qry.executeQuery():
        a = line.strip()

    sleep(1)
    
    list = a.split(b",")

    print((b" ".join(list)).decode('utf-8'))

