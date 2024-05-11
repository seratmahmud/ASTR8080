# J. Runnoe
# 04/11/2022
# ASTR 8080 tasks for week14 2022
import numpy as np
import pdb
import csv
import pandas as pd
import os

# for SDSS SQL queries
import mechanize
from io import BytesIO

def get_coord_string(ras,decs):
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   Take inputs of RA and DEC and return a search string for crossID search.
#
# CALLING SEQUENCE:
#   coord_str = get_coord_string(ras,decs)  
#
# INPUTS:
#   ras    - an array of RAs
#   decs   - an array of decs of the same length
#
# OUTPUTS:
#   coord_str     - a single string formatted for the crossID website SQL search box 
#
# NOTES:
#
# Based on code here: http://balbuceosastropy.blogspot.com/2013/10/an-easy-way-to-make-sql-queries-from.html
#
# MODIFICATION HISTORY:
#	 JCR, 15 April 2022: VERSION 1.00
#           -  Started writing the first version of the code. 
#-
#------------------------------------------------------------- 
    coord_str = "  name  ra       dec\r\n"
    for i,ra,dec in zip(range(len(ras)),ras,decs): coord_str+="  {}   {}   {}\r\n".format(i,ra,dec)
    return coord_str

def sdss_crossid_sql_query(sql,ra,dec,FROM_FILE=True):
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   Take an input string and run an SDSS SQL query over the internet.
#   This is limited to 150 entries.
#
# CALLING SEQUENCE:
#   result = sdss_sql_query(sql,FROM_FILE=True)   
#
# INPUTS:
#   sql       - a string containing an SQL query.
#   ra,dec    - arrays of N RAs and DECS
#   FROM_FILE - defuault True saves coords to a file and uploads to query, good for lists up to 10000
#               if False, then it is a string of coordinates formatted by get_coord_string() which can fail for shorter lists
#
# OUTPUTS:
#   result - a recarray with the query results
#
# NOTES:
#
# Based on code here: http://balbuceosastropy.blogspot.com/2013/10/an-easy-way-to-make-sql-queries-from.html
#
# Modified for CrossID website. Need to determine the form name and field name, done at the above URL.
#   Use the following to figure out form name:
#        url = "http://skyserver.sdss.org/dr16/en/tools/crossid/crossid.aspx"
#        br  = mechanize.Browser()
#        br.open(url)
#        for f in br.forms(): print(f.name) #--> returns only crossid 
#        for f in br.forms(): [print('\t',c.name, '\t', c.type) for c in f.controls] #-->uquery is the second text box 
#
# MODIFICATION HISTORY:
#	 JCR, 7 December 2020: VERSION 1.00
#           -  Started writing the first version of the code. 
#-
#------------------------------------------------------------- 
    # write to a file or create a string to upload to webform text box
    if FROM_FILE:
        with open('sql_coord_input.tab', 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["NUM","RA","        DEC"])
            writer.writerows(zip(np.arange(len(ra)),ra,dec))
    else:
        coord_str  = get_coord_string(ra,dec)

    # set up the webform
    url = "http://skyserver.sdss.org/dr16/en/tools/crossid/crossid.aspx"
    br  = mechanize.Browser()
    br.open(url)
    br.select_form(name="crossid")
    br['uquery'] = sql              # enter your query
    br['format'] = ['csv']          # return in csv format
    br['searchType'] = ['spectro']  # query specobjall
    if FROM_FILE:
        br.form.add_file(open('sql_coord_input.tab'),'text/plain','sql_coord_input.tab') # read coordinates from file
    else:
        br['paste'] = coord_str     # enter coordinates into webform text box
    response = br.submit()          # submit form
    file_like = BytesIO(response.get_data())
    
    # clean up
    if FROM_FILE: os.system('rm sql_coord_input.tab')

    return pd.read_csv(file_like,  skiprows=1).to_records()

def gen_wget(mjds,plates,fibers,run2ds):
#-------------------------------------------------------------
#+
# PURPOSE:                                                 
#   Takes input MJD PLATE FIBER and creates wget commands.
#
# CALLING SEQUENCE:
#   com = gen_wget(mjd,plate,fiber,SCRIPT=True)   
#
# INPUTS:
#   mjd      - An array of MJDs.
#   plate    - An array of plate numbers.
#   fiber    - An array of fiber ids.
#
# OUTPUTS:
#   com    - An array of string with wget commands in them. 
#   files  - An array of the filenames.
#
# NOTES:
#
#
# MODIFICATION HISTORY:
#	 JCR, 15 April 2022: VERSION 1.00
#           -  Started writing the first version of the code. 
#-
#------------------------------------------------------------- 
    # data are stored in three different places depending on the spectrograph and run2d 
    com   = [] 
    files = [] 
    spectrograph = 'eboss'
    for mjd,plate,fiber,run2d in zip(mjds,plates,fibers,run2ds): 
        if run2d=='v5_13_0':
            spectrograph = 'eboss'
        else:
            spectrograph = 'sdss'
        com.append("wget https://data.sdss.org/sas/dr16/{0:}/spectro/redux/{1:}/spectra/lite/{2:}/spec-{2:}-{3:}-{4:04d}.fits".format(spectrograph,run2d,plate,mjd,fiber))
        files.append("spec-{}-{}-{:04d}.fits".format(plate,mjd,fiber))

    return np.array(com),np.array(files)

if __name__=="__main__":

    # JCR run a search for 3 test stars
    ra = [181.57977,182.60141,182.74392]       
    dec = [30.61939,30.33945,30.27135]

    # JCR write an SQL query
    sql_query = 'SELECT s.ra, s.dec,s.mjd,s.plate,s.fiberid,s.run2d,p.psfMag_g\
    FROM #upload u\
          JOIN #x x ON x.up_id = u.up_id\
          JOIN SpecObjAll s ON s.specObjID = x.specObjID\
          JOIN PhotoTag p ON p.objID = s.bestObjID'
    
    # JCR run the SQL query and download data 
    sql_data = sdss_crossid_sql_query(sql_query,ra,dec,FROM_FILE=True)


    # JCR create a wget file from the return
    #     note that you need to query for RUN2D, MJD, PLATE, FIBER
    wget_coms,spec_files = gen_wget(sql_data['mjd'],sql_data['plate'],sql_data['fiberid'],sql_data['run2d'])

    # JCR you can run the commands with the following
    #     the spectra will download into your current dir
    #for com in wget_coms:
    #    os.system(com)

    import pdb;pdb.set_trace()
