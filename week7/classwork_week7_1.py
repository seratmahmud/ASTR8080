# S. Saad
# ASTR 8080 week7, classwork1

# IMPORT BLOCK
###############################
###############################


from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from io import StringIO
import mechanize
import pandas as pd



# FUNCTIONS
###############################
###############################


def SDSS_select(sql):
    url = "https://skyserver.sdss.org/dr16/en/tools/search/sql.aspx"
    br = mechanize.Browser()
    br.open(url)
    br.select_form(name="sql")
    br['cmd'] = sql
    #br['format']=['fits']
    br['format']=['csv']
    response = br.submit()
    print(response.get_data())
    file_like = StringIO(str(response.get_data()))
    
    
    return pd.read_csv(file_like, skiprows=1)




def main():
    sql_q = "SELECT p.ObjID, p.ra, p.dec, p.g FROM photoObj p JOIN dbo.fGetNearbyObjEq(300,-1,2) n ON n.objID = p.objID"
    

    
    t = SDSS_select(sql_q)
    
    
    print(t)
    
    
    
    

    
    ra = t['ra']
    dec = t['dec']
    
    g = t['g']
    print(g)
    a = np.where(g < 16)[0]
    b = np.where((g < 19) & (g >= 16))
    c = np.where(g>=19)
    
    s = np.zeros(len(g))
    
    s[a] = 10
    s[b] = 5
    s[c] = 1
    
    

    plt.figure(figsize=(10,6))
    plt.scatter(ra, dec, s = s)
    plt.xlabel('Right Ascension (degrees)')
    plt.xlabel('Declination (degrees)')
    plt.show()








# MAIN
###############################
###############################
if __name__ == '__main__':
    main()








