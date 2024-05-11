from stripe82query import sdssQuery
import matplotlib.pyplot as plt
import numpy as np

def stripe(ra, dec):
    qry = sdssQuery()
    qry.query = f"""SELECT p.ra, p.dec, p.psfmag_g, f.mjd_g, n.distance*60 FROM fGetNearbyObjEq({ra},{dec},0.03) n, PhotoPrimary p, Field f WHERE n.objID = p.objID and f.fieldID = p.fieldID"""

    mjd_g_values = []
    psfmag_g_values = []
    
    for line in qry.executeQuery():
        processed_line = line.strip().decode('utf-8')
        try: 
            if len(processed_line) >= 4:
                processed_line = processed_line.split(",")
                mjd_g, psfmag_g = float(processed_line[3]), float(processed_line[2])
                mjd_g_values.append(mjd_g)
                psfmag_g_values.append(psfmag_g)
        except:
            continue
    a = np.where(np.array(mjd_g_values) != 0)[0]
    mjd_g_values = np.array(mjd_g_values)[a]
    psfmag_g_values = np.array(psfmag_g_values)[a]
    
    plt.figure(figsize=(10,6))
    plt.scatter(mjd_g_values, psfmag_g_values)
    plt.xlim(min(mjd_g_values)-100, max(mjd_g_values)+100)
    plt.ylim(max(psfmag_g_values) + 0.2, min(psfmag_g_values) - 0.2)
    plt.ylabel("g")
    plt.xlabel("mjd")
    plt.show()
    plt.savefig(f"{ra}_{dec}.png")
    
def main():
    ra = [29.2256832, 35.3756676, 45.299833, 58.175468, 60.829041]
    dec = [0.4208970,  0.0017000, -0.55386111,  0.218697, -1.240793]
    
    for i in range(len(ra)):
        ra_i = ra[i]
        dec_i = dec[i]
        stripe(ra_i, dec_i)

def test():
    ra = 29.2256832
    dec = 0.4208970
    qry = sdssQuery()
    qry.query = f"""SELECT p.ra, p.dec, p.psfmag_g, f.mjd_g, n.distance*60 FROM fGetNearbyObjEq({ra},{dec},0.03) n, PhotoPrimary p, Field f WHERE n.objID = p.objID and f.fieldID = p.fieldID"""

    mjd_g_values = []
    psfmag_g_values = []
    
    for line in qry.executeQuery():
        processed_line = line.strip().decode('utf-8')
        try: 
            if len(processed_line) >= 4:
                processed_line = processed_line.split(",")
                mjd_g, psfmag_g = float(processed_line[3]), float(processed_line[2])
                mjd_g_values.append(mjd_g)
                psfmag_g_values.append(psfmag_g)
        except:
            continue
    # Removing the noise

    plt.figure(figsize=(10,6))
    plt.scatter(mjd_g_values, psfmag_g_values)
    plt.xlim(50000, )
    plt.ylim(17.8, 17)
    plt.ylabel("g")
    plt.xlabel("mjd")
    plt.show()
    plt.savefig("{ra}_{dec}.png")


if __name__ == '__main__':
    main()




















