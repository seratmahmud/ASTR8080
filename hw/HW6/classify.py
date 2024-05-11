import numpy as np
import os


def classify_objects(ra, dec):
    """
    NAME: classify_objects

    PURPOSE: Returning the type of source

    CALLING SEQUENCE: from the Python prompt
   
      from classify import classify_objects
      type = classify_objects(ra, dec)

    INPUTS:

      ra and dec arrays for positions in the sky

    OUTPUTS:

      array - an array that has a total length of data_array, the elements
      in that array are as follows:
      0: supernova
      1: quasar
      2: variable star
      3: standard star

    REVISION HISTORY:

    v1.0: version Serat Saad, April 22, 2024 to submit as a homework for
    ASTR 8080 class at Vanderbilt University
    """
  
    # SS Initializing an empty list to store data
    data_list = []

    # SS Looping over each coordinate pairs
    for r, d in zip(ra, dec):
        # SS Filename for the output data and running the command to query
        data_file_name = f"query_result_{r}_{d}.data"
        command = f"python stripe82query.py {r} {d} > {data_file_name}"
        os.system(command)
        
        # SS Reading data from the file
        with open(data_file_name, 'r') as file:
            data = file.readlines()[2:]  # SS skipping header
            # SS Converting data into numpy array if data is present
            if data: 
                rows = np.array([row.split() for row in data], dtype=float)
                data_list.append(rows)
            else:
                data_list.append(np.array([]))
        # SS Removing the data file for cleaning up
        os.remove(data_file_name)
    
    # SS Defaulting the classification type to 3
    obj_type = np.zeros(len(data_list)) + 3
    
    for i, data in enumerate(data_list):
        if data.size == 0:
            continue  # Skipping if there is no data for that coordinate
        
        # SS Getting psf magnitude and mjd
        psfMag = data[:, 2:7]
        mjd = data[:, 7:12]
             
        # SS Removing invalid mjds
        if not np.any(mjd[:, 4] > 0):
            continue
        
        # SS Using g magnitude for classification, getting mean and distribution of g magnitude
        g_mag = psfMag[mjd[:, 4] > 0, 1]
        g_mean = np.mean(g_mag)
        g_std = np.std(g_mag)
        top_mean = np.mean(np.sort(g_mag)[:3])  # SS Mean of the top 3 magnitudes
        
        # SS Suprnova(0) for significant mean deviation in top magnitudes
        if top_mean/g_mean < 0.95:
            obj_type[i] = 0  
        # SS Variable stars(2) for high standard deviation, which indicates variability
        elif g_std > 0.8:
            obj_type[i] = 2 
        else:
            # SS Calculate corrected magnitudes and color indices
            # SS seems like only u_g ang g_r is enough?
            corrected_mag = np.mean(psfMag[mjd[:, 4] > 0], axis=0)
            u_g = corrected_mag[0] - corrected_mag[1]
            g_r = corrected_mag[1] - corrected_mag[2]
            # SS Assigning type based on color indices ranges from HW5
            # SS QSO(1) for objects within specific color index ranges
            if -0.2 < u_g < 0.6 or -0.1 < g_r < 0.4:
                obj_type[i] = 1  
    return obj_type.astype(int)  # SS Returning the array entries as integers

