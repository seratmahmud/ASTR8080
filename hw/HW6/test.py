import numpy as np
from classify import classify_objects
import time
import matplotlib.pyplot as plt

def test_classify(ra_arr, dec_arr):
    """
    NAME: test_classify

    PURPOSE: Returning the type of source

    CALLING SEQUENCE: from the Python prompt
   
        test_classify(ra_arr, dec_arr)

    INPUTS:

      ra and dec arrays for positions in the sky

    OUTPUTS:

      Doesn't return anything but will print the type of source and will calculate
      total time taken to run the program. This will help to score the homoework.
      
    REVISION HISTORY:

    v1.0: version Serat Saad, April 22, 2024 to submit as a homework for
    ASTR 8080 class at Vanderbilt University
    """
    print("Test")

    # SS Calling the function with time measurement
    start = time.time()
    test_types = classify_objects(ra_arr, dec_arr)
    end = time.time()
    
    # SS Printing the total time taken to run
    print('Time taken to run the function: ', end-start, 's.')
    
    # SS Printing the types of sources
    for i in range(len(test_types)):
        if test_types[i] == 0:
            print(f"{i+1}. Supernova")
        elif test_types[i] == 1:
            print(f"{i+1}. Qso")
        elif test_types[i] == 2:
            print(f"{i+1}. Variable Star")
        elif test_types[i] == 3:
            print(f"{i+1}. Standard Star")
    

if __name__ == '__main__':
    # SS ra and dec arrays, these can be changed, but they always should be arrays
    ra_arr = [29.2256832, 35.3756676, 45.299833, 58.175468, 60.829041]
    dec_arr = [0.4208970, 0.0017000, -0.55386111, 0.218697, -1.240793]
    
    # SS Calling the test_classify function
    test_classify(ra_arr, dec_arr)
