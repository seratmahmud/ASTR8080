ASTR 8080 Final (See Comments at the end)

## Author: S. Saad 

### Overview
1. This Python script will find the SDSS stars and galaxies in a certian area of a sky.
2. It will then draw mangle mask for the stars less than a certain magnitude for that area in the sky. Then it will return the galaxies that will be with 5'' radius of each of these stars.
3. It will then calculate the density of galaxies in that area both with and without considering the masked region.
4. It will also plot an aitoff projection of the selected area along with the mask

### Requirements
- Python
- Astropy
- Mangle
- Matplotlib
- Numpy

### Setup and Execution
1. Ensure you have Python and the required libraries installed.
2. Ensure you have the `.py` script and `HW3quasarfile.dat` file to your local machine.
3. Open a terminal and execute the file. Write 'python hw.py'.

### Sample Output:

    (base) [saadsm@tomservo final]$ python final.py
    Total area: 0.3833009065188101 square radians
    Total mask area: 3.833009065188101e-06 square radians
    Total number of galaxies in the mask: 1127
    Number density total: 2509631.4243984018 galaxies per square radians
    Number density with mask: 2506716.2428402216 galaxies per square radians
    Total time taken to run the file: 148.61586689949036 seconds.


### Comments
1. I know that the time it's taking to run is a lot. I think it's taking time when I'm using mask.contains method to find out the sources that are in the mask. This method works very smoothly when the total number of polygons in the `.ply` file is not a lot. But for my case, I had a total of 7000+ polygons. So it's taking a lot of time. Not sure if there's any other way to find out what sources are there inside a mangle mask.

2. I also calculated the uncertainty in the mask area calculation using the random number. This is calculating the area in 10 iterations and then claculating the uncertainty. This is not always necessary, specially for the cases when we are using a lrge number of random points to claculate the area. So I commented out that part of the code as it takes a lot of time to run and calculate the uncertainty.

3. I saved the two files at first - one of them is the fits file for all the stars in that rectangluar segment that has r_magnitude less than 10, another one is all the galaxies that are in that rectangular segment that has r_corrected_magnitude less than 19. If we generate these files everytime, it takes a few minutes to run.

4. I created the `sweep_data_rectangular` function that will basically do the same thing as `sdss_sweep_data_index.py` file. But it will get sweep objects for a rectangular region rather than a circular region.

5. Definitely one of the best classes I have ever taken. Learned a lot. So grateful!