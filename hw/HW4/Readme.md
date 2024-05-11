ASTR 8080 Homework 4

## Author: S. Saad

### Overview
This python script does the following:

1. Identifies FIRST radio survey sources within specified regions of the sky.
2. Matches these FIRST sources with SDSS primary point sources and extracts PSFFLUXES.
3. Retrieves WISE W1 and W2 magnitudes for these matched sources.
4. Determines the brightest FIRST source in the WISE W1 band within the survey region and plots its 7 fluxes (5 SDSS + 2 WISE) as a function of wavelength.
5. Retrieves and displays the SDSS optical image and spectrum for the brightest source.

### Requirements
- Python
- Astropy
- Matplotlib
- Numpy

### Setup
1. Ensure you have Python and the required libraries installed.
2. Ensure you have the `.py` script and `sdss_sweep_data_index.py` file to the same directory of your local machine.

## Usages
The script is executed as follows in the terminal:

	python hw4.py
