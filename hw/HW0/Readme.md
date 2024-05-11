# Line Fitting Module

## Descriptions
The script generates a set of data points based on a linear equation with a specified slope (m) and y-intercept (b), adds Gaussian noise to these points, and then fits a line to the data. The output includes a plot showing the generated data points, the original line, and the best fitted line. The plot is also saved as a PNG file. The script also prints the time taken for the whole process.

## Requirements
- Python 3.x
- Numpy
- Matplotlib

## Usages
To run this script, you need to pass two command-line arguments: the slope (m) and y-intercept (b) of the original line. The script is executed as follows:

	python line_fit_module.py <m> <b>

## Examples
For example, to run the script with a slope of 2 and a y-intercept of 1, use the following:
	python line_fit_module.py 2 1

## Concerns
You need to ensure that there's no file with the same name (line_fit_plot.png) in the directory.


## Author
S. Saad
