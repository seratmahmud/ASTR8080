# S. Saad
# v1 1/19/2024
# ASTR 8080 HW0

# IMPORT BLOCK
###############################
###############################
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import time

# FUNCTIONS
###############################
###############################
def generate_line_data(m, b):
    # SS Generating x and y data
    x = np.random.random(10) * 10
    y = m * x + b

    # SS Adding random Gaussian noise
    y_offsets = np.random.normal(loc=0, scale=0.5, size=x.shape)
    y_with_offset = y + y_offsets

    # SS Adding the uncertainty
    y_err = np.full(x.shape, 0.5)

    return x, y_with_offset, y_err

def fit_line(x, y, y_err):
    # SS Fitting a line
    w = 1 / y_err**2
    m2, b2 = np.polyfit(x, y, 1, w=w)

    return m2, b2

def plot_data(x, y, y_err, m, b, m2, b2):
    # SS Setting the plot style for better plots
    mpl.rcParams['lines.linewidth'] = 2
    mpl.rcParams['axes.labelsize'] = 14
    mpl.rcParams['axes.titlesize'] = 16
    mpl.rcParams['xtick.labelsize'] = 12
    mpl.rcParams['ytick.labelsize'] = 12
    mpl.rcParams['legend.fontsize'] = 12

    plt.errorbar(x, y, yerr=y_err, fmt='o', color='blue', ecolor='lightgray', elinewidth=2, capsize=5, label='Data with uncertainty')
    
    # SS Plotting the original line
    x_line = np.linspace(min(x), max(x), 100)
    y_line = m * x_line + b
    plt.plot(x_line, y_line, color='green', linestyle='-', label='Original line (y=mx+b)')

    # SS Plotting the best-fitting line
    y_fit_line = m2 * x_line + b2
    plt.plot(x_line, y_fit_line, color='red', linestyle='--', label='Fitted line (y=m2x+b2)')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Line Fit Plot')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.show()
    plt.savefig('line_fit_plot.png', format='png', dpi=300)


# MAIN
###############################
###############################
if __name__ == '__main__':
    start_time = time.time()

    # SS Checking the command line arguments for m and b
    if len(sys.argv) != 3:
        print("Usage: python <script_name.py> <m> <b>")
        sys.exit(1)

    m = float(sys.argv[1])
    b = float(sys.argv[2])

    # SS Generating data
    x, y_with_offset, y_err = generate_line_data(m, b)
    print(f"Data generated. Time taken: {time.time() - start_time:.2f} seconds")

    # SS Fitting line
    start_time = time.time()
    m2, b2 = fit_line(x, y_with_offset, y_err)
    print(f"Line fitted. Time taken: {time.time() - start_time:.2f} seconds")

    # SS Plotting data
    start_time = time.time()
    plot_data(x, y_with_offset, y_err, m, b, m2, b2)
    print(f"Data plotted. Time taken: {time.time() - start_time:.2f} seconds")
