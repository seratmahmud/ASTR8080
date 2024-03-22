# J. Runnoe
# Example of how to make a plot of a histogram.

# JCR import block
import numpy as np
from numpy.random import random
import matplotlib.pyplot as plt

if __name__=="__main__":
    # generate some random data
    data = 15*random(100000)

    # bin them
    # use the -0.5 offset to make bins centered on integer values
    # take it away to make bins centered on half values
    # also note how many bins this makes: histogram takes a
    # list of bin edges, for N bins you have to feed it n+1 edges
    hist, bin_edges = np.histogram(data,bins=np.arange(13)-0.5)

    # make a plot of the histogram
    # the reason the 0 bin is so low is because it covers
    # -0.5 to 0.5, and there are no negative numbers in my array
    fig, ax = plt.subplots(figsize=(8,6))
    ax.hist(data,bins=np.arange(13)-0.5,color="blue", align="mid",edgecolor='black', linewidth=1.2)
    ax.set_xlabel("Value",fontsize=22)
    ax.set_ylabel("Number",fontsize=22)
    fig.tight_layout()
    plt.show()
    plt.close(fig)
