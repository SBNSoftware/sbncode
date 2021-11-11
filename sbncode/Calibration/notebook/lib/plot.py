import numpy as np

def makehist(plt, N, bins, **kwargs):
    plt.hist(bins[:-1], bins, weights=N, **kwargs)
