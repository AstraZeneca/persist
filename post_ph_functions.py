""" Functions for use once a CoSS score is computed for each gene. """

import kneed
import numpy as np
import pandas as pd

# Take in a list of norms and compute a cutoff

def find_knee(v):
    """ Takes in a Series v of norms and returns a cutoff point to use for deciding 
        which genes are spatially variable.
    
        - v [1 x n pandas.Series]: Series of CoSS scores. """
    
    elbow = kneed.KneeLocator(x=np.linspace(1, v.size, v.size), 
                              y=v.sort_values(ascending=False), 
                              curve="convex", direction="decreasing", S=1)
    
    return elbow.knee