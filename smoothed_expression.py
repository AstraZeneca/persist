""" Functions related to computing a smoothed expression value. """

import numpy as np

# Compute the dtm from a point to the weighted empirical measure

def distance_to_measure_point(weights, distances, m):
    """ Takes in a sorted vector of weights and associated squared distances for a single point, and
        computes the distance to measure for that point. 
        
        Designed to be called upon by distance_to_measure_weighted()
        
        - weights [n x 1 numpy.array]: Weights, sorted by distance from central vertex and normalised. 
        - distances [n x 1 numpy.array]: Squared distances from central vertex to surrounding vertices. """
    
    # Create indexes to sort via distance (distances will come in order of weight so this will break ties)
    indices = distances.argsort()
    # Sort the weights via distance (ascending), ties will end up broken via weight (descending)
    weights = weights[indices]
    
    # Set up counters
    mass = 0 # current mass
    nn = 0 # number of nodes added
    
    # Compute dtm
    while mass < m:
        mass += weights[nn]
        nn += 1
        
    res = distances[indices[:nn]].sum()
    
    return res


# Compute the weighted dtm for all points on a mesh

def distance_to_measure_weighted(weights, dmat, m):
    """ Takes in a vector of expression weights, a well-well distance matrix, and an ordered list of 
        distances from the central vertex in a radial network, and computes the associated
        distance to measure for each vertex. 
        
        Designed to be called upon by run_persistence()
        
        - weights [n x 1 numpy.array]: Expression in each well, should be normalised.
        - dmat [n x n numpy.array]: Well-well distance matrix, dmat[i,j] = d(well_i, well_j).
        - network_distances [n x 1 numpy.array]: Sorted array of discrete distances, distances[i] = d(central vertex, ith nearest other vertex).
        - m [float in (0,1)]: Threshold. """
     
    # Create vector to store dtms
    res = np.zeros(dmat.shape[0])
    
    for i in range(dmat.shape[0]):        
        # Compute dtm for well i
        res[i] = distance_to_measure_point(weights, dmat[i,:], m)
    
    # Invert dtm before returning
    res = res.max() - res
    
    return res