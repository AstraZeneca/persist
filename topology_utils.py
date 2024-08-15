""" Functions related to Topological Data Analysis - Simiplicial Complicies, Persistent Homology, and Persistence Diagrams. """

import numpy as np
import dionysus as d

# Convert a dionysus diagram object to a numpy array

def diagram_to_array(diagram, dimension=0):
    """ Takes in a dionysus diagram object and returns a m x 2 numpy array
        for the persistent features of the specified dimension, 
        with rows of the form
        [birth time] [death time] 
        
        - diagram [list of dionysus._dionysus.Diagram objects]: Dionysus diagram list.
        - dimension [int]: Which dimension to take the persistent features from. """
    
    res = np.zeros(( len(diagram[dimension]), 2 ))
    
    for i, pt in enumerate(diagram[dimension]):
        res[i,0] = pt.birth
        res[i,1] = pt.death
        
    return res 


# Compute the p-norm of a diagram 

def p_norm(diagram, p=2):
    """ Computes the p-norm of a diagram in array form.
    
        - diagram [m X 2 numpy.array]: Array with the birth and death times of features from a persistence diagram. Columns should be in order 'birth', 'death'.
        - p [int]: Specifies which norm to compute. """
    
    # We first remove features with no death time 
    # In particular, for a 0-dim persistence module, there will always be one feature
    # (the single cc that is the point cloud) with infinite lifetime
    features = diagram[diagram[:,1]!=np.inf]
    
    # Compute the liftime of each feature
    lifetimes = np.abs(features[:,1] - features[:,0])
    
    # Compute the p-norm of the vector of lifetimes
    if len(lifetimes)==0:
        res = 0
    elif p==np.inf:
        res = lifetimes.max()
    else:
        res = (np.sum(lifetimes**p))**(1/p)
    
    return res


# Construct the filtration for a specific set of vertices, edges, and values defined on vertices

def function_filtration(values, edges):
    """ Takes in an adjacency structure and function values for a set of vertices, and computes the upper star
        filtration for this function.
        
        index of each vertex = value of the function at that vertex
        index of an edge = value of the function on the lowest vertex of the edge
        Higher dimensional faces not included as they do not affect the 0D PH
    
        - values [n x 1 numpy.array]: Values of the function on each vertex 
        - edges [n x 2 numpy.array]: Which edges to add in the simplicial complex, each row [a,b] adds edge [a,b] to the vertex. """
    
    # Created filtration using list comprehension
    num_vertices = values.size
    f = d.Filtration( [ ([i], values[i]) for i in range(num_vertices) ] +
                      [ ([edges[i,0], edges[i,1]], np.min((values[edges[i,0]], values[edges[i,1]]))) for i in range(edges.shape[0]) ] )
    
    # Sort the simplices in descending order
    f.sort(reverse=True)
    
    return f