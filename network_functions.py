""" Functions related to the network structure of the co-ordinates in the input data. """

import numpy as np

def get_distances(num_wells, well_seperation, mesh_type):
    """ Computes an ordered list of the squared distances from the central vertex to the nearest other vertices in a radial 
        hexagonal or square grid with a specified distance between adjacent vertices. 
        
        - num_wells [int]: Number of vertices to compute distances to.
        - well_seperation [float]: Distance between adjacent vertices
        - mesh_type [str in ['hexagonal', 'square']]: What type of mesh the vertices lie on. """
    
    # compute how many "rings" needed to acheive the number of vertices
    R_max = 0
    included_wells = 1
    if mesh_type=="hexagonal":
        while included_wells < num_wells:
            R_max += 1
            included_wells += 6*R_max
    elif mesh_type=="square":
        while included_wells < num_wells:
            R_max += 1
            included_wells += 8*R_max
    
    # compute squared distances from central vertex
    distances = np.array([0]) # 0 = distance from central vertex to itself
    for ring in range(1, R_max+1):
        if mesh_type=="hexagonal":
            # distances computed using cosine rule
            indexes = np.linspace(0,ring-1,ring)
            distances = np.append( distances , 6*list(ring**2 + indexes**2 - ring*indexes) ) 
        elif mesh_type=="square":
            # central points
            distances = np.append( distances, 4*[ring**2] )
            # corner points
            distances = np.append( distances, 4*[2*(ring**2)] )
            # intermediate points, distances comptued using pythagoras' theorem
            if ring >= 2:
                indexes = np.linspace(1,ring-1,ring-1)
                distances = np.append( distances , 8*list(ring**2 + indexes**2) )
        
    # sort distances, multiply by squared vertex seperation, and return correct number of distances
    distances = np.sort(distances)[:num_wells]*(well_seperation**2)
    
    return distances 