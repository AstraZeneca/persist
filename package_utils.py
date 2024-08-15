# Package ready functions

import numpy as np
import pandas as pd
import dionysus as d
import pickle
import kneed
from scipy.sparse import csc_array
from sklearn.preprocessing import normalize



# == Network functions ==

def get_distances(num_wells, well_seperation, mesh_type):
    """ Computes an ordered list of the squared distances from the central vertex to nearest other vertices in a radial 
        hexagonal or square grid with a specified distance between adjacent vertices. 
        
        num_wells: integer, number of vertices
        well_seperation: float, distance between adjacent vertices
        mesh_type: either "hexagonal" or "square" 
    
        output: ordered list of squared distances from the central node
                to the nearest num_wells vertices in a radially expanding 
                network """
    
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
            distances = np.append( distances , 6*list(ring**2 + indexes**2 - ring*indexes) ) # np.append works fine with a list as the second argument :)
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


# == TDA Functions ==

# Convert a dionysus diagram object to a numpy array

def diagram_to_array(diagram, dimension=0):
    """ Takes in a dionysus diagram object and returns a m x 2 numpy array
        for the persistent features of the specified dimension, 
        with rows of the form
        [birth time] [death time] """
    
    res = np.zeros(( len(diagram[dimension]), 2 ))
    
    for i, pt in enumerate(diagram[dimension]):
        res[i,0] = pt.birth
        res[i,1] = pt.death
        
    # res = pd.DataFrame(res, columns=["dimension", "birth", "death"])
        
    return res 


# Compute the p-norm of a diagram 

def p_norm(diagram, p=2):
    """ takes in a diagram in array form and computes the p-norm.
    
        Array should be of the form
            birth death
              .     .
              .     .   
        with each row corresponding to a feature of the diagram. """
    
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
    
        values: n x 1 array, values of the function on each vertex 
        edges: n x 2 array indicating which edges to add in the simplicial complex """
    
    # Created filtration using list comprehension
    num_vertices = values.size
    f = d.Filtration( [ ([i], values[i]) for i in range(num_vertices) ] +
                      [ ([edges[i,0], edges[i,1]], np.min((values[edges[i,0]], values[edges[i,1]]))) for i in range(edges.shape[0]) ] )
    
    # Sort the simplices in descending order
    f.sort(reverse=True)
    
    return f

# === Distance to Measure Functions ===

# For use within run_persistence():
# Compute the dtm from a point to the weighted empirical measure

def distance_to_measure_point(weights, distances, m):
    """ Takes in a sorted vector of weights and associated squared distances for a single point, and
        computes the distance to measure for that point. 
        
        Designed to be called upon by distance_to_measure_weighted()
        
        weights: n x 1 array, sorted and normalised 
        distances: n x 1 of squared distances to add """
    
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
        distances from the central vertex in a radial hexagonal network, and computes the associated
        distance to measure. 
        
        Designed to be called upon by run_persistence()
        
        weights: n x 1 expression vector, should be normalised
        dmat: n x n distance matrix for wells, dmat[i,j] = d(well_i, well_j)
        network_distances: sorted vector of discrete distances, distances[i] = d(central vertex, ith nearest other vertex)
        m: in (0,1), threshold """
     
    # Create vector to store dtms
    res = np.zeros(dmat.shape[0])
    
    for i in range(dmat.shape[0]):        
        # Compute dtm for well i
        res[i] = distance_to_measure_point(weights, dmat[i,:], m)
    
    # Invert dtm before returning
    res = res.max() - res
    
    return res


# === Post norm computation functions ===

# Take in a list of norms and compute a cutoff

def find_knee(v):
    """ Takes in a Series v of norms and returns a cutoff point """
    
    elbow = kneed.KneeLocator(x=np.linspace(1, v.size, v.size), 
                              y=v.sort_values(ascending=False), 
                              curve="convex", direction="decreasing", S=1)
    
    return elbow.knee



# === Compute 0DPH ===

def run_persistence(data, p=2, m=0.1, mesh_type="hexagonal", sensitivity=1, 
                    metrics_storage_location=None, diagrams_storage_location=None, log_storage_location=None, 
                    notes=None, return_metrics=True, return_diagrams=True):
    
    """ Takes in expression and co-ordinate data for a set of wells from a single sample and computes the 0 dimensional persistence homology. 
        Optionally stores metrics (norms, ratios, ranks, and SVG calls) and diagrams in user-specified locations.
        
        data: pandas DataFrame of the form
               x , y , gene1 , gene2, ... , geneN
               .   .     .       .            .    
               .   .     .       .            .    
              where x, y = co-ordinates of each well
                    genei = expression of gene i in each well
        p: number from [1,np.inf] specifying which norm to compute from the 0D persistence diagram
        m: number in (0,1] speciying threshold to use in distance to measure computation
        
        return_metrics, return_diagrams: boolean flags specifiying what objects to return. If both are True,
                                         objects are returned in the order (metrics, diagrams) 
        metrics_storage_location: string specifying location to store metrics. Stored as a csv file
        diagrams_storage_location: string specifying where to store dictionary of persistence diagrams. Stored as a pkl file
        log_storage_location: Where to record details of persistence. Stored as a text file 
        notes: extra text to be appended on to the log 
        
        """
    
    # == Treat missing entries as zero expression ==
    
    data = data.fillna(0)
    
    # == Setup expression and co-ordinate data ==
    
    expression_array = normalize(np.array(data.iloc[:,2:]), axis=0, norm="l1") 
    expression_sparse = csc_array(expression_array)
    
    co_ordinates = np.array(data.iloc[:,:2])
    
    
    # == Network structure setup ==
    
    # Compute well-well distance matrix, so dmat[i,j] = d(well_i, well_j)
    # Also compute an upper triangular version for use in creating the filtration
    num_wells = data.shape[0]
    dmat_ut = np.zeros((num_wells,num_wells))
    dmat = np.zeros((num_wells,num_wells))
    for i in range(num_wells):
        dmat[i,:] = np.sqrt( ((co_ordinates[i,0] - co_ordinates[:,0])**2 + (co_ordinates[i,1] - co_ordinates[:,1])**2) )
        dmat_ut[i,i:] = np.sqrt( ((co_ordinates[i,0] - co_ordinates[i:,0])**2 + (co_ordinates[i,1] - co_ordinates[i:,1])**2) )
        
    # Compute the distance between adjacent wells as the smallest non-zero entry in dmat
    dunique = np.unique(dmat)
    well_sep = dunique[dunique>0].min()
    
    # Set distance threshold for treating wells as adjacent
    if mesh_type=="hexagonal":
        dthresh = (0.5*(1+np.sqrt(3))) * well_sep
    elif mesh_type=="square":
        dthresh = (0.5*(1+np.sqrt(2))) * well_sep
        
    # Compute an ordered list of the num_wells smallest squared distances from the central vertex that appear in a radial grid 
    grid_distances = get_distances(num_wells, well_sep, mesh_type)
    
    # Create a version of dmat with the actual distances replaced by the network distances 
    dmat_comp = np.zeros((num_wells,num_wells))
    for i in range(dmat.shape[0]):
        dmat_comp[i,dmat[i,:].argsort()] = grid_distances
    
    # Determine which pairs of wells are adjacent (i.e. which edges to include in the simplicial complex),
    # edge [i,j] is in the SC if it appears as a row in edges
    edges_to_add = np.where((dmat_ut>0) & (dmat_ut <= dthresh)) # Using upper triangular distance matrix ensures no duplication of edges (i.e. including both [i,j] and [j,i])
    edges = np.zeros((edges_to_add[0].size, 2))
    edges[:,0] = np.array(edges_to_add[0])
    edges[:,1] = np.array(edges_to_add[1])
    edges = edges.astype(int)  
    
    
    
    # == Compute persistent homology ==
    
    # Create objects for storing results 
    gene_list = data.columns[2:]
    num_genes = expression_array.shape[1]
    svg_calls = np.array(["No"]*(num_genes), dtype="<U3")# dtype="<U3" to allow "Yes" to be written later
    artefact_calls = np.array(["No"]*(num_genes), dtype="<U3")
    norms = np.zeros(num_genes)
    ratios = np.zeros(num_genes)
    
    record_diagrams = not diagrams_storage_location==None
    compute_diagrams = record_diagrams | return_diagrams
    if compute_diagrams:
        diagrams = {}

    # Determine which genes to compute persistence on
    compute_ph = np.where(expression_sparse.sum(axis=0)>0)[0]
    
    # Compute 0dPH for each gene with expression
    for i in compute_ph:
        expressed_wells = expression_array[:,i]>0
        
        # Determine which wells have expression (we will be skipping over the zero wells in the smoothing)
        weights = expression_sparse[expressed_wells,[i]].flatten()
        # Determine which indices sort the wells in descending order (- sign is as .argsort() by default sorts in ascending order)
        indices = (-weights).argsort()
        # Sort the weights by this
        weights = weights[indices]
        # Filter down dmat_comp columns to those wells with expression
        # Order columns by expression, to ensure that when we sort wells by distance to the well we are computing dtm for later on,
        # ties are broken by expression (by default .argsort() breaks ties by the order elements appear in the array pre-sorting)
        dists = dmat_comp[:,expressed_wells][:,indices]
        
        # Compute distance to measure
        dtm = distance_to_measure_weighted(weights=weights, dmat=dists, m=m)
        
        # Create filtration
        f = function_filtration(values=dtm, edges=edges)

        # Compute the 0 dimensional persistent homology
        ph = d.homology_persistence(f)
        dgm = d.init_diagrams(ph, f)
        dgm_array = diagram_to_array(dgm, 0)
        
        # Record the specified norm of persistence diagram, the persistence diagram itself, 
        # and the L_infinity / L_1 norm ratio for artifact detection
        norms[i] = p_norm(dgm_array, p=2)
        try: # try except clause is to handle the case where there are no features and so ||diagram||_1 = 0
            ratios[i] = p_norm(dgm_array, p=np.inf) / p_norm(dgm_array, p=1)
        except ZeroDivisionError:
            ratios[i] = 0
        if compute_diagrams:
            diagrams[gene_list[i]] = dgm
    
    # Call SVGs, using maximum curvature (via the kneed.KneeLocator() function) to determine the cutoff for SVG calling
    cutoff = find_knee(pd.Series(norms))
    ranks = np.zeros(norms.size)
    ranks[norms.argsort()] = np.linspace(norms.size, 1, norms.size) # think about it...
    if not cutoff==None:
        cutoff *= sensitivity
        svg_calls[ranks <= cutoff] = "Yes"
        
    artefact_calls[ratios > 0.9] = "Yes"
    
    # Create data frame of norm, ratio, rank and SVG call
    metrics = pd.DataFrame({"gene":pd.Series(gene_list), "CoSS":norms, "ratio":ratios, "gene_rank":ranks, "possible_artefact":artefact_calls, "svg":svg_calls}) 
    
    
    # == Store results if requested ==
    
    # Save persistence diagrams and norms
    if record_diagrams:
        with open(diagrams_storage_location, "wb") as f:
            pickle.dump(diagrams, f)
    
    if not metrics_storage_location==None:
        print("saving metrics")
        metrics.to_csv(metrics_storage_location, index=False)
        
    # Create log of details
    if not log_storage_location==None:
        with open(log_storage_location, "w") as f:
            f.write("p: " + str(p) + "\n")
            f.write("m: " + str(m) + "\n")  
            if not notes==None:
                f.write(notes + "\n")
    
    
    # == Return desired objects ==
    
    if return_metrics:
        if return_diagrams:
            return metrics, diagrams
        else:
            return metrics
    elif return_diagrams:
        return diagrams

    


