
# PREAMBLE:

import numpy as np
from numpy.linalg import *
import networkx

def get_paths(correlation_matrix,source,sink,number_of_paths):
        """
        """
        G = networkx.Graph(data = correlation_matrix)
        print 'Calculating paths...'
        shortest_length, shortest_path = get_shortest_path_length(correlation_matrix,sources,sinks,G)

        num_paths = 1

        path = [shortest_length]
        path.extend(shortest_path)
        paths = [path]

        cutoff = shortest_length

        cutoff_yields_max_num_paths_below_target = 0 
        cutoff_yields_min_num_paths_above_target = 1000000.0 

        while num_paths < number_of_paths:
                



