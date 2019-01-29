#!/home/rbdavid/bin/python

### NOTE: THE INTERPRETATION OF CENTRALITY METRICS IS WHOLLY 
### DEPENDENT ON THE ADJACENCY MATRIX USED TO DEFINE THE
### CONNECTIVITY OF YOUR NETWORK. 

### THIS CODE (AND ANALYSES WITHIN) DOES NOT SUPPLANT THE
### USER'S CAREFUL CONSIDERATION OF THEIR DATA. 

# ----------------------------------------
# PREAMBLE:
# ----------------------------------------

import sys
import os
import importlib
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------
# VARIABLE DECLARATION: 
# ----------------------------------------

simply_formatted_pathways = sys.argv[1]
adjacency_matrix = np.loadtxt(sys.argv[2])
contact_map = sys.argv[3]

# ----------------------------------------
# FUNCTIONS: 
# ----------------------------------------

def degree_centrality(adjacency_matrix):
        """
        """
        return np.sum(adjacency_matrix,axis=0)

def betweenness_centrality(pathways_file,nNodes):
        """
        """
        bc_array = np.zeros(nNodes,dtype=np.float64)
        nPathways = 0
        with open(pathways_file) as File:
                for line in File:
                        if line[0] == '#':
                                continue
                        nPathways += 1
                        data = line.split()
                        path_nodes = data[1:]
                        for i in path_nodes:
                                bc_array[i] += 1

        bc_array /= nPathways

        return bc_array

def eigenvector_centrality(adjacency_matrix, contact_map = None, damping_factor = np.inf):
        """
        """
        # apply adjacency matrix damping, to look at local correlations.
        if contact_map != None and damping_factor != np.inf:
                adjacency_matrix = np.multiply(adjacency_matrix,np.exp(-contact_map[i,j]/damping_factor))

        # diagonalize the Adjacency matrix
        eigenvalues, eigenvectors = np.linalg.eigh(adjacency_matrix)
        idx = eigenvalues.argsort()[::-1]   # sort by largest to smallest
        eigenvalues = eigenvalues[idx]      # apply sorting
        eigenvectors = eigenvectors[:,idx]    # apply sorting
        # compute centrality vector, c
        c = np.dot(adjacency_matrix,eigenvectors[:,0])/eigenvalues[0]   # only calculating centrality from the first eigenvector

        return c, eigenvalues, eigenvectors

# ----------------------------------------
# MAIN: 
# ----------------------------------------



