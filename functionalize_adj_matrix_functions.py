
# PREAMBLE:

import MDAnalysis
from MDAnalysis.analysis.align import rotation_matrix
import numpy as np
from numpy.linalg import *

# ----------------------------------------
def func_pearson_correlation(adjacency_matrix,output_directory, contact_map = None, Lambda = None):
	"""

        NOTES about using the binary contact map... no distance weighting and, so, will get trivial/unphysical pathways...
	"""
        
        functionalized_correlation_file_name = output_directory + 'func_pearson_correlation.dat'
        functionalized_correlation = -np.log(np.fabs(adjacency_matrix))
        
        if contact_map is not None:
                functionalized_correlation *= contact_map

        np.savetxt(functionalized_correlation_file_name,functionalized_correlation)

        return functionalized_correlation

# ----------------------------------------
def generalized_correlation_coefficient_calc(adjacency_matrix,output_directory, contact_map = None, Lambda = None):
	"""
	"""

        rMI_file_name = output_directory + 'gen_corr_coefficients.dat'
        functionalized_rMI_file_name = output_directory + 'func_gen_corr_coefficients.dat'

        rMI = np.sqrt(1.0 - np.exp(-2.0*adjacency_matrix/3))
       
        if contact_map is not None and Lambda is not None:
                print 'Applying a node-node distance-dependent exponential decay to the generalized correlation coefficient matrix (r_{MI}). The exponential damping factor (%.3f angstroms) is being used. This damping propogates into the functionalized adjacency matrix that will be used for pathway-finding analysis.'
                for i in range(rMI.shape[0]):
                        for j in range(rMI.shape[0]):
                                rMI[i,j] *= np.exp(-contact_map[i,j]/Lambda)

        functionalized_rMI = -np.log(rMI)
       
        np.savetxt(rMI_file_name,rMI)
        np.savetxt(functionalized_rMI_file_name,functionalized_rMI)
        return functionalized_rMI

# ----------------------------------------
def func_k_matrix(adjacency_matrix,output_directory, contact_map = None, Lambda = None):
	"""
	"""
        functionalized_k_matrix_file_name = output_directory + 'func_k_matrix.dat'
        
        nNodes = adjacency_matrix.shape[0]
        nNodes_range = range(nNodes)
        sum_k_matrix = np.sum(adjacency_matrix,axis=0)
        
        functionalized_k_matrix = np.zeros((nNodes,nNodes),dtype=np.float64)
        
        for i in nNodes_range:
                for j in nNodes_range:
                        functionalized_k_matrix[i,j] = -np.log(adjacency_matrix[i,j]/np.sqrt(sum_k_matrix[i]*sum_k_matrix[j]))

        np.savetxt(functionalized_k_matrix_file_name,functionalized_k_matrix)
        return functionalized_k_matrix

# ----------------------------------------
def do_nothing(adjacency_matrix,output_directory, contact_map = None, Lambda = None):
	"""
	"""
        return adjacency_matrix

# ----------------------------------------
def func_hessian(adjacency_matrix,output_directory, contact_map = None, Lambda = None):
	"""
	"""
        
        functionalized_correlation_file_name = output_directory + 'func_hessian.dat'
        nNodes = adjacency_matrix.shape[0]
        nNodes_range = range(nNodes)
        for i in nNodes_range[:-1]:
                for j in nNodes_range[i+1:]:
                        adjacency_matrix[i,j] /= np.sqrt(adjacency_matrix[i,i]*adjacency_matrix[j,j])
                        adjacency_matrix[j,i] = adjacency_matrix[i,j]
                adjacency_matrix[i,i] = 1.
        adjacency_matrix[-1,-1] = 1.

        functionalized_correlation = -np.log(np.fabs(adjacency_matrix))
        
        if contact_map is not None:
                functionalized_correlation *= contact_map

        np.savetxt(functionalized_correlation_file_name,functionalized_correlation)

        return functionalized_correlation

