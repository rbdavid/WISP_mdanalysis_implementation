
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
        
        if contact_map != None:
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
       
        if contact_map != None and Lambda != None:
                print 'Applying a node-node distance-dependent exponential decay to the generalized correlation coefficient matrix (r_{MI}). The exponential damping factor (%.3f angstroms) is being used. This damping propogates into the functionalized adjacency matrix that will be used for pathway-finding analysis.'
                for i in range(rMI.shape[0]):
                        for j in range(rMI.shape[0]):
                                rMI[i,j] *= np.exp(-contact_map[i,j]/Lambda)

        functionalized_rMI = -np.log(rMI)
       
        np.savetxt(rMI_file_name,rMI)
        np.savetxt(functionalized_rMI_file_name,functionalized_rMI)
        return functionalized_rMI

# ----------------------------------------
def do_nothing(adjacency_matrix,output_directory, contact_map = None):
	"""
	"""
        return adjacency_matrix

