
# PREAMBLE:

import MDAnalysis
from MDAnalysis.analysis.align import rotation_matrix
import numpy as np
from numpy.linalg import *

# ----------------------------------------
def func_pearson_correlation(adjacency_matrix,output_directory, contact_map = None):
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
def generalized_correlation_coefficient_calc(adjacency_matrix,output_directory, contact_map = None):
	"""
	"""
        rMI_file_name = output_directory + 'gen_corr_coefficients.dat'
        functionalized_rMI_file_name = output_directory + 'func_gen_corr_coefficients.dat'
        rMI = np.sqrt(1.0 - np.exp(-2.0*adjacency_matrix/3))
        functionalized_rMI = -np.log(rMI)
        np.savetxt(rMI_file_name,rMI)
        np.savetxt(functionalized_rMI_file_name,functionalized_rMI)
        return functionalized_rMI

# ----------------------------------------
def do_nothing(adjacency_matrix,output_directory, contact_map = None):
	"""
	"""
        return adjacency_matrix

