
# PREAMBLE:

import MDAnalysis
from MDAnalysis.analysis.align import rotation_matrix
import numpy as np
from numpy.linalg import *

# ----------------------------------------
def pearson_correlation_analysis(cartesian_covariance_matrix,output_directory):
	"""
	"""

        # ----------------------------------------
        # IO NAMING VARIABLES
        # ----------------------------------------
        variance_file_name = output_directory + 'node_variance.dat'  
        covariance_file_name = output_directory + 'node_covariance.dat'  
        correlation_file_name = output_directory + 'node_correlation.dat'  

        # ----------------------------------------
        # ASSIGNING TRAJECTORY VARIABLES
        # ----------------------------------------
        nCartCoords = cartesian_covariance_matrix.shape[0]
        if nCartCoords%3 != 0:
                print 'The user has read in a covariance array that does not have the expected number of elements. Please check the file.'
                sys.exit()
        nNodes = nCartCoords/3
        nNodes_range = range(nNodes)

        # ----------------------------------------
        # CALC VARIANCE AND COVARIANCE OF NODE DISPLACEMENT AWAY FROM AVERAGE POSITIONS
        # ----------------------------------------
        print 'Beginning Pearson correlation coefficient analysis.'
        node_variance = np.zeros((nNodes),dtype=np.float64)
        node_covariance = np.zeros((nNodes,nNodes),dtype=np.float64)
        node_correlation = np.zeros((nNodes,nNodes),dtype=np.float64)

        cartesian_variance_matrix = np.diagonal(cartesian_covariance_matrix)
        for i in nNodes_range:
                iIndex = i*3
                node_variance[i] = np.sum(cartesian_variance_matrix[iIndex:iIndex+3])
        
        for i in nNodes_range:
                iIndex = i*3
                for j in nNodes_range[i:]:
                        jIndex = j*3
                        node_covariance[i,j] = node_covariance[j,i] = np.trace(cartesian_covariance_matrix[iIndex:iIndex+3,jIndex:jIndex+3])
                        node_correlation[i,j] = node_correlation[j,i] = node_covariance[i,j]/np.sqrt(node_variance[i]*node_variance[j])

        # ----------------------------------------
        # SAVE VARIANCE, COVARIANCE, CORRELATION OF NODE DISPLACEMENTS TO FILE
        # ----------------------------------------
        np.savetxt(variance_file_name,node_variance)
        np.savetxt(covariance_file_name,node_covariance)
        np.savetxt(correlation_file_name,node_correlation)

        # ----------------------------------------
        # RETURNING CORRELATION ARRAY
        # ----------------------------------------
        return node_correlation

# ----------------------------------------
def linear_mutual_information_analysis(cartesian_covariance_matrix,output_directory):
	"""
	"""

        # ----------------------------------------
        # IO NAMING VARIABLES
        # ----------------------------------------
        linear_mutual_information_file_name = output_directory + 'linear_mutual_information.dat'

        # ----------------------------------------
        # ASSIGNING TRAJECTORY VARIABLES
        # ----------------------------------------
        nCartCoords = cartesian_covariance_matrix.shape[0]
        if nCartCoords%3 != 0:
                print 'The user has read in a covariance array that does not have the expected number of elements. Please check the file.'
                sys.exit()
        nNodes = nCartCoords/3
        nNodes_range = range(nNodes)
        nCartCoords_range = range(nCartCoords)

        # ----------------------------------------
        # CALC LINEAR MUTUAL INFORMATION BETWEEN NODES BASED ON THE COVARIANCE OF THEIR CARTESIAN COORDINATES
        # ----------------------------------------
        MI = np.zeros((nNodes,nNodes),dtype=np.float64)
        for i in nNodes_range[:-1]:
                iIndex = i*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                for j in nNodes_range[i+1:]:
                        jIndex = j*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                        temp = np.linalg.det(cartesian_covariance_matrix[iIndex:iIndex+3,iIndex:iIndex+3])*np.linalg.det(cartesian_covariance_matrix[jIndex:jIndex+3,jIndex:jIndex+3])  # compute numerator in argument of log of linear MI equation
                        idx = np.append(np.arange(iIndex,iIndex+3,1),np.arange(jIndex,jIndex+3,1))  # make list of indeces for the 2d X 2d C_ij matrix
                        MI[i,j] = MI[j,i] = 0.5*np.log(temp/np.linalg.det(cartesian_covariance_matrix[np.ix_(idx,idx)]))    # compute linear MI (eq 10 of Grubmuller 2005)

        MI += np.diagflat(np.full(nNodes,np.inf))

        # ----------------------------------------
        # SAVE LINEAR MUTUAL INFORMATION TO FILE
        # ----------------------------------------
        np.savetxt(linear_mutual_information_file_name,MI)

        # ----------------------------------------
        # RETURNING LINEAR MUTUAL INFORMATION MATRIX
        # ----------------------------------------
        return MI

def hENM_analysis(cartesian_covariance_matrix,output_directory):
        """
        """
        
        ####
        return cartesian_covariance

