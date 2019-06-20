
# PREAMBLE:

import MDAnalysis
from MDAnalysis.analysis.align import rotation_matrix
import numpy as np
from numpy.linalg import *

# ----------------------------------------
def pearson_correlation_analysis(cartesian_covariance_matrix,average_node_positions,output_directory,distance_cutoff = 9999.,threshold=1e-4,max_iterations=100,guess=None,alpha=0.1):
	"""
	"""

        # ----------------------------------------
        # IO NAMING VARIABLES
        # ----------------------------------------
        variance_file_name = output_directory + 'node_variance.dat'  
        covariance_file_name = output_directory + 'node_covariance.dat'  
        correlation_file_name = output_directory + 'node_correlation.dat'  

        # ----------------------------------------
        # ASSIGNING MATRIX VARIABLES
        # ----------------------------------------
        nCartCoords = cartesian_covariance_matrix.shape[0]
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
                iIndex3 = iIndex+3
                for j in nNodes_range[i:]:
                        jIndex = j*3
                        jIndex3 = jIndex+3
                        node_covariance[i,j] = node_covariance[j,i] = np.trace(cartesian_covariance_matrix[iIndex:iIndex3,jIndex:jIndex3])
                        node_correlation[i,j] = node_correlation[j,i] = node_covariance[i,j]/np.sqrt(node_variance[i]*node_variance[j])

        print 'Pearson correlation coefficients have been calculated from the cartesian covariance.'
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
def linear_mutual_information_analysis(cartesian_covariance_matrix,average_node_positions,output_directory,distance_cutoff = 9999.,threshold=1e-4,max_iterations=100,guess=None,alpha=0.1):
	"""
	"""

        # ----------------------------------------
        # IO NAMING VARIABLES
        # ----------------------------------------
        linear_mutual_information_file_name = output_directory + 'linear_mutual_information.dat'

        # ----------------------------------------
        # ASSIGNING MATRIX VARIABLES
        # ----------------------------------------
        nCartCoords = cartesian_covariance_matrix.shape[0]
        nNodes = nCartCoords/3
        nNodes_range = range(nNodes)

        # ----------------------------------------
        # CALC LINEAR MUTUAL INFORMATION BETWEEN NODES BASED ON THE COVARIANCE OF THEIR CARTESIAN COORDINATES
        # ----------------------------------------
        MI = np.zeros((nNodes,nNodes),dtype=np.float64)
        for i in nNodes_range[:-1]:
                iIndex = i*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                iIndex3 = iIndex+3
                for j in nNodes_range[i+1:]:
                        jIndex = j*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                        jIndex3 = jIndex+3
                        temp = np.linalg.det(cartesian_covariance_matrix[iIndex:iIndex3,iIndex:iIndex3])*np.linalg.det(cartesian_covariance_matrix[jIndex:jIndex3,jIndex:jIndex3])  # compute numerator in argument of log of linear MI equation
                        idx = np.append(np.arange(iIndex,iIndex3,1),np.arange(jIndex,jIndex3,1))  # make list of indeces for the 2d X 2d C_ij matrix
                        MI[i,j] = MI[j,i] = 0.5*np.log(temp/np.linalg.det(cartesian_covariance_matrix[np.ix_(idx,idx)]))    # compute linear MI (eq 10 of Grubmuller 2005)

        MI += np.diagflat(np.full(nNodes,np.inf))

        print 'LMI has been calculated from the cartesian covariance.'
        # ----------------------------------------
        # SAVE LINEAR MUTUAL INFORMATION TO FILE
        # ----------------------------------------
        np.savetxt(linear_mutual_information_file_name,MI)

        # ----------------------------------------
        # RETURNING LINEAR MUTUAL INFORMATION MATRIX
        # ----------------------------------------
        return MI

def hENM_analysis(cartesian_covariance_matrix,average_node_positions,output_directory,distance_cutoff = 9999.,threshold=1e-4,max_iterations=100,guess=None,alpha=0.1):
        """
        """

        # ----------------------------------------
        # IO NAMING VARIABLES
        # ----------------------------------------
        hessian_file_name = output_directory + 'hessian.dat'

        # ----------------------------------------
        # ASSIGNING MATRIX VARIABLES
        # ----------------------------------------
        nCartCoords = cartesian_covariance_matrix.shape[0]
        nNodes = nCartCoords/3
        nNodes_range = range(nNodes)

        # ----------------------------------------
        # PREP THE INITIAL GUESS OF THE HESSIAN MATRIX; units: kcal mol^-1 \AA^-2
        # ----------------------------------------
        if guess != None and guess.shape == (nNodes,nNodes):
                hessian = guess
        else:
                hessian = -10.0*np.ones((nNodes,nNodes),dtype=np.float64)
                for i in nNodes_range:
                        hessian[i,i] = -np.sum(hessian[i,:])

        # ----------------------------------------
        # CALCULATING THE NODE PAIR DISTANCE VECTORS
        # ----------------------------------------
        dispVecs = np.empty((nNodes,nNodes,3),dtype=np.float64)
        for i in nNodes_range[:-1]:
                for j in nNodes_range[i+1:]:
                        diff = average_node_positions[i,:] - average_node_positions[j,:]
                        dispVecs[i,j,:] = dispVecs[j,i,:] = diff/np.linalg.norm(diff)

        # ----------------------------------------
        # CALCULATING THE INITIAL RESIDUAL FROM THE ORIGINAL COVARIANCE MATRIX
        # ----------------------------------------
        targetResidual = np.zeros((nNodes,nNodes),dtype=np.float64)
        for i in nNodes_range[:-1]:
                iIndex = i*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                iIndex3 = iIndex+3
                for j in nNodes_range[i+1:]:
                        jIndex = j*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                        jIndex3 = jIndex+3
                        temp = cartesian_covariance_matrix[iIndex:iIndex3,iIndex:iIndex3] + cartesian_covariance_matrix[jIndex:jIndex3,jIndex:jIndex3] - cartesian_covariance_matrix[iIndex:iIndex3,jIndex:jIndex3] - cartesian_covariance_matrix[iIndex:iIndex3,jIndex:jIndex3].T
                        targetResidual[i,j] = targetResidual[j,i] = 1.0/(np.dot(dispVecs[i,j].T,np.dat(temp,dispVecs[i,j])))
        
        # ----------------------------------------
        # ITERATIVELY IMPROVE THE HESSIAN TO MATCH THE RESIDUAL
        # ----------------------------------------
        step = 0 
        dev = thresh + 9999.
        converged = 'False'
        while step < max_iterations and converged == 'False':
                # project hessian in 3Nx3N
                hessian3N = np.zeros((nCartCoords,nCartCoords),dtype=np.float64)
                for i in nNodes_range[:-1]:
                        iIndex = i*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                        iIndex3 = iIndex+3
                        for j in nNodes_range[i+1:]:
                                jIndex = j*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                                jIndex3 = jIndex+3
                                hessian3N[iIndex:iIndex3,jIndex:jIndex3] = hessian3N[jIndex:jIndex3,iIndex:iIndex3] = hessian[i,j] * np.outer(dispVecs[i,j,:],dispVecs[i,j,:])
                for i in nNodes_range:
                        iIndex = i*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                        iIndex3 = iIndex+3
                        for j in nNodes_range:
                                jIndex = j*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                                if j != i:
                                        hessian3N[iIndex:iIndex3,iIndex:iIndex3] -= hessian3N[iIndex:iIndex3,jIndex:jIndex+3]

                # invert hessian using psuedo inverse; multiplying by thermal energy (kcal mol^-1) to remove energy units
                covar3N = np.linalg.pinv(hessian3N,rcond=1e-12)*0.593
                
                # take difference of tensor covars and project into separation vector
                diffMatrix = np.zeros((nNodes,nNodes),dtype=np.float64)
                for i in nNodes_range[:-1]:
                        iIndex = i*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                        iIndex3 = iIndex+3
                        for j in nNodes_range[i+1:]:
                                jIndex = j*3    # i index assuming that each node has 3 values sequentially populating the covar matrix 
                                jIndex3 = jIndex+3
                                temp = covar3N[iIndex:iIndex3,iIndex:iIndex3] + covar3N[jIndex:jIndex3,jIndex:jIndex3] - covar3N[iIndex:iIndex3,jIndex:jIndex3] - covar3N[iIndex:iIndex3,jIndex:jIndex3].T
                                diffMatrix[i,j] = diffMatrix[j,i] = (1.0/np.dot(dispVecs[i,j,:],np.dot(temp,dispVecs[i,j,:]))) - targetResidual[i,j]
                
                # check to see if converged
                dev = np.linalg.norm(diffMatrix)
                if dev < threshold:
                        converged = 'True'
                else: # update Hessian
                        for i in nNodes_range[:-1]:
                                for j in nNodes_range[i+1:]:
                                        hessian[i,j] += alpha * diffMatrix[i,j]
                                        if hessian[i,j] > 0.0:
                                                hessian[i,j] = hessian[j,i] = 0.0
                                        else:
                                                hessian[j,i] = hessian[i,j]
                        for i in nNodes_range:
                                hessian[i,i] = -np.sum(hessian[i,:])
                
                step += 1
                print step, dev

        print 'Hessian has been iteratively converged to match the original cartesian covariance.'
                        
        # ----------------------------------------
        # SAVE LINEAR MUTUAL INFORMATION TO FILE
        # ----------------------------------------
        np.savetxt(hessian_file_name,hessian)

        # ----------------------------------------
        # RETURNING LINEAR MUTUAL INFORMATION MATRIX
        # ----------------------------------------
        return hessian

