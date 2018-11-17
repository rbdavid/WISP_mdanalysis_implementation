
# PREAMBLE:

import MDAnalysis
from MDAnalysis.analysis.align import rotation_matrix
import numpy as np
from numpy.linalg import *

def euclid_dist(x,y):
        """ Calculates the Euclidian Distance between two arrays of the same size
        Usage: dist,dist2 = euclid_dist(x,y)
            
        Arguments:
        x, y: numpy arrays with the same size
        """

        dist2 = np.sum(np.square(x-y))
        dist = dist2**0.5   # the MSD should never be negative, so using **0.5 rather than np.sqrt is safe
        return dist, dist2

def RMSD(x,y,n):
	""" Calculates the Root Mean Squared Distance between two arrays of the same size

	Usage: rmsd = RMSD(x,y,n)

	Arguments:
	x, y: numpy arrays with the same shape (n X 3)
	n: number of particles being summed over; ex: number of atoms in the atom selection being analyzed;
		if n = 1, this function calculates the distance between x and y arrays

	"""
	
	return (np.sum(np.square(x-y))/n)**0.5	# the MSD should never be negative, so using **0.5 rather than np.sqrt is safe

def correlation_trajectory_analysis(universe, alignment_selection, selection_list, trajectory_list, convergence_threshold = 1E-5, maximum_num_iterations = 100, average_file_name=None, variance_file_name=None, correlation_file_name=None):
	"""
	"""
	
        # ----------------------------------------
	# CREATING ALIGNMENT SELECTIONS
        # ----------------------------------------
        u_all = universe.select_atoms('all')
        u_alignment = universe.select_atoms(alignment_selection)
        nAlign = u_alignment.n_atoms
        nNodes = len(selection_list)
        nNodes_range = range(nNodes)

        print nAlign, nNodes

        # ----------------------------------------
	# ANALYZE TRAJECTORIES TO COLLECT THE NECESSARY POSITION DATA
        # ----------------------------------------
        all_pos_Align = []
        all_pos_Nodes = []
        nSteps = 0
        for traj in trajectory_list:
                print 'Loading trajectory', traj
                universe.load_new(traj)
                nSteps += len(universe.trajectory)
                for ts in universe.trajectory:
                        u_all.translate(-u_alignment.center_of_mass())
                        all_pos_Align.append(u_alignment.positions)
                        temp = []
                        for i in nNodes_range:
                                temp.append(selection_list[i].center_of_mass())
                        all_pos_Nodes.append(temp)

        print 'Analyzed', nSteps, 'frames. Does this match up with expectation?'
        
        all_pos_Align = np.array(all_pos_Align)
        all_pos_Nodes = np.array(all_pos_Nodes)

        avg_pos_Align = np.sum(all_pos_Align,axis=0)/nSteps
        avg_pos_Nodes = np.sum(all_pos_Nodes,axis=0)/nSteps

        # ----------------------------------------
	# ITERATIVE ALIGNMENT TO AVERAGE ALIGNMENT POSITIONS
        # ----------------------------------------
        iteration = 0 
        residual = convergence_threshold + 9999.
        nSteps_range = range(nSteps)
        print 'Beginning the iterative process of aligning to the average alignment positions, calculating new positions, and recalculating the average positions'
        while residual > convergence_threshold and iteration < maximum_num_iterations:
                temp_avg_pos_Align = np.zeros((nAlign,3),dtype=np.float32)
                temp_avg_pos_Nodes = np.zeros((nNodes,3),dtype=np.float32)

                for ts in nSteps_range:
                        R, d = rotation_matrix(all_pos_Align[ts,:,:],avg_pos_Align)      # calculate the rotation matrix (and distance) between frame i's alignment postions to the average alignment positions
                        all_pos_Align[ts,:,:] = np.dot(all_pos_Align[ts,:,:],R.T)       # take the dot product between frame i's alignment positions and the calculated rotation matrix; overwrite frame i's positions with the rotated postions
                        all_pos_Nodes[ts,:,:] = np.dot(all_pos_Nodes[ts,:,:],R.T) # take the dot product between frame i's analysis positions and the calculated rotation matrix; overwrite frame i's positions with the rotated postions
                        temp_avg_pos_Align += all_pos_Align[ts,:,:]      # running sum of alignment positions to calculate a new average
                        temp_avg_pos_Nodes += all_pos_Nodes[ts,:,:]      # running sum of analysis positions to calculate a new average
               
                temp_avg_pos_Align /= nSteps
                temp_avg_pos_Nodes /= nSteps
                residual = RMSD(avg_pos_Align,temp_avg_pos_Align,nAlign)
                analysis_RMSD = RMSD(avg_pos_Nodes,temp_avg_pos_Nodes,nNodes)
                iteration += 1
                avg_pos_Align = temp_avg_pos_Align
                avg_pos_Nodes = temp_avg_pos_Nodes
                print 'Iteration ', iteration, ': RMSD btw alignment landmarks: ', residual,', RMSD btw Node Positions: ', analysis_RMSD
        
        print 'Finished calculating the average structure using the iterative averaging approach. Outputting the average node positions to file.'

        # ----------------------------------------
        # SAVE AVERAGE NODE POSITIONS TO FILE
        # ----------------------------------------
        if average_file_name != None:
                np.savetxt(average_file_name,avg_pos_Nodes)

        # ----------------------------------------
        # CALC CORRELATIONS OF DISTANCES OF NODES FROM THEIR AVERAGE POSITIONS
        # ----------------------------------------
        print 'Beginning correlation analysis.'
        xyz_node_variance = np.zeros((nNodes,3),dtype=np.float64)
        xyz_node_covariance = np.zeros((nNodes,nNodes,3),dtype=np.float64)

        for ts in nSteps_range:
                for i in nNodes_range:
                        xyz_node_variance[i] += all_pos_Nodes[ts,i,:]**2
                        for j in nNodes_range[i:]:      # looking at diagonal and top triangle elements of the covariance matrix
                                xyz_node_covariance[i,j,:] += all_pos_Nodes[ts,i,:]*all_pos_Nodes[ts,j,:]

        xyz_node_variance /= nSteps
        xyz_node_variance -= avg_pos_Nodes**2
        xyz_node_covariance /= nSteps

        for i in nNodes_range:
                for j in nNodes_range[i:]:
                        xyz_node_covariance[i,j] -= avg_pos_Nodes[i]*avg_pos_Nodes[j]

        distance_node_variance = np.zeros((nNodes),dtype=np.float64)
        distance_node_correlation = np.zeros((nNodes,nNodes),dtype=np.float64)
        for i in nNodes_range:
                distance_node_variance[i] = np.sum(xyz_node_variance[i,:])
        
        for i in nNodes_range:
                for j in nNodes_range[i:]:
                        distance_node_correlation[i,j] = np.sum(xyz_node_covariance[i,j])/np.sqrt(distance_node_variance[i]*distance_node_variance[j])
                        distance_node_correlation[j,i] = distance_node_correlation[i,j]

        # ----------------------------------------
        # SAVE VARIANCE AND CORRELATION OF NODE-NODE POSITIONS TO FILE
        # ----------------------------------------
        if variance_file_name != None:
                np.savetxt(variance_file_name,distance_node_variance)
        if correlation_file_name != None:
                np.savetxt(correlation_file_name,distance_node_correlation)

        # ----------------------------------------
        # RETURNING THE VARIANCE AND CORRELATION ARRAYS TO BE USED AGAIN
        # ----------------------------------------
        return distance_node_correlation, distance_node_variance, all_pos_Nodes, avg_pos_Nodes
        

def calc_contact_map(trajectory_data,distance_cutoff,binary_contact_map_file_name=None,avg_node_node_distance_file_name=None):
        """
        """
	
        # ----------------------------------------
	# MEASURING THE AVERAGE DISTANCE BETWEEN NODES
        # ----------------------------------------
        nSteps = len(trajectory_data)
        nSteps_range = range(nSteps)
        nNodes = len(trajectory_data[0])
        nNodes_range = range(nNodes)
        avg_node_node_distance_array = np.zeros((nNodes,nNodes),dtype=np.float64)
        print 'Starting to calculate the average distance between nodes.'
        for ts in nSteps_range:
                for i in nNodes_range[:-1]:
                        node1_pos = trajectory_data[ts,i]
                        for j in nNodes_range[i+1:]:
                                dist,dist2 = euclid_dist(node1_pos,trajectory_data[ts,j])
                                avg_node_node_distance_array[i,j] += dist
       
        avg_node_node_distance_array /= nSteps
        
        # ----------------------------------------
	# CREATING THE BINARY DISTANCE CONTACT MAP
        # ----------------------------------------
        binary_node_node_distance_array = np.zeros((nNodes,nNodes),np.int)
        for i in nNodes_range[:-1]:
                for j in nNodes_range[i+1:]:
                        avg_node_node_distance_array[j,i] = avg_node_node_distance_array[i,j]
                        if avg_node_node_distance_array[i,j] < distance_cutoff:
                                binary_node_node_distance_array[i,j] = 1
                                binary_node_node_distance_array[j,i] = 1

        if binary_contact_map_file_name != None:
                np.savetxt(binary_contact_map_file_name,binary_node_node_distance_array)
        if avg_node_node_distance_file_name != None:
                np.savetxt(avg_node_node_distance_file_name,avg_node_node_distance_array)

        return binary_node_node_distance_array, avg_node_node_distance_array

#def average_structure_calc(pdb,alignment_selection,analysis_selection,traj_list,avg_structure_file_name, convergence_threshold = 1E-5, maximum_num_iterations = 100, wrapping_boolean = False):
#        """
#        """
#
#        # ----------------------------------------
#        # LOAD IN AND CREATE ATOM SELECTIONS IN THE ANALYSIS UNIVERSE OBJECT
#        u = MDAnalysis.Universe(pdb)
#        u_all = u.select_atoms('all')
#        u_align = u.select_atoms(alignment_selection)  # MDAnalysis atom selection string formatting required.
#        u_analysis = u.select_atoms(analysis_selection)  # MDAnalysis atom selection string formatting required.
#        nAtoms_align = u_align.n_atoms
#        nAtoms_analysis = u_analysis.n_atoms
#
#        # ----------------------------------------
#        # ANALYZE TRAJECTORIES TO COLLECT THE NECESSARY POSITION DATA
#        all_pos_align = []
#        all_pos_analysis = []
#        nSteps = 0
#        nResidues_range = range(u_analysis.n_residues)
#        for i in traj_list:
#                print 'Loading trajectory', i
#                u.load_new(i)
#                nSteps += len(u.trajectory)
#                for ts in u.trajectory:
#                        u_all.translate(-u_align.center_of_mass())
#                        
#                        if wrapping_boolean:
#                                dimensions = u.dimensions[:3]
#                                for j in nResidues_range:
#                                        COM = u_analysis.residues[j].center_of_mass()
#                                        t = wrapping(COM,dimensions)
#                                        u_analysis.residues[j].atoms.translate(t)
#                        
#                        all_pos_align.append(u_align.positions)
#                        all_pos_analysis.append(u_analysis.positions)
#        
#        print 'Analyzed', nSteps, 'frames. Does this match up with expectation?'
#
#        all_pos_align = np.array(all_pos_align)
#        all_pos_analysis = np.array(all_pos_analysis)
#
#        avg_pos_align = np.sum(all_pos_align,axis=0)/nSteps
#        avg_pos_analysis = np.sum(all_pos_analysis,axis=0)/nSteps
#
#        # ----------------------------------------
#        # ITERATIVE ALIGNMENT TO AVERAGE ALIGNMENT POSITIONS
#        iteration = 0 
#        residual = convergence_threshold + 9999.
#        nSteps_range = range(nSteps)
#        print 'Beginning the iterative process of aligning to the average alignment positions, calculating new positions, and recalculating the average positions'
#        while residual > convergence_threshold and iteration < maximum_num_iterations:
#                temp_avg_pos_align = np.zeros((nAtoms_align,3),dtype=np.float32)
#                temp_avg_pos_analysis = np.zeros((nAtoms_analysis,3),dtype=np.float32)
#
#                for i in nSteps_range:
#                        R, d = rotation_matrix(all_pos_align[i,:,:],avg_pos_align)      # calculate the rotation matrix (and distance) between frame i's alignment postions to the average alignment positions
#                        all_pos_align[i,:,:] = dot_prod(all_pos_align[i,:,:],R.T)       # take the dot product between frame i's alignment positions and the calculated rotation matrix; overwrite frame i's positions with the rotated postions
#                        all_pos_analysis[i,:,:] = dot_prod(all_pos_analysis[i,:,:],R.T) # take the dot product between frame i's analysis positions and the calculated rotation matrix; overwrite frame i's positions with the rotated postions
#                        temp_avg_pos_align += all_pos_align[i,:,:]          # running sum of alignment positions to calculate a new average
#                        temp_avg_pos_analysis += all_pos_analysis[i,:,:]    # running sum of analysis positions to calculate a new average
#               
#                temp_avg_pos_align /= nSteps
#                temp_avg_pos_analysis /= nSteps
#                residual = RMSD(avg_pos_align,temp_avg_pos_align,nAtoms_align)
#                analysis_RMSD = RMSD(avg_pos_analysis,temp_avg_pos_analysis,nAtoms_analysis)
#                iteration += 1
#                avg_pos_align = temp_avg_pos_align
#                avg_pos_analysis = temp_avg_pos_analysis
#                print 'Iteration ', iteration, ': RMSD btw alignment landmarks: ', residual,', RMSD btw analysis atoms: ', analysis_RMSD
#        
#        print 'Finished calculating the average structure using the iterative averaging approach. Outputting the average structure to a pdb now.'
#        # ----------------------------------------
#        # LOAD IN AND CREATE ATOM SELECTIONS IN THE RESULTS UNIVERSE OBJECT
#        avg = MDAnalysis.Universe(pdb)
#        avg_analysis = avg.select_atoms(analysis_selection)
#        avg_analysis.positions = avg_pos_analysis
#        avg_analysis.write(avg_structure_file_name)
#
#        return avg_pos_analysis, all_pos_analysis

