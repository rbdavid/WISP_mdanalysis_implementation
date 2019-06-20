
# PREAMBLE:

import MDAnalysis
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.analysis.distances import distance_array
import numpy as np
from numpy.linalg import *

# ----------------------------------------
def euclid_dist(x,y):
        """ Calculates the Euclidian Distance between two arrays of the same size
        Usage: dist,dist2 = euclid_dist(x,y)
            
        Arguments:
        x, y: numpy arrays with the same size
        """

        dist2 = np.sum(np.square(x-y))
        dist = dist2**0.5   # the MSD should never be negative, so using **0.5 rather than np.sqrt is safe
        return dist, dist2

# ----------------------------------------
def RMSD(x,y,n):
	""" Calculates the Root Mean Squared Distance between two arrays of the same size

	Usage: rmsd = RMSD(x,y,n)

	Arguments:
	x, y: numpy arrays with the same shape (n X 3)
	n: number of particles being summed over; ex: number of atoms in the atom selection being analyzed;
		if n = 1, this function calculates the distance between x and y arrays

	"""
	
	return (np.sum(np.square(x-y))/n)**0.5	# the MSD should never be negative, so using **0.5 rather than np.sqrt is safe

# ----------------------------------------
def traj_alignment_and_averaging(universe, alignment_selection, selection_list, node_definition, trajectory_list, output_directory, convergence_threshold = 1E-5, maximum_num_iterations = 100, step = 1):
	"""
	"""
	
        # ----------------------------------------
        # IO NAMING VARIABLES
        # ----------------------------------------
        average_file_name = output_directory + 'average_node_positions.dat'
        
        # ----------------------------------------
	# CREATING ALIGNMENT SELECTIONS
        # ----------------------------------------
        u_all = universe.select_atoms('all')
        u_alignment = universe.select_atoms(alignment_selection)
        nAlign = u_alignment.n_atoms
        nNodes = len(selection_list)
        nNodes_range = range(nNodes)

        # ----------------------------------------
	# ANALYZE TRAJECTORIES TO COLLECT THE NECESSARY POSITION DATA
        # ----------------------------------------
        all_pos_Align = []
        all_pos_Nodes = []
        nSteps = 0
        for traj in trajectory_list:
                print 'Loading trajectory', traj
                universe.load_new(traj)
                if len(universe.trajectory)%step != 0:
                        print 'User defined step size is not a factor of the number of frames in ', traj, '. The user is not getting the desired step size.'
                        sys.exit()
                nSteps += len(universe.trajectory)/step
                for ts in universe.trajectory[::step]:
                        u_all.translate(-u_alignment.center_of_mass())  # removing center of mass (of the alignment landmark) translation
                        all_pos_Align.append(u_alignment.positions)
                        temp = []
                        for i in nNodes_range:
                                if node_definition.upper() == 'COM':
                                        temp.append(selection_list[i].center_of_mass())
                                elif node_definition.upper() == 'ATOMIC':
                                        temp.append(selection_list[i].position)
                        all_pos_Nodes.append(temp)

        print 'Analyzed', nSteps, 'frames.'
        
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
                temp_avg_pos_Align = np.zeros((nAlign,3),dtype=np.float64)
                temp_avg_pos_Nodes = np.zeros((nNodes,3),dtype=np.float64)

                for ts in nSteps_range:
                        
                        # calculate the rotation matrix (and distance) between frame i's alignment postions to the average alignment positions
                        R, d = rotation_matrix(all_pos_Align[ts,:,:],avg_pos_Align)      
                        
                        # take the dot product between frame i's alignment positions and the calculated rotation matrix; overwrite frame i's positions with the rotated postions
                        all_pos_Align[ts,:,:] = np.dot(all_pos_Align[ts,:,:],R.T)       
                        all_pos_Nodes[ts,:,:] = np.dot(all_pos_Nodes[ts,:,:],R.T)

                        # running sum of alignment positions to calculate a new average
                        temp_avg_pos_Align += all_pos_Align[ts,:,:]      
                        temp_avg_pos_Nodes += all_pos_Nodes[ts,:,:]
              
                # finish calculating the new averages
                temp_avg_pos_Align /= nSteps
                temp_avg_pos_Nodes /= nSteps

                # calculate the difference between old and new averages
                residual = RMSD(avg_pos_Align,temp_avg_pos_Align,nAlign)
                analysis_RMSD = RMSD(avg_pos_Nodes,temp_avg_pos_Nodes,nNodes)
                
                # increment the iteration number and assign new averages to old averages variables
                iteration += 1
                avg_pos_Align = np.copy(temp_avg_pos_Align)
                avg_pos_Nodes = np.copy(temp_avg_pos_Nodes)

                print 'Iteration ', iteration, ': RMSD btw alignment landmarks: ', residual,', RMSD btw Node Positions: ', analysis_RMSD
        
        print 'Finished calculating the average structure using the iterative averaging approach. Outputting the average node positions to file.'

        # ----------------------------------------
        # SAVE AVERAGE NODE POSITIONS TO FILE
        # ----------------------------------------
        np.savetxt(average_file_name,avg_pos_Nodes,header='Shape: nNodes x 3 (%d x 3)'%(nNodes))
      
        # ----------------------------------------
        # RETURNING THE ALIGNED NODE TRAJECTORY AND AVERAGE NODE POSITIONS
        # ----------------------------------------
        return all_pos_Nodes, avg_pos_Nodes

# ----------------------------------------
def cartesian_covariance(node_trajectory, avg_node_positions, output_directory):
	"""
	"""
        # ----------------------------------------
        # IO NAMING VARIABLES
        # ----------------------------------------
        variance_file_name = output_directory + 'cartesian_variance.dat'
        covariance_file_name = output_directory + 'cartesian_covariance.dat'
        
        # ----------------------------------------
        # CALCULATING THE COVARIANCE OF NODE CARTESIAN COORDINATES
        # ----------------------------------------
        nSteps = node_trajectory.shape[0]
        nSteps_range = range(nSteps)
        nNodes = node_trajectory.shape[1]
        nCartCoords = 3*nNodes
        nCartCoords_range = range(nCartCoords)

        # ----------------------------------------
        # CALC VARIANCE AND COVARIANCE OF CARTESIAN COORDINATES
        # ----------------------------------------
        print 'Beginning cartesian covariance analysis.'
        xyz_node_covariance = np.zeros((nCartCoords,nCartCoords),dtype=np.float64)
        for ts in nSteps_range:
                flatten_positions = all_pos_Nodes[ts].flatten()
                xyz_node_covariance += np.dot(flatten_positions.reshape(nCartCoords,1),flatten_positions.reshape(1,nCartCoords))

        xyz_node_covariance /= nSteps
        flatten_avg_positions = avg_pos_Nodes.flatten()
        xyz_node_covariance -= np.dot(flatten_avg_positions.reshape(nCartCoords,1),flatten_avg_positions.reshape(1,nCartCoords))
        xyz_node_variance = np.diag(xyz_node_covariance)
        
        print 'Finished calculating the variance and covariance of node cartesian coordinates. Outputting the two arrays to file. The covariance matrix file can be reused in subsequent analyses of various adjacency matrices, assuming studying the same range of frames, selections, etc.'

        # ----------------------------------------
        # SAVE THE COVARIANCE OF NODE CARTESIAN COORDINATES TO FILE
        # ----------------------------------------
        np.savetxt(variance_file_name,xyz_node_variance)
        np.savetxt(covariance_file_name,xyz_node_covariance)
       
        # ----------------------------------------
        # COVARIANCE OF NODE CARTESIAN COORDINATES
        # ----------------------------------------
        return xyz_node_covariance

# ----------------------------------------
def calc_contact_map(trajectory_data, output_directory, distance_cutoff=9999.):
        """
        """
	
        # ----------------------------------------
        # IO NAMING VARIABLES
        # ----------------------------------------
        binary_contact_map_file_name = output_directory + 'binary_contact_map.dat'
        avg_node_node_distance_file_name = output_directory + 'avg_distance_contact_map.dat'
        
        # ----------------------------------------
	# MEASURING THE AVERAGE DISTANCE BETWEEN NODES - add error analysis?
        # ----------------------------------------
        nSteps = trajectory_data.shape[0]
        nSteps_range = range(nSteps)
        nNodes = trajectory_data.shape[1]
        nNodes_range = range(nNodes)
        avg_node_node_distance_array = np.zeros((nNodes,nNodes),dtype=np.float64)
        temp_avg_node_node_distance_array = np.zeros((nNodes,nNodes),dtype=np.float64)
        print 'Starting to calculate the average distance between nodes.'
        for ts in nSteps_range:
                avg_node_node_distance_array += distance_array(trajectory[ts,:,:],trajectory[ts,:,:],result=temp_avg_node_node_distance_array)
                #for i in nNodes_range[:-1]:
                #        node1_pos = trajectory_data[ts,i]
                #        for j in nNodes_range[i+1:]:
                #                dist,dist2 = euclid_dist(node1_pos,trajectory_data[ts,j])
                #                avg_node_node_distance_array[i,j] += dist
       
        avg_node_node_distance_array /= nSteps
        
        # ----------------------------------------
	# CREATING THE BINARY DISTANCE CONTACT MAP
        # ----------------------------------------
        binary_node_node_distance_array = np.zeros((nNodes,nNodes),np.int)
        for i in nNodes_range[:-1]:
                for j in nNodes_range[i+1:]:
                        #avg_node_node_distance_array[j,i] = avg_node_node_distance_array[i,j]
                        if avg_node_node_distance_array[i,j] < distance_cutoff:
                                binary_node_node_distance_array[i,j] = binary_node_node_distance_array[j,i] = 1

        # ----------------------------------------
        # SAVING THE AVERAGE AND BINARY CONTACT MAPS
        # ----------------------------------------
        np.savetxt(binary_contact_map_file_name,binary_node_node_distance_array,fmt='%i')
        np.savetxt(avg_node_node_distance_file_name,avg_node_node_distance_array)

        # ----------------------------------------
        # RETURNING BINARY AND AVERAGE DISTANCE CONTACT MAPS
        # ----------------------------------------
        return binary_node_node_distance_array, avg_node_node_distance_array

