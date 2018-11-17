#!/home/rbdavid/bin/python

# ----------------------------------------
# PREAMBLE:
# ----------------------------------------

import sys
import os
import importlib
import numpy as np
import MDAnalysis

# ----------------------------------------
# VARIABLE DECLARATION: 
# ----------------------------------------

config_file = sys.argv[1]

necessary_parameters = ['output_directory','node_selection_file','trajectory_functions_file','trajectory_analysis_boolean','pdb','substrate_node_definition','substrate_selection_string','source_selection_string_list','sink_selection_string_list']  # ,'network_functions_file'
all_parameters = ['output_directory','node_selection_file','trajectory_functions_file','trajectory_analysis_boolean','pdb','substrate_node_definition','substrate_selection_string','selection_output_filename','nonstandard_substrates_selection','homemade_selections','alignment_selection','traj_list','summary_boolean','summary_filename','user_input_correlation_data','user_input_avg_node_node_distance_data','user_input_binary_contact_map_data']        # ,'network_functions_file'

# ----------------------------------------
# FUNCTIONS: 
# ----------------------------------------

def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
        """ Function to take config file and create/fill the parameter dictionary (created before function call). 
        
        Usage: 
            parameters = {}     # initialize the dictionary to be filled with keys and values
            config_parser(config_file)

        Arguments:
            config_file: string object that corresponds to the local or global position of the config file to be used for this analysis.

        """
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''
	
	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['selection_output_filename'] = 'selections.txt'
        parameters['nonstandard_substrates_selection'] = None
        parameters['homemade_selections'] = None
        parameters['alignment_selection'] = 'all'
        parameters['traj_list'] = None
        parameters['contact_map_distance_cutoff'] = 99999.9
	parameters['summary_boolean'] = False 
	parameters['summary_filename'] = 'summary.txt'

        parameters['user_input_correlation_data'] = None
        parameters['user_input_avg_node_node_distance_data'] = None
        parameters['user_input_binary_contact_map_data'] = None

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
			sys.exit()

def summary(summary_filename):
        """ Function to create a text file that holds important information about the analysis that was just performed. Outputs the version of MDAnalysis, how to rerun the analysis, and the parameters used in the analysis.

        Usage:
            summary(summary_filename)

        Arguments:
            summary_filename: string object of the file name to be written that holds the summary information.

        """
	with open(summary_filename,'w') as f:
		f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
		f.write('\nAtom selections analyzed have been written out to %s\n' %(parameters['selection_output_filename']))
		f.write('To recreate this analysis, run this line:\n')
		for i in range(len(sys.argv)):
			f.write('%s ' %(sys.argv[i]))
		f.write('\n\n')
		f.write('Parameters used:\n')
                for key, value in parameters.iteritems():
                        if key == '__builtins__':
                                continue
                        if type(value) == int or type(value) == float:
			        f.write("%s = %s\n" %(key,value))
                        else:
			        f.write("%s = '%s'\n" %(key,value))

def main():

        # ----------------------------------------
        # LOAD IN THE ANALYSIS UNIVERSE
        # ----------------------------------------
        u = MDAnalysis.Universe(parameters['pdb'])
        
        # ----------------------------------------
        # MAKE ATOM SELECTIONS FOR USER-DEFINED NODE DEFINITIONS
        # ----------------------------------------
        selection_list, source_indices, sink_indices = make_selections(u,parameters['output_directory'] + parameters['selection_output_filename'],parameters['substrate_node_definition'],parameters['substrate_selection_string'],parameters['source_selection_string_list'],parameters['sink_selection_string_list'], nonstandard_substrates_selection = parameters['nonstandard_substrates_selection'], homemade_selections = parameters['homemade_selections'])
        print 'source node indices:', source_indices, '     sink node indices:', sink_indices

        # ----------------------------------------
        # TRAJECTORY ANALYSIS
        # ----------------------------------------
        if parameters['trajectory_analysis_boolean']:
                
                # ----------------------------------------
                # CALC THE AVERAGE NODE POSITIONS, VARIANCE OF NODE POSITIONS AWAY FROM THE AVERAGE POSITIONS, AND NODE-NODE CORRELATIONS 
                # ----------------------------------------
                correlation_array, variance_array, Node_trajectory, avg_Node_positions = correlation_trajectory_analysis(u,parameters['alignment_selection'],selection_list,parameters['traj_list'], average_file_name = parameters['output_directory']+ 'avg_node_coordinates.dat', variance_file_name = parameters['output_directory']+ 'node_node_variance.dat', correlation_file_name = parameters['output_directory']+ 'node_node_correlation.dat')

                # ----------------------------------------
                # CONTACT MAP CALCULATION
                # ----------------------------------------
                binary_contact_map, avg_node_node_distances = calc_contact_map(Node_trajectory,parameters['contact_map_distance_cutoff'],binary_contact_map_file_name = parameters['output_directory']+'binary_contact_map.dat', avg_node_node_distance_file_name = parameters['output_directory']+'avg_node_node_distance.dat')

        # ----------------------------------------
        # USER SPECIFIED DATA INPUT
        # ----------------------------------------
        else:
                correlation_array = np.loadtxt(parameters['user_input_correlation_data'])
                if parameters['weight_by_avg_distance_boolean']:
                        avg_node_node_distances = np.loadtxt(parameters['user_input_avg_node_node_distance_data'])
                else:
                        binary_contact_map = np.loadtxt(parameters['user_input_binary_contact_map_data'])

        # ----------------------------------------
        # WEIGHTING THE CORRELATION MATRIX BY THE CONTACT MAP
        # ----------------------------------------
        if parameters['weight_by_avg_distance_boolean']:
                max_distance = np.max(avg_node_node_distances)
                correlation_array *= avg_node_node_distances/max_distance
        else:
                correlation_array *= binary_contact_map

        # ----------------------------------------
        # FUNCTIONALIZE THE CORRELATION MATRIX
        # ----------------------------------------
        func_corr_array = -np.log(np.fabs(correlation_array))
        np.savetxt(parameters['output_directory'] + parameters['functionalized_weighted_correlation_matrix_file_name'],func_corr_array)

        # ----------------------------------------
        # COMPUTE THE ENSEMBLE OF PATHS THAT CONNECT THE SOURCE AND SINK RESIDUES 
        # ----------------------------------------
        


        # ----------------------------------------
        # CREATE THE VISUALIZATION STATES TO BE USED IN VMD 
        # ----------------------------------------



        # ----------------------------------------
        # OUTPUTTING SUMMARY INFORMATION
        # ----------------------------------------
        if parameters['summary_boolean']:
                summary(parameters['output_directory'] + parameters['summary_filename'])

# ----------------------------------------
# CREATING PARAMETER DICTIONARY
# ----------------------------------------
parameters = {}
config_parser(config_file)

# ----------------------------------------
# LOADING IN NECESSARY FUNCTIONS FROM MODULE FILES
# ----------------------------------------

make_selections = importlib.import_module(parameters['node_selection_file'].split('.')[0],package=None).make_selections

if parameters['trajectory_analysis_boolean']:
        correlation_trajectory_analysis = importlib.import_module(parameters['trajectory_functions_file'].split('.')[0],package=None).correlation_trajectory_analysis
        calc_contact_map = importlib.import_module(parameters['trajectory_functions_file'].split('.')[0],package=None).calc_contact_map

#make_selections = importlib.import_module(parameters['network_functions_file'].split('.')[0],package=None).make_selections

# ----------------------------------------
# MAIN
# ----------------------------------------
if parameters['output_directory'][-1] != os.sep:
        parameters['output_directory'] += os.sep

if os.path.exists(parameters['output_directory']):
        print 'The output directory, ', parameters['output_directory'], 'already exists. Please select a different directory name for output.'
        sys.exit()
else:
        os.mkdir(parameters['output_directory'])

# ----------------------------------------
# MAIN
# ----------------------------------------
if __name__ == '__main__':
	main()










#def get_average_coordinates(avg_universe,node_type,alignment_selection_string,user_defined_selection,nNodes):
#	"""
#	"""
#	if node_type == 'atomic':
#		node_average_coordinates = avg_universe.select_atoms(user_defined_selection).positions
#		align_average_coordinates = avg_universe.select_atoms(alignment_selection_string).positions
#		return node_average_coordinates, align_average_coordinates
#	elif node_type == 'residue_COM':
#		align_average_coordinates = avg_universe.select_atoms(alignment_selection_string).positions
#		selection = avg_universe.select_atoms(user_defined_selection)	# does not necessitate only protein selection; assumes only one node per residue
#		
#		node_average_coordinates = np.zeros((nNodes,3),dtype=np.float32)
#		for i in range(nNodes):
#			node_average_coordinates[i] = selection.residues[i].center_of_mass()
#		return node_average_coordinates, align_average_coordinates
#
#	elif node_type == 'backbone_COM':
#		align_average_coordinates = avg_universe.select_atoms(alignment_selection_string).positions
#		selection = avg_universe.select_atoms('(backbone or nucleicbackbone) and '+user_defined_selection)	# does not necessitate only protein selection; assumes only one node per residue
#		
#		node_average_coordinates = np.zeros((nNodes,3),dtype=np.float32)
#		for i in range(nNodes):
#			node_average_coordinates[i] = selection[selection.resids == selection.residues[i].resid].center_of_mass()
#		return node_average_coordinates, align_average_coordinates
#
#	elif node_type == 'sidechain_COM':
#		align_average_coordinates = avg_universe.select_atoms(alignment_selection_string).positions
#		selection = avg_universe.select_atoms('not (backbone or nucleicbackbone) and '++user_defined_selection)	# does not necessitate only protein selection; assumes only one node per residue
#		
#		node_average_coordinates = np.zeros((nNodes,3),dtype=np.float32)
#		for i in range(nNodes):
#			node_average_coordinates[i] = selection[selection.resids == selection.residues[i].resid].center_of_mass()
#		return node_average_coordinates, align_average_coordinates
#
#	### STILL TO DO:
#	# MORE COMPLEX, USER DEFINED ATOM SELECTIONS; good example: multiple sites per residue (i.e. backbone and sidechain COM for each residue of interest, with the correct node ordering based on resid)

#def iterative_avg(u,trajectory_list,u_all,u_align,u_analysis_selection_list,node_type,nNodes,nAlign,nSteps,trajectory_stepping_value=1,threshold=1E-5,maxIter=100):
#	"""
#	"""
#	# ----------------------------------------
#	# INITIALIZING ARRAYS
#	all_coord = zeros((nSteps,nNodes,3),dtype=np.float32)
#	avg_coord = zeros((nNodes,3),dtype=np.float32)
#	all_align = zeros((nSteps,nAlign,3),dtype=np.float32)
#	avg_align = zeros((nAlign,3),dtype=np.float32)
#
#	# ----------------------------------------
#	# ANALYZING TRAJECTORY
#	print 'Beginning trajectory analysis'
#	temp = 0
#	if node_type == 'atomic':
#		for a in trajectory_list:
#			u.load_new(a)
#			for ts in u.trajectory[::trajectory_stepping_value]:
#				u_all.translate(-u_align.center_of_mass())
#				avgAlign += u_align.positions
#				all_align[tem] = u_align.positions
#				for i in range(len(u_analysis_selection_list)):
#					avgCoord += u_analysis_selection_list[i].pos
#					all_coord[temp] = u_analysis_selection_list[i].pos
#				temp += 1
#	elif 'COM' in node_type:
#		for a in trajectory_list:
#			u.load_new(a)
#			for ts in u.trajectory[::trajectory_stepping_value]:
#				u_all.translate(-u_align.center_of_mass())
#				avgAlign += u_align.positions
#				all_align[tem] = u_align.positions
#				for i in range(len(u_analysis_selection_list)):
#					avgCoord += u_analysis_selection_list[i].center_of_mass()
#					all_coord[temp] = u_analysis_selection_list[i].center_of_mass()
#				temp += 1
#	
#	avg_align /= nSteps
#	avg_coord /= nSteps
#	
#	iteration = 0 
#	residual = threshold + 100.
#	print 'Beginning iterative process of calculating average coordinates and aligning to the new average.'
#	while residual > threshold and iteration < maxIter:
#		temp_avg_coord = np.zeros((nNodes,3),dtype=np.float32)
#		temp_avg_align = np.zeros((nAlign,3),dtype=np.float32)
#		for i in range(nSteps):
#			R, d = rotation_matrix(all_align[i,:,:],avg_align)
#			all_align[i,:,:] = dot(all_align[i,:,:],R.T)
#			all_coord[i,:,:] = dot(all_coord[i,:,:],R.T)
#			temp_avg_align += all_align[i,:,:]
#			temp_avg_coord += all_coord[i,:,:]
#		temp_avg_align /= nSteps
#		temp_avg_coord /= nSteps
#		residual = RMSD(avg_align, temp_avg_align, nAlign)
#		rmsd_nodes = RMSD(avg_coord, temp_avg_coord, nNodes)
#		avg_align = temp_avg_align
#		avg_coord = temp_avg_coord
#		print 'Finished with Iteration step', iteration, ', RMSD of alignment landmark atoms is', residual, ' RMSD of Node selections is', rmsd_nodes
#		iteration += 1
#	
#	return avg_coord, avg_align











#        # ----------------------------------------
#        # LOAD IN THE AVERAGE UNIVERSE
#        avg = MDAnalysis.Universe(parameters['average_pdb'])
#        avg_align = u.select_atoms(parameters['alignment_selection'])
#        avg_all = u.select_atoms('all')
#        avg_all.translate(-avg_align.center_of_mass())
#        pos0 = avg_align.positions
#
#        # ----------------------------------------
#        # CALCULATING THE CARTESIAN COVARIANCE, AVERAGE, AND VARIANCE ARRAYS
#        if selections[1] == 'COM':
#                print 'Performing a covariance analysis of the cartesian coordinates of the center of mass of residues defined in the selection_file.'
#                covariance_array, average_array, variance_array = calc_cart_covar_matrix(u,parameters['traj_list'],parameters['alignment_selection'],pos0,selection_list,nNodes)
#               
#                np.savetxt('cart_covar.' + parameters['data_output_filename'],covariance_array)
#                np.savetxt('cart_var.' + parameters['data_output_filename'],variance_array)
#                np.savetxt('cart_avg.' + parameters['data_output_filename'],average_array)
#
#        else:
#                print 'Requested to calc the covariance array of something other than the center of mass of residues defined in the selection_file... Currently, this code is unable to do anything other than COM'
#                sys.exit()
#
#        # ----------------------------------------
#        # CALCULATING THE DISTANCE CORRELATION MATRIX OF NODE-NODE PAIRS
#        distance_correlation_matrix = np.zeros((nNodes,nNodes),dtype=np.float64)
#        distance_variance_array = np.zeros((nNodes),dtype=np.float64)
#        for i in nNodes_range:
#                dim1 = i*3
#                distance_variance_array[i] = sum(variance_array[dim1:dim1+3])   # summing the variances of the cartesian dimensions of the node i; <r_{i}(t)**2>
#                for j in nNodes_range[i:]:
#                        dim2 = j*3
#                        distance_correlation_matrix[i,j] = covariance_array[dim1,dim2] + covariance_array[dim1+1,dim2+1] + covariance_array[dim1+2,dim2+2]      # taking the trace of the cartesian covariance array for the desired dimensions of nodes i and j; <r_{i}(t) dot r_{j}(t)> - <r_{i}(t)> dot <r_{j}(t)> 
#
#        for i in nNodes_range:
#                for j in nNodes_range[i:]:
#                        distance_correlation_matrix[i,j] /= np.sqrt(distance_variance_array[i]*distance_variance_array[j])     # finishing the distance correlation matrix by dividing by the product of the variances of nodes i,j;  <r_{i}(t) dot r_{j}(t)> - <r_{i}(t)> dot <r_{j}(t)> / sqrt((<x_{i}(t)**2> - <x_{i}(t)>**2)*(<x_{j}(t)**2> - <x_{j}(t)>**2))
#                        distance_correlation_matrix[j,i] = distance_correlation_matrix[i,j] # filling in the bottom traingle of this matrix
#
#        np.savetxt('dist_corr.' + parameters['data_output_filename'],distance_correlation_matrix)
#        np.savetxt('dist_var.' + parameters['data_output_filename'],distance_variance_array)
#
#        # ----------------------------------------
#        # FUNCTIONALIZE THE DISTANCE CORRELATION MATRIX (USING WISP EQUATION)
#        if parameters['functionalize_distance_correlation_boolean']:
#                func_distance_correlation_matrix = np.zeros((nNodes,nNodes),dtype=np.float64)
#                print 'Beginning to functionalize the distance correlation matrix.'
#                for i in nNodes_range:
#                        for j in nNodes_range[i:]:
#                                func_distance_correlation_matrix[i,j] = -np.log(np.fabs(distance_correlation_matrix[i,j]))
#                                func_distance_correlation_matrix[j,i] = func_distance_correlation_matrix[i,j]
#                np.savetxt('func_dist_corr.' + parameters['data_output_filename'],func_distance_correlation_matrix)
#        else:
#                print 'No functionalization of the distance correlation matrix is being performed.'
#
#
#def main():
#
#####code breakdown
## 1) calc iterative, average structure using an alignment landmark if required; only necessary for the atom selection to be analyzed
## 2) calc the distance between node(t) and <node> for every t; produces d matrix
## 3) calc the correlation between d(i) and d(j)
## 4) read in contact map...
## 5) 
#
#	# ----------------------------------------
#	# INITIALIZING MDANALYSIS UNIVERSES 
#	u = MDAnalysis.Universe(parameters['pdb'])
#	u_all = u.select_atoms('all')
#	u_align = u.select_atoms(parameters['alignment_selection'])
#	u_analysis_selection, u_analysis_selection_list, nNodes = make_selection(u,parameters['node_type'],parameters['user_defined_selection'],parameters['selection_output_filename'])
#	nAlign = u_align.n_atoms
#
#	# ----------------------------------------
#	# INITIALIZING PARAMETER VARIABLES
#	if parameters['trajectory_file_list'] == None:
#		start = int(parameters['start'])
#		end = int(parameters['end'])
#
#		trajectory_file_list = []
#		temp = start
#		while temp <= end:
#			trajectory_file_list.append(parameters['trajectory_location_string']%(temp,temp))	# HUGE SOURCE OF ANNOYANCE... NOT EVERYONE ORGANIZES TRAJECTORIES IN THIS SPECIFIC FORMAT/LOCATION/ETC... potential solution: read in a file that contains the global location for all trajectory files to be analyzed
#	else: print 'User has read in a trajectory file list.'
#
#	# ----------------------------------------
#	# DETERMINING nSTEPS TO BE ANALYZED 
#	if parameters['nSteps'] == None:
#		nSteps = 0
#		for i in trajectory_file_list:
#			u.load_new(i)
#			nSteps += len(range(0,u.trajectory.n_frames,parameters['trajectory_stepping_value']))
#	else:
#		nSteps = int(parameters['nSteps'])
#
#	# ----------------------------------------
#	# CALCULATING ITERATIVE AVERAGE STRUCTURE
#	if parameters['user_defined_average_structure'] == None:
#		node_average_coordinates, align_average_coordinates = iterative_avg(u,trajectory_file_list,u_all,u_align,u_analysis_selection_list,parameters['node_type'],nNodes,nAlign,nSteps,trajectory_stepping_value=parameters['trajectory_stepping_value'])
#	else:
#		avg = MDAnalysis.Universe(parameters['user_defined_average_structure'])		# assumes a pdb file of the alignment and analysis selections, at the very least.
#		node_average_coordinates, align_average_coordinates = get_average_coordinates(avg,parameters['node_type'],parameters['alignment_selection'],parameters['user_defined_selection'],nNodes)
#	
#	# ----------------------------------------
#	# OUTPUTTING AVG COORDINATES OF ALIGNMENT LANDMARK AND ANALYSIS SELECTION
#	###
#
#	# ----------------------------------------
#	# READ IN USER DEFINED AVG DISTANCE MATRIX BETWEEN EACH NODE PAIR; TO BE USED TO SIMPLIFY THE CORRELATION MATRIX
#
#	# ----------------------------------------
#	# CALCULATE DISTANCE BETWEEN AVERAGE NODE POSITION AND NODE POSITION AT TIME t
#



