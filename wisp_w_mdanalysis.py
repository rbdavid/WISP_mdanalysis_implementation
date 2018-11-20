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

necessary_parameters = ['output_directory','node_selection_file','trajectory_functions_file','network_functions_file','visualization_functions_file','trajectory_analysis_boolean','pdb','visualization_frame_pdb','substrate_node_definition','substrate_selection_string','source_selection_string_list','sink_selection_string_list','number_of_paths']

all_parameters = ['output_directory','node_selection_file','trajectory_functions_file','network_functions_file','visualization_functions_file','trajectory_analysis_boolean','pdb','visualization_frame_pdb','substrate_node_definition','substrate_selection_string','source_selection_string_list','sink_selection_string_list','number_of_paths','nonstandard_substrates_selection','homemade_selections','alignment_selection','traj_list','summary_boolean','user_input_correlation_data','user_input_avg_node_node_distance_data','user_input_binary_contact_map_data','weight_by_avg_distance_boolean','node_sphere_radius','node_sphere_rgb','shortest_path_rgb','longest_path_rgb','shortest_path_radius','longest_path_radius','node_sphere_color_index','VMD_color_index_range','VMD_resolution','VMD_spline_smoothness']

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
        parameters['nonstandard_substrates_selection'] = None
        parameters['homemade_selections'] = None
        parameters['alignment_selection'] = 'all'
        parameters['traj_list'] = None
        parameters['contact_map_distance_cutoff'] = 99999.9
        parameters['user_input_correlation_data'] = None
        parameters['user_input_avg_node_node_distance_data'] = None
        parameters['user_input_binary_contact_map_data'] = None
        parameters['weight_by_avg_distance_boolean'] = False
        parameters['node_sphere_radius'] = 0.0
        parameters['node_sphere_rgb'] = (1.0,1.0,1.0)
        parameters['shortest_path_rgb'] = (1.0,1.0,1.0)
        parameters['longest_path_rgb'] = (1.0,1.0,1.0)
        parameters['node_sphere_color_index'] = 34
        parameters['VMD_color_index_range'] = (35,1056)
        parameters['summary_boolean'] = False 

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
			sys.exit()

def summary(summary_file_name):
        """ Function to create a text file that holds important information about the analysis that was just performed. Outputs the version of MDAnalysis, how to rerun the analysis, and the parameters used in the analysis.

        Usage:
            summary(summary_file_name)

        Arguments:
            summary_file_name: string object of the file name to be written that holds the summary information.

        """
	with open(summary_file_name,'w') as f:
		f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
		f.write('\nAtom selections analyzed have been written out to node_selections.txt\n')
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
        selection_list, source_indices, sink_indices = make_selections(u,parameters['output_directory'] + 'node_selections.txt',parameters['substrate_node_definition'],parameters['substrate_selection_string'],parameters['source_selection_string_list'],parameters['sink_selection_string_list'], nonstandard_substrates_selection = parameters['nonstandard_substrates_selection'], homemade_selections = parameters['homemade_selections'])
        print 'source node indices:', source_indices, '     sink node indices:', sink_indices

        # ----------------------------------------
        # TRAJECTORY ANALYSIS
        # ----------------------------------------
        if parameters['trajectory_analysis_boolean']:
                
                # ----------------------------------------
                # CALC THE AVERAGE NODE POSITIONS, VARIANCE OF NODE POSITIONS AWAY FROM THE AVERAGE POSITIONS, AND NODE-NODE CORRELATIONS 
                # ----------------------------------------
                correlation_array, variance_array, Node_trajectory, avg_Node_positions = correlation_trajectory_analysis(u,parameters['alignment_selection'],selection_list,parameters['traj_list'], average_file_name = parameters['output_directory']+ 'avg_node_coordinates.dat', variance_file_name = parameters['output_directory']+ 'node_node_variance.dat', correlation_file_name = parameters['output_directory']+ 'node_node_correlation.dat')
                print 'Finished calculating the correlation matrix.'

                # ----------------------------------------
                # CONTACT MAP CALCULATION
                # ----------------------------------------
                binary_contact_map, avg_node_node_distances = calc_contact_map(Node_trajectory,parameters['contact_map_distance_cutoff'],binary_contact_map_file_name = parameters['output_directory']+'binary_contact_map.dat', avg_node_node_distance_file_name = parameters['output_directory']+'avg_node_node_distance.dat')
                print 'Finished calculating the node-node contact map.'

        # ----------------------------------------
        # USER SPECIFIED DATA INPUT
        # ----------------------------------------
        else:
                print 'Loading in the user specified data files.'
                correlation_array = np.loadtxt(parameters['user_input_correlation_data'])

                if parameters['weight_by_avg_distance_boolean']:
                        avg_node_node_distances = np.loadtxt(parameters['user_input_avg_node_node_distance_data'])
                else:
                        binary_contact_map = np.loadtxt(parameters['user_input_binary_contact_map_data'])

        # ----------------------------------------
        # FUNCTIONALIZE THE CORRELATION MATRIX AND WEIGHTING THE FUNCTIONALIZED CORRELATION MATRIX BY THE CONTACT MAP
        # ----------------------------------------
        if parameters['weight_by_avg_distance_boolean']:
                func_corr_array = -np.log(np.fabs(correlation_array)) * avg_node_node_distances
        else:
                func_corr_array = -np.log(np.fabs(correlation_array)) * binary_contact_map

        np.savetxt(parameters['output_directory'] + 'func_correlation_matrix.dat',func_corr_array)
        
        print 'Finished functionalizing and weighting the correlation matrix.'
        
        # ----------------------------------------
        # COMPUTE THE ENSEMBLE OF PATHS THAT CONNECT THE SOURCE AND SINK RESIDUES 
        # ----------------------------------------
        paths = get_paths(func_corr_array, source_indices, sink_indices, parameters['number_of_paths'])
        print 'Finished calculating the pathways that connect the source and sink selections.'
        ### NEED TO OUTPUT THE PATHS

        # ----------------------------------------
        # CREATE THE VISUALIZATION STATES TO BE USED IN VMD 
        # ----------------------------------------
        u.load_new(parameters['visualization_frame_pdb'])

        create_vis_state(parameters['visualization_frame_pdb'],selection_list,paths,parameters['output_directory'] + 'wisp_pathways_vis_state.vmd', node_sphere_radius = parameters['node_sphere_radius'], node_sphere_rgb = parameters['node_sphere_rgb'], shortest_path_radius = parameters['shortest_path_radius'], shortest_path_rgb = parameters['shortest_path_rgb'], longest_path_radius = parameters['longest_path_radius'], longest_path_rgb = parameters['longest_path_rgb'], node_sphere_color_index = parameters['node_sphere_color_index'], VMD_color_index_range = parameters['VMD_color_index_range'], VMD_resolution = parameters['VMD_resolution'], VMD_spline_smoothness = parameters['VMD_spline_smoothness'])
        print 'Finished creating the vis state file.'
        
        # ----------------------------------------
        # OUTPUTTING SUMMARY INFORMATION
        # ----------------------------------------
        if parameters['summary_boolean']:
                summary(parameters['output_directory'] + 'summary.txt')

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

get_paths = importlib.import_module(parameters['network_functions_file'].split('.')[0],package=None).get_paths

create_vis_state = importlib.import_module(parameters['visualization_functions_file'].split('.')[0],package=None).create_vis_state

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

