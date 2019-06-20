
# ----------------------------------------
# PREAMBLE:
# ----------------------------------------

import importlib
import MDAnalysis
import sys

# ----------------------------------------
# FUNCTIONS: 
# ----------------------------------------
def config_parser(config_file,parameters):	# Function to take config file and create/fill the parameter dictionary 
        """ Function to take config file and create/fill the parameter dictionary (created before function call). 
        
        Usage: 
            parameters = {}     # initialize the dictionary to be filled with keys and values
            config_parser(config_file)

        Arguments:
            config_file: string object that corresponds to the local or global position of the config file to be used for this analysis.

        """
        necessary_parameters = ['output_directory','node_selection_file','trajectory_functions_file','adjacency_matrix_functions_file','func_adjacency_matrix_functions_file','network_functions_file','visualization_functions_file','trajectory_analysis_boolean','adjacency_matrix_style','pdb','visualization_frame_pdb','substrate_node_definition','substrate_selection_string','source_selection_string_list','sink_selection_string_list','number_of_paths']

        all_parameters = ['output_directory','node_selection_file','trajectory_functions_file','adjacency_matrix_functions_file','func_adjacency_matrix_functions_file','network_functions_file','visualization_functions_file','trajectory_analysis_boolean','adjacency_matrix_style','pdb','visualization_frame_pdb','substrate_node_definition','substrate_selection_string','source_selection_string_list','sink_selection_string_list','number_of_paths','calc_contact_map_boolean','summary_boolean','nonstandard_substrates_selection','homemade_selections','alignment_selection','traj_list','trajectory_step','contact_map_distance_cutoff','which_contact_map','user_input_cartesian_covariance_matrix','user_input_adjacency_matrix','user_input_contact_map','node_sphere_radius','node_sphere_rgb','shortest_path_rgb','longest_path_rgb','shortest_path_radius','longest_path_radius','node_sphere_color_index','VMD_color_index_range','VMD_resolution','VMD_spline_smoothness']

	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''
	
	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
        parameters['calc_contact_map_boolean'] = False
        parameters['summary_boolean'] = False 
        parameters['nonstandard_substrates_selection'] = None
        parameters['homemade_selections'] = None
        parameters['alignment_selection'] = 'protein and name CA'
        parameters['traj_list'] = None
        parameters['trajectory_step'] = 1
        parameters['contact_map_distance_cutoff'] = 99999.9
        parameters['which_contact_map'] = 'average contact map'
        parameters['user_input_cartesian_covariance_matrix'] = None
        parameters['user_input_adjacency_matrix'] = None
        parameters['user_input_contact_map'] = None
        parameters['node_sphere_radius'] = 0.0
        parameters['node_sphere_rgb'] = (1.0,1.0,1.0)
        parameters['shortest_path_rgb'] = (1.0,1.0,1.0)
        parameters['longest_path_rgb'] = (1.0,1.0,1.0)
        parameters['node_sphere_color_index'] = 34
        parameters['VMD_color_index_range'] = (35,1056)
        parameters['VMD_resolution'] = 10
        parameters['VMD_spline_smoothness'] = 1

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
			sys.exit()

def summary(summary_file_name,arguments,parameters):
        """ Function to create a text file that holds important information about the analysis that was just performed. Outputs the version of MDAnalysis, how to rerun the analysis, and the parameters used in the analysis.

        Usage:
            summary(summary_file_name,arguments,parameters)

        Arguments:
            summary_file_name: string object of the file name to be written that holds the summary information.

        """
	with open(summary_file_name,'w') as f:
		f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
		f.write('\nAtom selections analyzed have been written out to node_selections.txt\n')
		f.write('To recreate this analysis, run this line:\n')
		for i in range(len(arguments)):
			f.write('%s ' %(arguments[i]))
		f.write('\n\n')
		f.write('Parameters used:\n')
                for key, value in parameters.iteritems():
                        if key == '__builtins__':
                                continue
                        if type(value) == int or type(value) == float:
			        f.write("%s = %s\n" %(key,value))
                        else:
			        f.write("%s = '%s'\n" %(key,value))



