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

necessary_parameters = ['output_directory','node_selection_file','trajectory_functions_file','adjacency_matrix_functions_file','func_adjacency_matrix_functions_file','network_functions_file','visualization_functions_file','trajectory_analysis_boolean','adjacency_matrix_style','pdb','visualization_frame_pdb','substrate_node_definition','substrate_selection_string','source_selection_string_list','sink_selection_string_list','number_of_paths']

all_parameters = ['output_directory','node_selection_file','trajectory_functions_file','adjacency_matrix_functions_file','func_adjacency_matrix_functions_file','network_functions_file','visualization_functions_file','trajectory_analysis_boolean','adjacency_matrix_style','pdb','visualization_frame_pdb','substrate_node_definition','substrate_selection_string','source_selection_string_list','sink_selection_string_list','number_of_paths','weight_by_contact_map_boolean','summary_boolean','nonstandard_substrates_selection','homemade_selections','alignment_selection','traj_list','trajectory_step','contact_map_distance_cutoff','which_contact_map','user_input_cartesian_covariance_matrix','user_input_contact_map','node_sphere_radius','node_sphere_rgb','shortest_path_rgb','longest_path_rgb','shortest_path_radius','longest_path_radius','node_sphere_color_index','VMD_color_index_range','VMD_resolution','VMD_spline_smoothness']

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
        parameters['weight_by_contact_map_boolean'] = False
        parameters['summary_boolean'] = False 
        parameters['nonstandard_substrates_selection'] = None
        parameters['homemade_selections'] = None
        parameters['alignment_selection'] = 'protein and name CA'
        parameters['traj_list'] = None
        parameters['trajectory_step'] = 1
        parameters['contact_map_distance_cutoff'] = 99999.9
        parameters['which_contact_map'] = 'average contact map'
        parameters['user_input_cartesian_covariance_matrix'] = None
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
        # FILE NAMING VARIABLES
        # ----------------------------------------
        node_selection_file_name = parameters['output_directory'] + 'node_selections.txt'
        simply_formatted_paths_file_name = parameters['output_directory'] + 'simply_formatted_paths.txt'
        pathways_vis_state_file_name = parameters['output_directory'] + 'pathways_vis_state.vmd'
        summary_file_name = parameters['output_directory'] + 'summary.txt'

        # ----------------------------------------
        # LOAD IN THE ANALYSIS UNIVERSE
        # ----------------------------------------
        u = MDAnalysis.Universe(parameters['pdb'])
        
        # ----------------------------------------
        # MAKE ATOM SELECTIONS FOR USER-SPECIFIED NODE DEFINITIONS
        # ----------------------------------------
        selection_list, source_indices, sink_indices = make_selections(u,node_selection_file_name,parameters['substrate_node_definition'],parameters['substrate_selection_string'],parameters['source_selection_string_list'],parameters['sink_selection_string_list'], nonstandard_substrates_selection = parameters['nonstandard_substrates_selection'], homemade_selections = parameters['homemade_selections'])
        print 'source node indices:', source_indices, '     sink node indices:', sink_indices

        # ----------------------------------------
        # TRAJECTORY ANALYSIS
        # ----------------------------------------
        if parameters['trajectory_analysis_boolean']:
                
                # ----------------------------------------
                # ALIGN THE NODE TRAJECTORY TO THE AVERAGE NODE POSITIONS; ITERATIVE AVERAGE METHOD
                # ----------------------------------------
                print 'Beginning trajectory analysis.'
                Node_trajectory, avg_Node_positions, Node_cart_covariance = alignment_averaging_and_covariance_analysis(u,parameters['alignment_selection'],selection_list,parameters['traj_list'],parameters['output_directory'], step = parameters['trajectory_step'], convergence_threshold = 1E-5)

                # ----------------------------------------
                # CONTACT MAP CALCULATION
                # ----------------------------------------
                if parameters['weight_by_contact_map_boolean']:
                        print 'Beginning node pair distance calculation.'
                        binary_contact_map, avg_node_node_distances = calc_contact_map(Node_trajectory,parameters['contact_map_distance_cutoff'], output_directory = parameters['output_directory'])
                        if parameters['which_contact_map'].upper() == 'AVERAGE CONTACT MAP':
                                contact_map = avg_node_node_distances
                        elif parameters['which_contact_map'].upper() == 'BINARY CONTACT MAP':
                                contact_map = binary_contact_map
                        else:
                                print "User has not defined which contact map should be used in subsequent weighting of the adjacency matrix. Acceptable values for the 'which_contact_map' parameter are 'AVERAGE CONTACT MAP' or 'BINARY CONTACT MAP'."
                                sys.exit()

        # ----------------------------------------
        # USER SPECIFIED DATA INPUT
        # ----------------------------------------
        else:
                print 'Loading in the user specified data files.'
                Node_cart_covariance = np.loadtxt(parameters['user_input_cartesian_covariance_matrix'])

                if parameters['weight_by_contact_map_boolean']:
                        contact_map = np.loadtxt(parameters['user_input_contact_map'])  # can be either binary or average distances
                        if Node_cart_covariance.shape != contact_map.shape:
                                print 'User has read in a contact map that does not have the same shape as the cartesian covariance data.'
                                sys.exit()

        # ----------------------------------------
        # CALC THE DESIRED ADJACENCY MATRIX 
        # ----------------------------------------
        adjacency_matrix = adjacency_matrix_analysis(Node_cart_covariance,parameters['output_directory'])
        print 'Finished calculating the adjacency matrix (' + parameters['adjacency_matrix_style'] + ').'

        # ----------------------------------------
        # FUNCTIONALIZE THE ADJACENCY MATRIX
        # ----------------------------------------
        if parameters['weight_by_contact_map_boolean']:
                funcionalized_adjacency_matrix = functionalize_adjacency_matrix(adjacency_matrix,parameters['output_directory'], contact_map = contact_map, Lambda = parameters['lambda'])
        else:
                funcionalized_adjacency_matrix = functionalize_adjacency_matrix(adjacency_matrix,parameters['output_directory'])
       
        print 'Finished functionalizing and weighting the adjacency matrix.'
        
        ## ----------------------------------------
        ## PLOTTING DATA - CORRELATION AND CONTACT MAPS
        ## ----------------------------------------
        #if parameters['plotting_boolean']:
        #        matrix_plot = parameters['output_directory'] + 'node_node_correlation_matrix.png'
        #        plot_square_matrix(correlation_array,matrix_plot,axes_label = 'Node Index', cbar_label = 'Correlation', v_range = [-1.0,1.0], minor_ticks = 5, major_ticks = 50)

        #        contact_map_plot = parameters['output_directory'] + 'node_node_contact_map.png'
        #        plot_square_matrix(contact_map,contact_map_plot,axes_label = 'Node Index', cbar_label = r'Distance ($\AA$)', minor_ticks = 5, major_ticks = 50)

        #        func_matrix_plot = parameters['output_directory'] + 'functionalized_node_node_correlation_matrix.png'
        #        if parameters['weight_by_avg_distance_boolean']:
        #                plot_square_matrix(func_corr_array,func_matrix_plot,axes_label = 'Node Index', cbar_label = r'Correlation Distance ($\AA$)', minor_ticks = 5, major_ticks = 50)
        #        else:
        #                plot_square_matrix(func_corr_array,func_matrix_plot,axes_label = 'Node Index', cbar_label = 'Correlation', minor_ticks = 5, major_ticks = 50)

        # ----------------------------------------
        # COMPUTE THE ENSEMBLE OF PATHS THAT CONNECT THE SOURCE AND SINK RESIDUES 
        # ----------------------------------------
        paths = get_paths(funcionalized_adjacency_matrix,source_indices,sink_indices,parameters['number_of_paths'])
        print 'Finished calculating the pathways that connect the source and sink selections.'
        
        with open(simply_formatted_paths_file_name,'w') as W:
                W.write('# Pathway_Length Node_Indices (zero indexed)\n')
                for path in paths:
                        path_string = '%f ' %(path[0])
                        for node in path[1:]:
                                path_string += '%d ' %(node)
                        path_string += '\n'

                        W.write(path_string)

        # ----------------------------------------
        # CREATE THE VISUALIZATION STATES TO BE USED IN VMD 
        # ----------------------------------------
        u.load_new(parameters['visualization_frame_pdb'])

        create_vis_state(parameters['visualization_frame_pdb'],selection_list,paths,pathways_vis_state_file_name, node_sphere_radius = parameters['node_sphere_radius'], node_sphere_rgb = parameters['node_sphere_rgb'], shortest_path_radius = parameters['shortest_path_radius'], shortest_path_rgb = parameters['shortest_path_rgb'], longest_path_radius = parameters['longest_path_radius'], longest_path_rgb = parameters['longest_path_rgb'], node_sphere_color_index = parameters['node_sphere_color_index'], VMD_color_index_range = parameters['VMD_color_index_range'], VMD_resolution = parameters['VMD_resolution'], VMD_spline_smoothness = parameters['VMD_spline_smoothness'])
        print 'Finished creating the vis state file.'

        ## ----------------------------------------
        ## PLOTTING DATA - PATHWAYS SAMPLING BAR PLOTS
        ## ----------------------------------------
        #if parameters['plotting_boolean'] and len(paths) != 1:
        #        length_frequency_plot = parameters['output_directory'] + 'path_length_frequency.png'
        #        node_frequency_plot = parameters['output_directory'] + 'node_frequency.png'
        #        paths_analysis_plotting(paths,source_indices,sink_indices,len(func_corr_array),length_frequency_plot,node_frequency_plot)

        # ----------------------------------------
        # OUTPUTTING SUMMARY INFORMATION
        # ----------------------------------------
        if parameters['summary_boolean']:
                summary(summary_file_name)

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
        alignment_averaging_and_covariance_analysis = importlib.import_module(parameters['trajectory_functions_file'].split('.')[0],package=None).alignment_averaging_and_covariance_analysis
        if parameters['weight_by_contact_map_boolean']:
                calc_contact_map = importlib.import_module(parameters['trajectory_functions_file'].split('.')[0],package=None).calc_contact_map

# which adjacency matrix?
if parameters['adjacency_matrix_style'].upper() == 'PEARSON CORRELATION':
        adjacency_matrix_analysis = importlib.import_module(parameters['adjacency_matrix_functions_file'].split('.')[0],package=None).pearson_correlation_analysis
        functionalize_adjacency_matrix = importlib.import_module(parameters['func_adjacency_matrix_functions_file'].split('.')[0],package=None).func_pearson_correlation

elif parameters['adjacency_matrix_style'].upper() in ('LMI','LINEAR MUTUAL INFORMATION'):
        adjacency_matrix_analysis = importlib.import_module(parameters['adjacency_matrix_functions_file'].split('.')[0],package=None).linear_mutual_information_analysis
        functionalize_adjacency_matrix = importlib.import_module(parameters['func_adjacency_matrix_functions_file'].split('.')[0],package=None).generalized_correlation_coefficient_calc

elif parameters['adjacency_matrix_style'].upper() == 'MUTUAL INFORMATION':
        adjacency_matrix_analysis = importlib.import_module(parameters['adjacency_matrix_functions_file'].split('.')[0],package=None).mutual_information_analysis
        functionalize_adjacency_matrix = importlib.import_module(parameters['func_adjacency_matrix_functions_file'].split('.')[0],package=None).generalized_correlation_coefficient_calc

elif parameters['adjacency_matrix_style'].upper() in ('REACH','COVARIANCE HESSIAN'):
        adjacency_matrix_analysis = importlib.import_module(parameters['adjacency_matrix_functions_file'].split('.')[0],package=None).REACH_analysis
        functionalize_adjacency_matrix = importlib.import_module(parameters['func_adjacency_matrix_functions_file'].split('.')[0],package=None).func_k_matrix

elif parameters['adjacency_matrix_style'].upper() == 'HENM':
        adjacency_matrix_analysis = importlib.import_module(parameters['adjacency_matrix_functions_file'].split('.')[0],package=None).hENM_analysis
        functionalize_adjacency_matrix = importlib.import_module(parameters['func_adjacency_matrix_functions_file'].split('.')[0],package=None).func_k_matrix

else:
        print 'The user has not read in an accepted style of adjacency matrix. Current options are pearson correlation coefficient matrix (denoted by PEARSON CORRELATION), linear mutual information (denoted by LMI or LINEAR MUTUAL INFORMATION), covariance hessian or reach (denoted by COVARIANCE HESSIAN or REACH), or hetero elastic network model (denoted by HENM).'
        sys.exit()

get_paths = importlib.import_module(parameters['network_functions_file'].split('.')[0],package=None).get_paths

create_vis_state = importlib.import_module(parameters['visualization_functions_file'].split('.')[0],package=None).create_vis_state

#if parameters['plotting_boolean']:
#        plot_square_matrix = importlib.import_module(parameters['plotting_functions_file'].split('.')[0],package=None).plot_square_matrix
#        paths_analysis_plotting = importlib.import_module(parameters['plotting_functions_file'].split('.')[0],package=None).paths_analysis_plotting

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

