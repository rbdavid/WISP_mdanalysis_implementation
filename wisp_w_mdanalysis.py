
# ----------------------------------------
# USAGE AND REFERENCING
# ----------------------------------------
#   python code_name.py config_file_name
#
# ----------------------------------------
# CODE OUTLINE
# ----------------------------------------
#   Allosteric Paths in Proteins
#       1) Load in user defined parameters
#       2) Load in necessary functions from module files
#       3) Create the MDAnalysis.Universe object and desired atom selections
#       4) Analyze trajectories (or the user can read in the necessary matrices)
#           a) Fill a numpy array with node positions
#           b) Calculate the average node positions using iterative alignment
#               i) Calculate RMSD from average node positions to find a traj. frame that is most representative of the average structure
#           c) Calculate the node cartesian covariance matrix
#           d) Calculate the node pair distances (if desired)
#       5) Calculate the desired adjacency matrix
#       6) Functionalize the adjacency matrix into a cost matrix
#       7) Run the desired path-finding algorithm on the cost matrix
#       8) Calculate desired edge and node metrics. 
#       9) Create VMD visualization states
#           a) Adjacency Matrix
#           b) Cost Matrix
#           c) Paths
#           d) Edge Metrics
#           e) Node Metrics
#           f) Node sampling around the average
#      10) Plot the various results
#           a) Adjacency Matrix
#           b) Cost Matrix
#           c) Edge Metrics
#           d) Node Metrics
#           e) Node pair distances; contact map
#      11) Output summary information to document the analyses performed
#
# ----------------------------------------
# PREAMBLE:
# ----------------------------------------

import sys
import os
import importlib
import numpy as np
import MDAnalysis
from IO import *

# ----------------------------------------
# VARIABLE DECLARATION: 
# ----------------------------------------

config_file = sys.argv[1]

# ----------------------------------------
# FUNCTIONS: 
# ----------------------------------------

def main():

#        # ----------------------------------------
#        # FILE NAMING VARIABLES
#        # ----------------------------------------
#        simply_formatted_paths_file_name = parameters['output_directory'] + 'simply_formatted_paths.txt'
#        pathways_vis_state_file_name = parameters['output_directory'] + 'pathways_vis_state.vmd'
#        summary_file_name = parameters['output_directory'] + 'summary.txt'

        # ----------------------------------------
        # 3) CREATE THE MDAnalysis.Universe OBJECT AND DESIRED ATOM SELECTIONS
        # ----------------------------------------
        u = MDAnalysis.Universe(parameters['pdb'])
        selection_list, source_indices, sink_indices = make_selections(u,parameters['output_directory'] + 'node_selections.txt',parameters['substrate_node_definition'],parameters['substrate_selection_string'],parameters['source_selection_string_list'],parameters['sink_selection_string_list'], nonstandard_substrates_selection = parameters['nonstandard_substrates_selection'], homemade_selections = parameters['homemade_selections'])

        print 'Number of nodes:', len(selection_list), '\nsource node indices:', source_indices, '\nsink node indices:', sink_indices, '\nNode selections written out to', parameters['output_directory'] + 'node_selections.txt'

        # ----------------------------------------
        # 4) ANALYZE TRAJECTORIES
        # ----------------------------------------
        if parameters['trajectory_analysis_boolean']:
                print 'Beginning trajectory analysis.'
                
                # ----------------------------------------
                # 4a and 4b) FILL A NUMPY ARRAY WITH NODE POSITIONS, CALCULATE THE AVERAGE POSITIONS USING ITERATIVE ALIGNMENT
                # ----------------------------------------
                Node_trajectory, avg_Node_positions = traj_alignment_and_averaging(u,parameters['alignment_selection'],selection_list,parameters['substrate_node_definition'],parameters['traj_list'],parameters['output_directory'], step = parameters['trajectory_step'], convergence_threshold = 1E-5)
                
                # ----------------------------------------
                # 4c) CALCULATE THE NODE CARTESIAN COVARIANCE MATRIX 
                # ----------------------------------------
                Node_cart_covariance = cartesian_covariance(Node_trajectory,avg_Node_positions,parameters['output_directory'])

                # ----------------------------------------
                # 4d) CALCULATE THE NODE PAIR DISTANCES
                # ----------------------------------------
                print 'Beginning node pair distance calculation.'
                binary_contact_map, avg_node_node_distances = calc_contact_map(Node_trajectory,parameters['output_directory'],distance_cutoff = parameters['contact_map_distance_cutoff'])
                if parameters['which_contact_map'].upper() == 'AVERAGE CONTACT MAP':
                        contact_map = avg_node_node_distances
                elif parameters['which_contact_map'].upper() == 'BINARY CONTACT MAP':
                        contact_map = binary_contact_map
                else:
                        print "User has not defined which contact map should be used in subsequent weighting of the adjacency matrix. Acceptable values for the 'which_contact_map' parameter are 'AVERAGE CONTACT MAP' or 'BINARY CONTACT MAP'."
                        sys.exit()
        
        else:
                print 'Loading in the user specified cartesian covariance data file. Note: no node trajectory has been created.'
                Node_cart_covariance = np.loadtxt(parameters['user_input_cartesian_covariance_matrix'])
                if Node_cart_covariance.shape[0]%3 != 0 or Node_cart_covariance.shape[0] != Node_cart_covariance.shape[1]:
                        print 'User has read in a cartesian covariance matrix that does not have the correct shape (expected: %d x %d, got: %d x %d).'%(3*len(selection_list),3*len(selection_list),Node_cart_covariance.shape[0],Node_cart_covariance.shape[1])
                        sys.exit()
                
                print 'Loading in the user specified average node positions file.'
                avg_Node_positions = np.loadtxt(parameters['user_input_average_node_positions'])
                if avg_Node_positions.shape != (len(selection_list),3):
                        print 'User has read in an average node positions file that does not have the correct shape (expected: %d x %d, got: %d x %d).'%(len(selection_list),3,avg_Node_positions.shape[0],avg_Node_positions.shape[1])
                        sys.exit()
                
                print 'Loading in the user specified contact map data file.'
                contact_map = np.loadtxt(parameters['user_input_contact_map'])  # can be either binary or average distances
                if contact_map.shape != (len(selection_list),len(selection_list)):
                        print 'User has read in a contact map that does not have the correct shape (expected: %d x %d, got: %d x %d).'%(len(selection_list),len(selection_list),contact_map.shape[0],contact_map.shape[1])
                        sys.exit()

        # ----------------------------------------
        # 5) CALCULATE THE DESIRED ADJACENCY MATRIX 
        # ----------------------------------------
        if parameters['adjacency_matrix_analysis_boolean']:
                print 'Beginning the calculation of the adjacency matrix (' + parameters['adjacency_matrix_style'] + ').'
                adjacency_matrix = adjacency_matrix_analysis(Node_cart_covariance,parameters['output_directory'])
        else:
                print 'Loading in the user specified adjacency matrix data files.'
                adjacency_matrix = np.loadtxt(parameters['user_input_adjacency_matrix'])
                if adjacency_matrix.shape != (len(selection_list),len(selection_list)):
                        print 'User has read in an adjacency matrix that does not have the correct shape (expected: %d x %d, got: %d x %d).'%(len(selection_list),len(selection_list),adjacency_matrix.shape[0],adjacency_matrix.shape[1])
                        sys.exit()

        # ----------------------------------------
        # 6) FUNCTIONALIZE THE ADJACENCY MATRIX INTO A COST MATRIX
        # ----------------------------------------
        if parameters['weight_by_contact_map_boolean']:
                functionalized_adjacency_matrix = functionalize_adjacency_matrix(adjacency_matrix,parameters['output_directory'], contact_map = contact_map, Lambda = parameters['lambda'])
        else:
                functionalized_adjacency_matrix = functionalize_adjacency_matrix(adjacency_matrix,parameters['output_directory'])
       
        print 'Finished functionalizing and weighting the adjacency matrix.'
        
        # ----------------------------------------
        # 7) RUN THE DESIRED PATH-FINDING ALGORITHM ON THE COST MATRIX
        # ----------------------------------------
        
        os.mkdir(parameters['output_directory'] + parameters[''] + os.sep)
        
        paths = get_paths(functionalized_adjacency_matrix,source_indices,sink_indices,parameters['number_of_paths'])
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
        # 8) CALCULATE DESIRED EDGE AND NODE METRICS
        # ----------------------------------------

        # ----------------------------------------
        # 9) CREATE VMD VISUALIZATION STATES
        # ----------------------------------------
        
        u.load_new(parameters['visualization_frame_pdb'])

        create_vis_state(parameters['visualization_frame_pdb'],selection_list,paths,pathways_vis_state_file_name, node_sphere_radius = parameters['node_sphere_radius'], node_sphere_rgb = parameters['node_sphere_rgb'], shortest_path_radius = parameters['shortest_path_radius'], shortest_path_rgb = parameters['shortest_path_rgb'], longest_path_radius = parameters['longest_path_radius'], longest_path_rgb = parameters['longest_path_rgb'], node_sphere_color_index = parameters['node_sphere_color_index'], VMD_color_index_range = parameters['VMD_color_index_range'], VMD_resolution = parameters['VMD_resolution'], VMD_spline_smoothness = parameters['VMD_spline_smoothness'])
        print 'Finished creating the vis state file.'

        # ----------------------------------------
        # 10) PLOT THE VARIOUS RESULTS
        # ----------------------------------------
        #if parameters['plotting_boolean'] and len(paths) != 1:
        #        length_frequency_plot = parameters['output_directory'] + 'path_length_frequency.png'
        #        node_frequency_plot = parameters['output_directory'] + 'node_frequency.png'
        #        paths_analysis_plotting(paths,source_indices,sink_indices,len(func_corr_array),length_frequency_plot,node_frequency_plot)

        # ----------------------------------------
        # PLOTTING DATA - CORRELATION AND CONTACT MAPS
        # ----------------------------------------
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
        # 11) OUTPUT SUMMARY INFORMATION TO DOCUMENT THE ANALYSES PERFORMED
        # ----------------------------------------
        if parameters['summary_boolean']:
                summary(summary_file_name,sys.argv,parameters)

# ----------------------------------------
# 1) LOAD IN USER DEFINED PARAMETERS
# ----------------------------------------
parameters = {}
config_parser(config_file,parameters)

# ----------------------------------------
# SETTING UP THE OUTPUT DIRECTORY
# ----------------------------------------
if parameters['output_directory'][-1] != os.sep:
        parameters['output_directory'] += os.sep

if os.path.exists(parameters['output_directory']):
        print 'The output directory, ', parameters['output_directory'], 'already exists. Please select a different directory name for output.'
        sys.exit()
else:
        os.mkdir(parameters['output_directory'])

# ----------------------------------------
# 2) LOAD IN NECESSARY FUNCTIONS FROM MODULE FILES
# ----------------------------------------

###
make_selections = importlib.import_module(parameters['node_selection_file'].split('.')[0],package=None).make_selections

###
if parameters['trajectory_analysis_boolean']:
        traj_alignment_and_averaging = importlib.import_module(parameters['trajectory_functions_file'].split('.')[0],package=None).traj_alignment_and_averaging
        cartesian_covariance = importlib.import_module(parameters['trajectory_functions_file'].split('.')[0],package=None).cartesian_covariance
        calc_contact_map = importlib.import_module(parameters['trajectory_functions_file'].split('.')[0],package=None).calc_contact_map

###
if parameters['adjacency_matrix_style'].upper() in ('PEARSON','PEARSON CORRELATION'):
        adjacency_matrix_analysis = importlib.import_module(parameters['adjacency_matrix_functions_file'].split('.')[0],package=None).pearson_correlation_analysis
        functionalize_adjacency_matrix = importlib.import_module(parameters['func_adjacency_matrix_functions_file'].split('.')[0],package=None).func_pearson_correlation

elif parameters['adjacency_matrix_style'].upper() in ('LMI','LINEAR MUTUAL INFORMATION'):
        adjacency_matrix_analysis = importlib.import_module(parameters['adjacency_matrix_functions_file'].split('.')[0],package=None).linear_mutual_information_analysis
        functionalize_adjacency_matrix = importlib.import_module(parameters['func_adjacency_matrix_functions_file'].split('.')[0],package=None).generalized_correlation_coefficient_calc

elif parameters['adjacency_matrix_style'].upper() in ('HENM','HENM HESSIAN'):
        adjacency_matrix_analysis = importlib.import_module(parameters['adjacency_matrix_functions_file'].split('.')[0],package=None).hENM_hessian_analysis
        functionalize_adjacency_matrix = importlib.import_module(parameters['func_adjacency_matrix_functions_file'].split('.')[0],package=None).func_hessian

elif parameters['adjacency_matrix_style'].upper() in ('REACH','REACH HESSIAN'):
        adjacency_matrix_analysis = importlib.import_module(parameters['adjacency_matrix_functions_file'].split('.')[0],package=None).REACH_hessian_analysis
        functionalize_adjacency_matrix = importlib.import_module(parameters['func_adjacency_matrix_functions_file'].split('.')[0],package=None).func_hessian

else:
        print "The user has not read in an accepted style of adjacency matrix. Current options are pearson correlation coefficient matrix (denoted by 'PEARSON' or 'PEARSON CORRELATION'), linear mutual information (denoted by 'LMI' or 'LINEAR MUTUAL INFORMATION'), hetero elastic network model hessian (denoted by 'HENM' or 'HENM HESSIAN'), or REACH hessian (denoted by COVARIANCE HESSIAN or REACH), or ."
        sys.exit()
#functionalize_adjacency_matrix = importlib.import_module(parameters['func_adjacency_matrix_functions_file'].split('.')[0],package=None).do_nothing

###
get_paths = importlib.import_module(parameters['network_functions_file'].split('.')[0],package=None).get_paths

###
create_vis_state = importlib.import_module(parameters['visualization_functions_file'].split('.')[0],package=None).create_vis_state

#if parameters['plotting_boolean']:
#        plot_square_matrix = importlib.import_module(parameters['plotting_functions_file'].split('.')[0],package=None).plot_square_matrix
#        paths_analysis_plotting = importlib.import_module(parameters['plotting_functions_file'].split('.')[0],package=None).paths_analysis_plotting

# ----------------------------------------
# MAIN
# ----------------------------------------
if __name__ == '__main__':
	main()

