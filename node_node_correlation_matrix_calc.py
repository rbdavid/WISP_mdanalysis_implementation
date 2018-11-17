#!/home/rbdavid/bin/python
# ----------------------------------------
# USAGE:

# ./rmsd.ref.py config_file

# ----------------------------------------
# PREAMBLE:

import sys
import importlib
import numpy as np
import MDAnalysis

config_file = sys.argv[1]

# ----------------------------------------
# FUNCTIONS: 

necessary_parameters = ['average_pdb','pdb','traj_list','data_output_filename','selection_file','distance_functions_file']
all_parameters = ['average_pdb','pdb','traj_list','data_output_filename','selection_file','distance_functions_file','alignment_selection','substrates_selection','homemade_selections','functionalize_distance_correlation_boolean','summary_boolean','summary_filename','selection_output_filename']
def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['alignment_selection'] = 'protein'
	parameters['substrates_selection'] = None
	parameters['homemade_selections'] = None
	parameters['summary_boolean'] = False 
	parameters['summary_filename'] = None
	parameters['selection_output_filename'] = 'selections.txt'

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
			sys.exit()

def summary(summary_filename):
	with open(summary_filename,'w') as f:
		f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
		f.write('To recreate this analysis, run this line:\n')
		for i in range(len(sys.argv)):
			f.write('%s ' %(sys.argv[i]))
		f.write('\n\n')
		f.write('Parameters used:\n')
                for key, value in parameters.iteritems():
                        if key == '__builtins__':
                                continue
                        else:
			        f.write("%s = '%s'\n" %(key,value))
		f.write('\n')

def main():

        # ----------------------------------------
        # LOAD IN THE ANALYSIS UNIVERSE
        u = MDAnalysis.Universe(parameters['pdb'])
        
        # ----------------------------------------
        # MAKE ATOM SELECTIONS FOR SUBSEQUENT COVARIANCE ANALYSIS
        u_temp_selection = u.select_atoms(selections[0])
        nResidues = u_temp_selection.n_residues

        selection_list = []
        with open(parameters['selection_output_filename'],'w') as f:
                for i in range(nResidues):
                        u_temp = u_temp_selection.residues[i].atoms
                        selection_list.append(u_temp)
                        f.write("%02d   %s   '%s'\n" %(i,u_temp.resnames[0],u_temp.resids[0]))

                count = len(selection_list)
                if parameters['substrates_selection'] != None:
                        u_subs = u.select_atoms(parameters['substrates_selection'])
        		for i in range(u_subs.n_residues):
        			temp_resname = u_subs.residues[i].resname
        			temp_id = u_subs.residues[i].resid
        			if temp_resname in parameters['homemade_selections']:
        				make_selections(u,temp_resname,temp_id,f,selection_list,count)
                                        count = len(selection_list)
        			else:
        				u_temp = u.select_atoms('resname %s and resid %d' %(temp_resname,temp_id))
        				selection_list.append(u_temp)
        				f.write('%02d   %s   resid %d\n' %(count,temp_resname,temp_id))
                                        count = len(selection_list)
       
        nNodes = len(selection_list)
        nNodes_range = range(nNodes)

        # ----------------------------------------
        # LOAD IN THE AVERAGE UNIVERSE
        avg = MDAnalysis.Universe(parameters['average_pdb'])
        avg_align = u.select_atoms(parameters['alignment_selection'])
        avg_all = u.select_atoms('all')
        avg_all.translate(-avg_align.center_of_mass())
        pos0 = avg_align.positions

        # ----------------------------------------
        # CALCULATING THE CARTESIAN COVARIANCE, AVERAGE, AND VARIANCE ARRAYS
        if selections[1] == 'COM':
                print 'Performing a covariance analysis of the cartesian coordinates of the center of mass of residues defined in the selection_file.'
                covariance_array, average_array, variance_array = calc_cart_covar_matrix(u,parameters['traj_list'],parameters['alignment_selection'],pos0,selection_list,nNodes)
               
                np.savetxt('cart_covar.' + parameters['data_output_filename'],covariance_array)
                np.savetxt('cart_var.' + parameters['data_output_filename'],variance_array)
                np.savetxt('cart_avg.' + parameters['data_output_filename'],average_array)

        else:
                print 'Requested to calc the covariance array of something other than the center of mass of residues defined in the selection_file... Currently, this code is unable to do anything other than COM'
                sys.exit()

        # ----------------------------------------
        # CALCULATING THE DISTANCE CORRELATION MATRIX OF NODE-NODE PAIRS
        distance_correlation_matrix = np.zeros((nNodes,nNodes),dtype=np.float64)
        distance_variance_array = np.zeros((nNodes),dtype=np.float64)
        for i in nNodes_range:
                dim1 = i*3
                distance_variance_array[i] = sum(variance_array[dim1:dim1+3])   # summing the variances of the cartesian dimensions of the node i; <r_{i}(t)**2>
                for j in nNodes_range[i:]:
                        dim2 = j*3
                        distance_correlation_matrix[i,j] = covariance_array[dim1,dim2] + covariance_array[dim1+1,dim2+1] + covariance_array[dim1+2,dim2+2]      # taking the trace of the cartesian covariance array for the desired dimensions of nodes i and j; <r_{i}(t) dot r_{j}(t)> - <r_{i}(t)> dot <r_{j}(t)> 

        for i in nNodes_range:
                for j in nNodes_range[i:]:
                        distance_correlation_matrix[i,j] /= np.sqrt(distance_variance_array[i]*distance_variance_array[j])     # finishing the distance correlation matrix by dividing by the product of the variances of nodes i,j;  <r_{i}(t) dot r_{j}(t)> - <r_{i}(t)> dot <r_{j}(t)> / sqrt((<x_{i}(t)**2> - <x_{i}(t)>**2)*(<x_{j}(t)**2> - <x_{j}(t)>**2))
                        distance_correlation_matrix[j,i] = distance_correlation_matrix[i,j] # filling in the bottom traingle of this matrix

        np.savetxt('dist_corr.' + parameters['data_output_filename'],distance_correlation_matrix)
        np.savetxt('dist_var.' + parameters['data_output_filename'],distance_variance_array)

        # ----------------------------------------
        # FUNCTIONALIZE THE DISTANCE CORRELATION MATRIX (USING WISP EQUATION)
        if parameters['functionalize_distance_correlation_boolean']:
                func_distance_correlation_matrix = np.zeros((nNodes,nNodes),dtype=np.float64)
                print 'Beginning to functionalize the distance correlation matrix.'
                for i in nNodes_range:
                        for j in nNodes_range[i:]:
                                func_distance_correlation_matrix[i,j] = -np.log(np.fabs(distance_correlation_matrix[i,j]))
                                func_distance_correlation_matrix[j,i] = func_distance_correlation_matrix[i,j]
                np.savetxt('func_dist_corr.' + parameters['data_output_filename'],func_distance_correlation_matrix)
        else:
                print 'No functionalization of the distance correlation matrix is being performed.'

        # ----------------------------------------
        # OUTPUTTING SUMMARY INFORMATION
        if parameters['summary_boolean']:
                summary(parameters['summary_filename'])

# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

# ----------------------------------------
# LOADING IN NECESSARY FUNCTIONS FROM MODULE FILES
make_selections = importlib.import_module(parameters['selection_file'].split('.')[0],package=None).make_selections
selections = importlib.import_module(parameters['selection_file'].split('.')[0],package=None).selections
calc_cart_covar_matrix = importlib.import_module(parameters['distance_functions_file'].split('.')[0],package=None).calc_cart_covar_matrix

# ----------------------------------------
# MAIN
if __name__ == '__main__':
	main()

