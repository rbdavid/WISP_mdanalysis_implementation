#!
#!/home/rbdavid/bin/python

# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import MDAnalysis
from MDAnalysis.analysis.align import rotation_matrix
from distance_functions import *
from selection_list import *
from homemade_atom_selections import *

# ----------------------------------------
# VARIABLE DECLARATION: 
# ----------------------------------------

config_file = sys.argv[1]

zeros = np.zeros
sqrt = np.sqrt
sum = np.sum 
dot = np.dot
eigen = np.linalg.eig

necessary_parameters = ['ref_pdb','pdb','traj_loc','start','end','Wrapped','outputname','selection_file','resid_offset']
all_parameters = ['ref_pdb','pdb','traj_loc','start','end','Wrapped','outputname','selection_file','resid_offset','alignment','wrapping_selection','substrates','homemade_selections','write_summary','summary_filename','selection_output']

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
	parameters['alignment'] = 'protein'
	parameters['wrapping_selection'] = 'not (resname WAT or resname Na+ or resname Cl- or protein)'
	parameters['substrates'] = None
	parameters['homemade_selections'] = None
	parameters['write_summary'] = False 
	parameters['summary_filename'] = 'rmsd.summary'
	parameters['selection_output'] = 'selections.txt'
	
	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
			sys.exit()

	### ADD VARIABLE TESTS...

def summary(summary_filename):
        """ Function to create a text file that holds important information about the analysis that was just performed. Outputs the version of MDAnalysis, how to rerun the analysis, and the parameters used in the analysis.

        Usage:
            summary(summary_filename)

        Arguments:
            summary_filename: string object of the file name to be written that holds the summary information.

        """
	with open(summary_filename,'w') as f:
		f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
		f.write('\nAtom selections analyzed have been written out to %s\n' %(parameters['selection_output']))
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

def make_selection(u,node_type,user_defined_selection,selection_output_filename):
	selection_list = []

	if node_type == 'atomic':
		selection = u.select_atoms(user_defined_selection)	# atomic selection string; most dangerous selection string... 
		nAtoms = selection.n_atoms 
		temp = '%0'+'%d'%(int(np.log10(nAtoms))+1)+'d'	# for correct zero pad when numbering selections
		with open(selection_output_filename,'w') as W:
			W.write('# selection_num   Atom_num   Atom_name   Resid\n')
			for i in range(nAtoms):
				selection_list.append(selection.atoms[i])
				W.write(temp%(i)+'   %d   %s   %d\n' %(selection.atoms[i].number,selection.atoms[i].name,selection.atoms[i].resid))
		return selection, selection_list, nAtoms

	elif node_type == 'residue_COM':
		selection = u.select_atoms(user_defined_selection)	# does not necessitate only protein selection; assumes only one node per residue
		nResidues = selection.n_residues
		order_of_magnitude = int(np.log10(nResidues))+1		# for correct zero pad when numbering selections
		temp = '%0'+'%d'%(order_of_magnitude)+'d'
		with open(selection_output_filename,'w') as W:
			W.write('# selection_num   Res_name   Resid   Num_atoms ! CENTER OF MASS SELECTION\n')
			for i in range(nResidues):
				selection_list.append(selection.residues[i])	# want to select all atoms in residue i of selection;
				W.write(temp%(i)+'   %s   %d   %d\n'%(selection.residues[i].resname,selection.residues[i].resid,selection.residues[i].n_atoms))
		return selection, selection_list, nResidues

	elif node_type == 'backbone_COM':
		selection = u.select_atoms('(backbone or nucleicbackbone) and '+user_defined_selection)	# does not necessitate only protein selection; assumes only one node per residue.
		nResidues = selection.n_residues
		order_of_magnitude = int(np.log10(nResidues))+1		# for correct zero pad when numbering selections
		temp = '%0'+'%d'%(order_of_magnitude)+'d'
		with open(selection_output_filename,'w') as W:
			W.write('# selection_num   Res_name   Resid   Num_atoms ! CENTER OF MASS OF BACKBONE SELECTION\n')
			for i in range(nResidues):
				temp_atom_group = selection[selection.resids == selection.residues[i].resid]	# only want to include atoms in that SHOULD BE in the atom selection. 
				selection_list.append(temp_atom_group)
				W.write(temp%(i)+'   %s   %d   %d\n'%(temp_atom_group.atoms[0].resname, temp_atom_group.atoms[0].resid, temp_atom_group.n_atoms))
		return selection, selection_list, nResidues

	elif node_type == 'sidechain_COM':
		selection = u.select_atoms('not (backbone or nucleicbackbone) and '+user_defined_selection)	# does not necessitate only protein selection; assumes only one node per residue.
		nResidues = selection.n_residues
		order_of_magnitude = int(np.log10(nResidues))+1		# for correct zero pad when numbering selections
		temp = '%0'+'%d'%(order_of_magnitude)+'d'
		with open(selection_output_filename,'w') as W:
			W.write('# selection_num   Res_name   Resid   Num_atoms ! CENTER OF MASS OF SIDECHAIN SELECTION\n')
			for i in range(nResidues):
				temp_atom_group = selection[selection.resids == selection.residues[i].resid]	# only want to include atoms in that SHOULD BE in the atom selection. 
				selection_list.append(temp_atom_group)
				W.write(temp%(i)+'   %s   %d   %d\n'%(temp_atom_group.atoms[0].resname, temp_atom_group.atoms[0].resid, temp_atom_group.n_atoms))
		return selection, selection_list, nResidues
		
	### STILL TO DO:
	# MORE COMPLEX, USER DEFINED ATOM SELECTIONS; good example: multiple sites per residue (i.e. backbone and sidechain COM for each residue of interest, with the correct node ordering based on resid)

	else:
		print "Parameter value for 'node_type' does not equal an expected string. Please check the parameter and try again."
		sys.exit()

def get_average_coordinates(avg_universe,node_type,alignment_selection_string,user_defined_selection,nNodes):
	"""
	"""
	if PARAMETERS['node_type'] == 'atomic':
		node_average_coordinates = avg_universe.select_atoms(user_defined_selection).positions
		align_average_coordinates = avg_universe.select_atoms(alignment_selection_string).positions
		return node_average_coordinates, align_average_coordinates
	elif node_type == 'residue_COM':
		align_average_coordinates = avg_universe.select_atoms(alignment_selection_string).positions
		selection = avg_universe.select_atoms(user_defined_selection)	# does not necessitate only protein selection; assumes only one node per residue
		
		node_average_coordinates = np.zeros((nNodes,3),dtype=np.float32)
		for i in range(nNodes):
			node_average_coordinates[i] = selection.residues[i].center_of_mass()
		return node_average_coordinates, align_average_coordinates

	elif node_type == 'backbone_COM':
		align_average_coordinates = avg_universe.select_atoms(alignment_selection_string).positions
		selection = avg_universe.select_atoms('(backbone or nucleicbackbone) and '+user_defined_selection)	# does not necessitate only protein selection; assumes only one node per residue
		
		node_average_coordinates = np.zeros((nNodes,3),dtype=np.float32)
		for i in range(nNodes):
			node_average_coordinates[i] = selection[selection.resids == selection.residues[i].resid].center_of_mass()
		return node_average_coordinates, align_average_coordinates

	elif node_type == 'sidechain_COM':
		align_average_coordinates = avg_universe.select_atoms(alignment_selection_string).positions
		selection = avg_universe.select_atoms('not (backbone or nucleicbackbone) and '++user_defined_selection)	# does not necessitate only protein selection; assumes only one node per residue
		
		node_average_coordinates = np.zeros((nNodes,3),dtype=np.float32)
		for i in range(nNodes):
			node_average_coordinates[i] = selection[selection.resids == selection.residues[i].resid].center_of_mass()
		return node_average_coordinates, align_average_coordinates

	### STILL TO DO:
	# MORE COMPLEX, USER DEFINED ATOM SELECTIONS; good example: multiple sites per residue (i.e. backbone and sidechain COM for each residue of interest, with the correct node ordering based on resid)

def iterative_avg(u,trajectory_list,u_all,u_align,u_analysis_selection_list,node_type,nNodes,nAlign,nSteps,trajectory_stepping_value=1,threshold=1E-5,maxIter=100):
	"""
	"""
	# ----------------------------------------
	# INITIALIZING ARRAYS
	all_coord = zeros((nSteps,nNodes,3),dtype=np.float32)
	avg_coord = zeros((nNodes,3),dtype=np.float32)
	all_align = zeros((nSteps,nAlign,3),dtype=np.float32)
	avg_align = zeros((nAlign,3),dtype=np.float32)

	# ----------------------------------------
	# ANALYZING TRAJECTORY
	print 'Beginning trajectory analysis'
	temp = 0
	if node_type == 'atomic':
		for a in trajectory_list:
			u.load_new(a)
			for ts in u.trajectory[::trajectory_stepping_value]:
				u_all.translate(-u_align.center_of_mass())
				avgAlign += u_align.positions
				all_align[tem] = u_align.positions
				for i in range(len(u_analysis_selection_list)):
					avgCoord += u_analysis_selection_list[i].pos
					all_coord[temp] = u_analysis_selection_list[i].pos
				temp += 1
	elif 'COM' in node_type:
		for a in trajectory_list:
			u.load_new(a)
			for ts in u.trajectory[::trajectory_stepping_value]:
				u_all.translate(-u_align.center_of_mass())
				avgAlign += u_align.positions
				all_align[tem] = u_align.positions
				for i in range(len(u_analysis_selection_list)):
					avgCoord += u_analysis_selection_list[i].center_of_mass()
					all_coord[temp] = u_analysis_selection_list[i].center_of_mass()
				temp += 1
	
	avg_align /= nSteps
	avg_coord /= nSteps
	
	iteration = 0 
	residual = threshold + 100.
	print 'Beginning iterative process of calculating average coordinates and aligning to the new average.'
	while residual > threshold and iteration < maxIter:
		temp_avg_coord = np.zeros((nNodes,3),dtype=np.float32)
		temp_avg_align = np.zeros((nAlign,3),dtype=np.float32)
		for i in range(nSteps):
			R, d = rotation_matrix(all_align[i,:,:],avg_align)
			all_align[i,:,:] = dot(all_align[i,:,:],R.T)
			all_coord[i,:,:] = dot(all_coord[i,:,:],R.T)
			temp_avg_align += all_align[i,:,:]
			temp_avg_coord += all_coord[i,:,:]
		temp_avg_align /= nSteps
		temp_avg_coord /= nSteps
		residual = RMSD(avg_align, temp_avg_align, nAlign)
		rmsd_nodes = RMSD(avg_coord, temp_avg_coord, nNodes)
		avg_align = temp_avg_align
		avg_coord = temp_avg_coord
		print 'Finished with Iteration step', iteration, ', RMSD of alignment landmark atoms is', residual, ' RMSD of Node selections is', rmsd_nodes
		iteration += 1
	
	return avg_coord, avg_align
	
def main():

####code breakdown
# 1) calc iterative, average structure using an alignment landmark if required; only necessary for the atom selection to be analyzed
# 2) calc the distance between node(t) and <node> for every t; produces d matrix
# 3) calc the correlation between d(i) and d(j)
# 4) read in contact map...
# 5) 

	# ----------------------------------------
	# INITIALIZING MDANALYSIS UNIVERSES 
	u = MDAnalysis.Universe(parameters['pdb'])
	u_all = u.select_atoms('all')
	u_align = u.select_atoms(parameters['alignment'])
	u_analysis_selection, u_analysis_selection_list, nNodes = make_selection(u,PARAMETERS['node_type'],PARAMETERS['user_defined_selection'],PARAMETERS['selection_output_filename'])
	nAlign = u_align.n_atoms

	# ----------------------------------------
	# INITIALIZING PARAMETER VARIABLES
	if PARAMETERS['trajectory_file_list'] == None:
		start = int(parameters['start'])
		end = int(parameters['end'])

		trajectory_file_list = []
		temp = start
		while temp <= end:
			trajectory_file_list.append(PARAMETERS['trajectory_location_string']%(temp,temp))	# HUGE SOURCE OF ANNOYANCE... NOT EVERYONE ORGANIZES TRAJECTORIES IN THIS SPECIFIC FORMAT/LOCATION/ETC... potential solution: read in a file that contains the global location for all trajectory files to be analyzed
	else: print 'User has read in a trajectory file list.'

	# ----------------------------------------
	# DETERMINING nSTEPS TO BE ANALYZED 
	if PARAMETERS['nSteps'] == None:
		nSteps = 0
		for i in trajectory_file_list:
			u.load_new(i)
			nSteps += len(range(0,u.trajectory.n_frames,PARAMETERS['trajectory_stepping_value']))
	else:
		nSteps = int(PARAMETERS['nSteps'])

	# ----------------------------------------
	# CALCULATING ITERATIVE AVERAGE STRUCTURE
	if PARAMETERS['user_defined_average_structure'] == None:
		node_average_coordinates, align_average_coordinates = iterative_avg(u,trajectory_file_list,u_all,u_align,u_analysis_selection_list,PARAMATERS['node_type'],nNodes,nAlign,nSteps,trajectory_stepping_value=PARAMETERS['trajectory_stepping_value'])
	else:
		avg = MDAnalysis.Universe(PARAMETERS['user_defined_average_structure'])		# assumes a pdb file of the alignment and analysis selections, at the very least.
		node_average_coordinates, align_average_coordinates = get_average_coordinates(avg,PARAMETERS['node_type'],parameters['alignment'],PARAMETERS['user_defined_selection'],nNodes)
	
	# ----------------------------------------
	# OUTPUTTING AVG COORDINATES OF ALIGNMENT LANDMARK AND ANALYSIS SELECTION
	###






# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

euclid_dist = importlib.import_module(parameters['distance_functions_file'].split('.')[0],package=None).euclid_dist
#make_selection = importlib.import_module(parameters['selection_functions_file'].split('.')[0],package=None).make_selection

# ----------------------------------------
# MAIN
if __name__ == '__main__':
	main()

