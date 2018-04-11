#!
#!/home/rbdavid/bin/python

# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import MDAnalysis
from MDAnalysis.analysis.align import *
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
eigen = np.linalg.eig

necessary_parameters = ['ref_pdb','pdb','traj_loc','start','end','Wrapped','outputname','selection_file','resid_offset']
all_parameters = ['ref_pdb','pdb','traj_loc','start','end','Wrapped','outputname','selection_file','resid_offset','alignment','wrapping_selection','substrates','homemade_selections','write_summary','summary_filename','selection_output']

# ----------------------------------------
# FUNCTIONS: 
# ----------------------------------------

def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
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

def summary():
	with open(parameters['summary_filename'],'w') as f:
		f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
		f.write('To recreate this analysis, run this line:\n')
		for i in range(len(sys.argv)):
			f.write('%s ' %(sys.argv[i]))
		f.write('\n\n')
		f.write('Parameters used:\n')
		for i in all_parameters:
			f.write('%s = %s \n' %(i,parameters[i]))
		f.write('\n\n')
		f.write('output is written to:\n')
		f.write('	%s\n' %(parameters['outputname']))
		f.write('\nAtom selections analyzed have been written out to %s\n' %(parameters['selection_output']))


def make_selection(u,node_type,user_defined_selection,selection_output_filename):
	"""
	I need to write stuff here about what is going on...
	"""
	selection_list = []

	if node_type == 'atomic':
		selection = u.select_atoms(user_defined_selection)	# atomic selection string; MOST DANGEROUS SELECTION STRING... e.g. user_defined_selection = 'resid 20:651 and name CA CG P PA PB PG'
		nAtoms = selection.n_atoms 
		# if nAtoms > 2500:
		# 	print holy cow, too many nodes
		#	sys.exit because holy cow that is too many nodes
		order_of_magnitude = int(np.log10(nAtoms))+1		# for correct zero pad when numbering selections
		temp = '%0'+'%d'%(order_of_magnitude)+'d'		# ex: nAtoms = 150, order_of_magnitude = 3, temp = '%03d'
		with open(selection_output_filename,'w') as W:
			W.write('# selection_num   Atom_num   Atom_name   Resid\n')
			for i in range(nAtoms):
				selection_list.append(selection.atoms[i])
				W.write(temp%(i)+'   %d   %s   %d\n' %(selection.atoms[i].number,selection.atoms[i].name,selection.atoms[i].resid))
		return selection, selection_list, nAtoms

	elif node_type == 'residue_COM':
		selection = u.select_atoms(user_defined_selection)	# does not necessitate only protein selection...
		nResidues = selection.n_residues
		order_of_magnitude = int(np.log10(nResidues))+1		# for correct zero pad when numbering selections
		temp = '%0'+'%d'%(order_of_magnitude)+'d'
		with open(selection_output_filename,'w') as W:
			W.write('# selection_num   Res_name   Resid   Num_atoms ! CENTER OF MASS SELECTION\n')
			for i in range(nResidues):
				selection_list.append(selection.residues[i])
				W.write(temp%(i)+'   %s   %d   %d\n'%(selection.residues[i].resname,selection.residues[i].resid,selection.residues[i].n_atoms))
		return selection, selection_list, nResidues

	elif node_type == 'backbone_COM':
		selection = u.select_atoms('(backbone or nucleicbackbone) and '+user_defined_selection)	# does not necessitate only protein selection...
		nResidues = selection.n_residues
		order_of_magnitude = int(np.log10(nResidues))+1		# for correct zero pad when numbering selections
		temp = '%0'+'%d'%(order_of_magnitude)+'d'
		with open(selection_output_filename,'w') as W:
			W.write('# selection_num   Res_name   Resid Num_atoms! CENTER OF MASS OF BACKBONE SELECTION\n')
			for i in range(nResidues):
				temp_atom_group = selection[selection.resids == selection.residues[i].resid]	# only one backbone selection per residue, thus no bug here...
				selection_list.append(temp_atom_group)
				W.write(temp%(i)+'   %s   %d   %d\n'%(temp_atom_group.atoms[0].resname, temp_atom_group.atoms[0].resid, temp_atom_group.n_atoms))
		return selection, selection_list, nResidues

	elif node_type == 'sidechain_COM':
		selection = u.select_atoms('not (backbone or nucleicbackbone) and '+user_defined_selection)	# does not necessitate only protein selection...
		nResidues = selection.n_residues
		order_of_magnitude = int(np.log10(nResidues))+1		# for correct zero pad when numbering selections
		temp = '%0'+'%d'%(order_of_magnitude)+'d'
		with open(selection_output_filename,'w') as W:
			W.write('# selection_num   Res_name   Resid Num_atoms! CENTER OF MASS OF SIDECHAIN SELECTION\n')
			for i in range(nResidues):
				temp_atom_group = selection[selection.resids == selection.residues[i].resid]	# only one sidechain selection per residue, thus no bug here...
				selection_list.append(temp_atom_group)
				W.write(temp%(i)+'   %s   %d   %d\n'%(temp_atom_group.atoms[0].resname, temp_atom_group.atoms[0].resid, temp_atom_group.n_atoms))
		return selection, selection_list, nResidues
	
	elif node_type == 'BACKBONE_AND_SIDECHAIN_COMS':
		selection1 = u.select_atoms('(backbone or nucleicbackbone) and '+user_defined_selection)	# does not necessitate only protein selection...
		selection2 = u.select_atoms('not (backbone or nucleicbackbone) and '+user_defined_selection)	# does not necessitate only protein selection...
		all_selection = u.select_atoms(user_defined_selection)

	### STILL TO DO:
	# MORE COMPLEX, USER DEFINED ATOM SELECTIONS; good example: multiple sites per residue (i.e. backbone and sidechain COM for each residue of interest, with the correct node ordering based on resid)

	else:
		print "Parameter value for 'node_type' does not equal an expected string. Please check the parameter and try again."


def iterative_avg(u,u_align,u_analysis_selection,nNodes,thresh=1E-5,maxIter=100):
	



def main():

####code breakdown
# 1) load in analysis universe, make atom selections to be used later on
# 2) calc iterative, average structure using an alignment landmark if required; only necessary for the atom selection to be analyzed; add capability to read in user defined average structure.
# 3) calc the distance between node(t) and <node> for every t; produces d matrix (nSteps*nNodes)
#	--> also calculate or read in node-node contact map (square matrix)...
# 4) calc the covariance/correlation matrices between d(i;t) and d(j;t)
# 5) network analysis/graph stuff...

	# ----------------------------------------
	# INITIALIZING MDANALYSIS UNIVERSES 
	u = MDAnalysis.Universe(parameters['pdb'])
	u_align = u.select_atoms(parameters['alignment'])
	u_analysis_selection, u_analysis_selection_list, nNodes = make_selection(u,parameters['node_type'],parameters['user_defined_selection'],parameters['selection_output_filename'])

	ref = MDAnalysis.Universe(parameters['ref_pdb'])
	ref_align = ref.select_atoms(parameters['alignment'])
	ref_align_positions = ref_align.positions


#	if parameters['nSteps'] == None:
#		temp = int(parameters['start'])
#		nSteps = 0
#		while temp <= int(parameters['end']):
#			u.load_new(


# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

euclid_dist = importlib.import_module(parameters['distance_functions_file'].split('.')[0],package=None).euclid_dist
make_selection = importlib.import_module(parameters['selection_functions_file'].split('.')[0],package=None).make_selection

# ----------------------------------------
# MAIN
if __name__ == '__main__':
	main()

