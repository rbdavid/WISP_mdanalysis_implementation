
import MDAnalysis
import sys

# ----------------------------------------
# CLASSIFICATION OF NONSTANDARD RESIDUES:
# ----------------------------------------
nucleic = ['A5','A3','A','G5','G3','G','C5','C3','C','T5','T3','T','U5','U3','U']
triphosphate = ['atp','adp','PHX']
other = ['MG']

# ----------------------------------------
# HOMEMADE ATOM SELECTION STRINGS FOR THE NONSTANDARD RESIDUES
# ----------------------------------------
sugar = "name C5' C4' O4' C1' C3' C2' O2' " + " C5* C4* O4* C1* C3* O3* C2* O2* "		# NO HYDROGENS; DOES NOT INCLUDE THE O5' atom (which I will include in the phosphate atom selection string...; the atoms with * are found in triphosphates;
sugar_5= sugar + " O5'"		# NO HYDROGENS
sugar_3= sugar + " O3' "	# NO HYDROGENS
base = 'name N9 C8 N7 C5 C6 N6 N1 C2 N3 C4 O6 N4 C2 O2 O4'	# NO HYDROGENS; selection string that will select all appropriate atoms for any of the nucleic residues...
a_phos = 'name O5* O2A O1A PA O3A'
b_phos = 'name PB O1B O2B O3B'
g_phos = 'name PG O1G O2G O3G'
inorg_phos = 'name P O1 O2 O3 O4'	# NO HYDROGENS

# ----------------------------------------
# FUNCTION USED TO MAKE ANY OF THE HOMEMADE ATOM SELECTIONS FOR THE NONSTANDARD RESIDUES
# ----------------------------------------

def make_selections(analysis_universe,file_name,node_definition,selection_string,source_selection_strings,sink_selection_strings,nonstandard_substrates_selection=None,homemade_selections=None):
        """
        """
        selection_list = []
        count = 0
        source_index_list = []
        sink_index_list = []
        with open(file_name,'w') as f:

                # ----------------------------------------
                # SUBSTRATE SELECTION - CENTER OF MASS OF RESIDUES
                # ----------------------------------------
                if node_definition.upper() == 'COM':
                        substrate_selection = analysis_universe.select_atoms(selection_string)
                        nResidues_range = range(substrate_selection.n_residues)
                        for i in nResidues_range:
                                temp = substrate_selection.residues[i].atoms
                                selection_list.append(temp)
                                
                                temp_resname = temp.resnames[0]
                                temp_resid = temp.resids[0]
                                node_string = '%s_%s' %(temp_resname,temp_resid)
                                if node_string in source_selection_strings:
                                        f.write("%02d   %s   %s # SOURCE\n" %(count,temp_resname,temp_resid))
                                        source_index_list.append(count)
                                elif node_string in sink_selection_strings:
                                        f.write("%02d   %s   %s # SINK\n" %(count,temp_resname,temp_resid))
                                        sink_index_list.append(count)
                                else:
                                        f.write("%02d   %s   %s\n" %(count,temp_resname,temp_resid))
                                
                                count += 1

                # ----------------------------------------
                # SUBSTRATE SELECTION - ATOMS   ### NEED TO DEBUG/CHECK THIS CODE
                # ----------------------------------------
                elif node_definition.upper() == 'ATOMIC':
                        u_substrate_selection = analysis_universe.select_atoms(selection_string)
                        nAtoms_range = range(substrate_selection.n_atoms)
                        for i in nAtoms_range:
                                temp = u_substrate_selection.residues[i].atoms
                                selection_list.append(temp)

                                temp_resname = temp.resnames[0]
                                temp_resid = temp.resids[0]
                                node_string = '%s_%s' %(temp_resname,temp_resid)
                                if node_string in source_selection_strings:
                                        f.write("%02d   %s   %s # SOURCE\n" %(count,temp_resname,temp_resid))
                                        source_index_list.append(count)
                                elif node_string in sink_selection_strings:
                                        f.write("%02d   %s   %s # SINK\n" %(count,temp_resname,temp_resid))
                                        sink_index_list.append(count)
                                else:
                                        f.write("%02d   %s   %s\n" %(count,temp_resname,temp_resid))
                                
                                count += 1
                
                # ----------------------------------------
                # SUBSTRATE SELECTION - MORE COMPLEX STUFF... ###
                # ----------------------------------------
                else:
                        print 'The substrate_node_definition parameter is not understood.'
                        sys.exit()

                # ----------------------------------------
                # NONSTANDARD SUBSTRATES - USER DEVELOPED SELECTIONS
                # ----------------------------------------
                count = len(selection_list)
                if nonstandard_substrates_selection != None:
                        nonstandard_substrates = analysis_universe.select_atoms(nonstandard_substrates_selection)
                        nResidues_range = range(nonstandard_substrates.n_residues)
                	for i in nResidues_range:
                		temp_resname = nonstandard_substrates.residues[i].resname
                		temp_resid = nonstandard_substrates.residues[i].resid
                		if temp_resname in homemade_selections:
                			make_nonstandard_selections(analysis_universe,temp_resname,temp_resid,f,source_selection_strings,sink_selection_strings,selection_list,source_index_list,sink_index_list,count)
                                        count = len(selection_list)
                		else:
                			temp = analysis_universe.select_atoms('resname %s and resid %d' %(temp_resname,temp_resid))
                			selection_list.append(temp)
                                
                                        node_string = '%s_%s' %(temp_resname,temp_resid)
                                        if node_string in source_selection_strings:
                                                f.write("%02d   %s   resid %s # SOURCE\n" %(count,temp_resname,temp_resid))
                                                source_index_list.append(count)
                                        elif node_string in sink_selection_strings:
                                                f.write("%02d   %s   resid %s # SINK\n" %(count,temp_resname,temp_resid))
                                                sink_index_list.append(count)
                                        else:
                                                f.write("%02d   %s   resid %s\n" %(count,temp_resname,temp_resid))
                                        
                                        count = len(selection_list)

        return selection_list, source_index_list, sink_index_list

def make_nonstandard_selections(analysis_universe,resname,resid,output_file,source_selection_strings,sink_selection_strings,selection_list,source_index_list,sink_index_list,count):
	"""A function that takes in a residue name and creates a non-standard MDAnalysis atom selection
	
	Usage: make_nonstandard_selection(........)
	
	Arguments:
		analysis_universe: MDAnalysis Universe object to be used as the analysis universe.
		resname: string of the residue name;
		resid: int of the residue ID number;
		output_file: file object that is to be written to;
                ...
        """

	# ----------------------------------------
	# CREATING THE NUCLEIC SELECTIONS
	# ----------------------------------------
	if resname in nucleic:
		
                # CREATING THE SLECTION FOR THE BASE OF NUCLEIC RESIDUES
		sel_string = 'resname %s and resid %d and %s' %(resname,resid,base)
		temp = analysis_universe.select_atoms(sel_string)
		selection_list.append(temp)
                
                node_string = '%s_%s_Base' %(resname,resid)
                if node_string in source_selection_strings:
                        output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                        source_index_list.append(count)
                elif node_string in sink_selection_strings:
                        output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                        sink_index_list.append(count)
                else:
                        output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
                                        
		count +=1

		# CREATING THE SLECTION FOR THE SUGAR OF NUCLEIC RESIDUES
		if resname in ['A5','G5','C5','T5','C5']:
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar_5)
			temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(temp)

                        node_string = '%s_%s_Sugar' %(resname,resid)
                        if node_string in source_selection_strings:
                                output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                                source_index_list.append(count)
                        elif node_string in sink_selection_strings:
                                output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                                sink_index_list.append(count)
                        else:
                                output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
			
                        count +=1
			return

		elif resname in ['A3','U3','C3','G3']:
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar_3)
			temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(temp)

                        node_string = '%s_%s_Sugar' %(resname,resid)
                        if node_string in source_selection_strings:
                                output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                                source_index_list.append(count)
                        elif node_string in sink_selection_strings:
                                output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                                sink_index_list.append(count)
                        else:
                                output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
			
			count +=1

		else:
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar)
			temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(temp)

                        node_string = '%s_%s_Sugar' %(resname,resid)
                        if node_string in source_selection_strings:
                                output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                                source_index_list.append(count)
                        elif node_string in sink_selection_strings:
                                output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                                sink_index_list.append(count)
                        else:
                                output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
			
			count +=1

		# CREATING THE SLECTION FOR THE PHOSPHATE OF NUCLEIC RESIDUES
		sel_string = "(resname %s and resid %s and name P OP1 OP2 O5') or (resid %s and name O3')" %(resname,resid,resid-1)
		temp = analysis_universe.select_atoms(sel_string)
                print temp.n_atoms
		selection_list.append(temp)

                node_string = '%s_%s_Phosphate' %(resname,resid)
                if node_string in source_selection_strings:
                        output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                        source_index_list.append(count)
                elif node_string in sink_selection_strings:
                        output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                        sink_index_list.append(count)
                else:
                        output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
			
		count += 1
		return

	# ----------------------------------------
	# CREATING THE TRIPHOSPHATE ATOM SELECTIONS
	elif resname in triphosphate:
		if resname in ['atp','adp']:
                
                        # CREATING THE SLECTION FOR THE BASE OF TRIPHOSPHATES
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,base)
			temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(temp)
                
                        node_string = '%s_%s_Base' %(resname,resid)
                        if node_string in source_selection_strings:
                                output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                                source_index_list.append(count)
                        elif node_string in sink_selection_strings:
                                output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                                sink_index_list.append(count)
                        else:
                                output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
                                        
			count +=1

		        # CREATING THE SLECTION FOR THE SUGAR OF TRIPHOSPHATES
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,sugar)
			temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(temp)

                        node_string = '%s_%s_Sugar' %(resname,resid)
                        if node_string in source_selection_strings:
                                output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                                source_index_list.append(count)
                        elif node_string in sink_selection_strings:
                                output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                                sink_index_list.append(count)
                        else:
                                output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
			
			count +=1

		# CREATING THE SLECTION FOR THE PHOSPHATES OF TRIPHOSPHATES
		if resname == 'atp':
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,a_phos)
			temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(temp)

                        node_string = '%s_%s_A_Phosphate' %(resname,resid)
                        if node_string in source_selection_strings:
                                output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                                source_index_list.append(count)
                        elif node_string in sink_selection_strings:
                                output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                                sink_index_list.append(count)
                        else:
                                output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
			
			count +=1

			sel_string = 'resname %s and resid %d and %s' %(resname,resid,b_phos)
			temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(temp)

                        node_string = '%s_%s_B_Phosphate' %(resname,resid)
                        if node_string in source_selection_strings:
                                output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                                source_index_list.append(count)
                        elif node_string in sink_selection_strings:
                                output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                                sink_index_list.append(count)
                        else:
                                output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
			
			count +=1

			sel_string = 'resname %s and resid %d and %s' %(resname,resid,g_phos)
			temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(temp)

                        node_string = '%s_%s_G_Phosphate' %(resname,resid)
                        if node_string in source_selection_strings:
                                output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                                source_index_list.append(count)
                        elif node_string in sink_selection_strings:
                                output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                                sink_index_list.append(count)
                        else:
                                output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
			
			count +=1
			return

		elif resname == 'adp':
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,a_phos)
			temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(temp)

                        node_string = '%s_%s_A_Phosphate' %(resname,resid)
                        if node_string in source_selection_strings:
                                output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                                source_index_list.append(count)
                        elif node_string in sink_selection_strings:
                                output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                                sink_index_list.append(count)
                        else:
                                output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
			
			count +=1

			sel_string = 'resname %s and resid %d and %s' %(resname,resid,b_phos)
			temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(temp)

                        node_string = '%s_%s_B_Phosphate' %(resname,resid)
                        if node_string in source_selection_strings:
                                output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                                source_index_list.append(count)
                        elif node_string in sink_selection_strings:
                                output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                                sink_index_list.append(count)
                        else:
                                output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
			
			count +=1
			return

		# CREATING THE SLECTION FOR INORGANIC PHOSPHATE MOLECULE
		elif resname == 'PHX':
			sel_string = 'resname %s and resid %d and %s' %(resname,resid,inorg_phos)
			temp = analysis_universe.select_atoms(sel_string)
			selection_list.append(temp)

                        node_string = '%s_%s' %(resname,resid)
                        if node_string in source_selection_strings:
                                output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                                source_index_list.append(count)
                        elif node_string in sink_selection_strings:
                                output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                                sink_index_list.append(count)
                        else:
                                output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
			
			count +=1
			return

	# ----------------------------------------
	# CREATING ANY REMAINING ATOM SELECTIONS...
	elif resname in other:
		sel_string = 'resname %s and resid %d' %(resname,resid)
		temp = analysis_universe.select_atoms(sel_string)
		selection_list.append(temp)

                node_string = '%s_%s' %(resname,resid)
                if node_string in source_selection_strings:
                        output_file.write("%02d   %s   %s # SOURCE\n" %(count,resname,sel_string))
                        source_index_list.append(count)
                elif node_string in sink_selection_strings:
                        output_file.write("%02d   %s   %s # SINK\n" %(count,resname,sel_string))
                        sink_index_list.append(count)
                else:
                        output_file.write("%02d   %s   %s\n" %(count,resname,sel_string))
			
		count +=1
		return

# ----------------------------------------












## ----------------------------------------
## OLD FUNCTIONS
## ----------------------------------------
#
#def old_make_selection(u,node_type,user_defined_selection,selection_output_filename):
#	selection_list = []
#
#	if node_type == 'atomic':
#		selection = u.select_atoms(user_defined_selection)	# atomic selection string; most dangerous selection string... 
#		nAtoms = selection.n_atoms 
#		temp = '%0'+'%d'%(int(np.log10(nAtoms))+1)+'d'	# for correct zero pad when numbering selections
#		with open(selection_output_filename,'w') as W:
#			W.write('# selection_num   Atom_num   Atom_name   Resid\n')
#			for i in range(nAtoms):
#				selection_list.append(selection.atoms[i])
#				W.write(temp%(i)+'   %d   %s   %d\n' %(selection.atoms[i].number,selection.atoms[i].name,selection.atoms[i].resid))
#		return selection, selection_list, nAtoms
#
#	elif node_type == 'residue_COM':
#		selection = u.select_atoms(user_defined_selection)	# does not necessitate only protein selection; assumes only one node per residue
#		nResidues = selection.n_residues
#		order_of_magnitude = int(np.log10(nResidues))+1		# for correct zero pad when numbering selections
#		temp = '%0'+'%d'%(order_of_magnitude)+'d'
#		with open(selection_output_filename,'w') as W:
#			W.write('# selection_num   Res_name   Resid   Num_atoms ! CENTER OF MASS SELECTION\n')
#			for i in range(nResidues):
#				selection_list.append(selection.residues[i])	# want to select all atoms in residue i of selection;
#				W.write(temp%(i)+'   %s   %d   %d\n'%(selection.residues[i].resname,selection.residues[i].resid,selection.residues[i].n_atoms))
#		return selection, selection_list, nResidues
#
#	elif node_type == 'backbone_COM':
#		selection = u.select_atoms('(backbone or nucleicbackbone) and '+user_defined_selection)	# does not necessitate only protein selection; assumes only one node per residue.
#		nResidues = selection.n_residues
#		order_of_magnitude = int(np.log10(nResidues))+1		# for correct zero pad when numbering selections
#		temp = '%0'+'%d'%(order_of_magnitude)+'d'
#		with open(selection_output_filename,'w') as W:
#			W.write('# selection_num   Res_name   Resid   Num_atoms ! CENTER OF MASS OF BACKBONE SELECTION\n')
#			for i in range(nResidues):
#				temp_atom_group = selection[selection.resids == selection.residues[i].resid]	# only want to include atoms in that SHOULD BE in the atom selection. 
#				selection_list.append(temp_atom_group)
#				W.write(temp%(i)+'   %s   %d   %d\n'%(temp_atom_group.atoms[0].resname, temp_atom_group.atoms[0].resid, temp_atom_group.n_atoms))
#		return selection, selection_list, nResidues
#
#	elif node_type == 'sidechain_COM':
#		selection = u.select_atoms('not (backbone or nucleicbackbone) and '+user_defined_selection)	# does not necessitate only protein selection; assumes only one node per residue.
#		nResidues = selection.n_residues
#		order_of_magnitude = int(np.log10(nResidues))+1		# for correct zero pad when numbering selections
#		temp = '%0'+'%d'%(order_of_magnitude)+'d'
#		with open(selection_output_filename,'w') as W:
#			W.write('# selection_num   Res_name   Resid   Num_atoms ! CENTER OF MASS OF SIDECHAIN SELECTION\n')
#			for i in range(nResidues):
#				temp_atom_group = selection[selection.resids == selection.residues[i].resid]	# only want to include atoms in that SHOULD BE in the atom selection. 
#				selection_list.append(temp_atom_group)
#				W.write(temp%(i)+'   %s   %d   %d\n'%(temp_atom_group.atoms[0].resname, temp_atom_group.atoms[0].resid, temp_atom_group.n_atoms))
#		return selection, selection_list, nResidues
#		
#	### STILL TO DO:
#	# MORE COMPLEX, USER DEFINED ATOM SELECTIONS; good example: multiple sites per residue (i.e. backbone and sidechain COM for each residue of interest, with the correct node ordering based on resid)
#
#	else:
#		print "Parameter value for 'node_type' does not equal an expected string. Please check the parameter and try again."
#		sys.exit()
#
