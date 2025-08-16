"""
get_conformationally_equivalent_molecules.py, Geoffrey Weal, 10/11/23

This script is designed to obtain all the conformationally equivalent molecules details from file.
"""

def get_conformationally_equivalent_molecules(molecules_path, original_molecule_names):
	"""
	This method is designed to obtain all the conformationally equivalent molecules details from file.

	Parameters
	----------
	molecules_path : str.
		This is the path to all the molecules in the crystal.
	original_molecule_names : list of ints
		These are the names of the molecules in the crystal. 

	Returns
	------- 
	conformationally_equivalent_molecules : dict. 
		This dictionary contains informatino about which conformationally equivalent molecules are assigned to which conformatinoally unique molecules (and theirfore have the same reorganisation energy data assigned to them)
	"""

	# First, remove solvent data from molecule names in the original_molecule_names list.
	molecule_names = [int(str(molecule_name).replace('S','')) for molecule_name in original_molecule_names]

	# Second, collect the information about which conformationally equivalent molecules are the same as the conformationally unique molecules.
	conformationally_equivalent_molecules = {}
	with open(molecules_path+'/Conformationally_Unique_Molecule_Information.txt','r') as Conformationally_Unique_Molecule_InformationTXT:
		Conformationally_Unique_Molecule_InformationTXT.readline()
		for line in Conformationally_Unique_Molecule_InformationTXT:
			line = line.rstrip().split()
			symmetric_molecule_no = int(line[0].replace('S',''))
			unique_molecule_no    = int(line[2].replace('S',''))
			if (symmetric_molecule_no in molecule_names) and (unique_molecule_no in molecule_names):
				conformationally_equivalent_molecules[symmetric_molecule_no] = unique_molecule_no

	# Third, quick check to make sure Conformationally_Unique_Molecule_Information.txt is all good.
	# Values are unique molecules, Keys are equivalent molecules
	if len(set(conformationally_equivalent_molecules.keys()) & set(conformationally_equivalent_molecules.values())) > 0:
		raise Exception('Error: there are equivalent molecules that have also been assigned as unique in Conformationally_Unique_Molecule_Information.txt.\nCheck this out\nconformationally_equivalent_molecules = '+str(conformationally_equivalent_molecules))

	# Fourth, obtain a list of the unique molecules in the crystal.
	conformationally_unique_molecule_names = [molecule_name for molecule_name in molecule_names if (molecule_name not in conformationally_equivalent_molecules.keys())]

	# Fourth, include unique molecules that are not included as keys in conformationally_equivalent_molecules yet:
	for molecule_name in molecule_names:
		if molecule_name not in conformationally_equivalent_molecules:
			conformationally_equivalent_molecules[molecule_name] = molecule_name

	# Fifth, return conformationally_equivalent_molecules and conformationally_unique_molecule_names
	return conformationally_equivalent_molecules, conformationally_unique_molecule_names