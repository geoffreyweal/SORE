"""
get_EET_coupling_data.py, Geoffrey Weal, 19/4/22

This script will obtain the electronic coupling energies as calculated using the eet function in Gaussian between pairs of neighbouring molecules (dimers).
"""
import os

from ase.io import read

from EKMC.EKMC_Setup.EKMC_Only_Setup.get_RE_and_bandgap_data_methods.helpful_functions import get_reorganisation_energy, convert_hartree_to_eV

def get_RE_and_bandgap_data(reorganisation_and_bandgap_energy_details, crystal_name, functional_and_basis_set, conformationally_equivalent_molecules, original_molecule_names, conformationally_unique_molecule_names):
	"""
	This method will obtain the ground and excited sstructure energy gaps requierd for obtaining reorganisation energies between molecules in crystal.

	Parameters
	----------
	reorganisation_and_bandgap_energy_details : str.
		These are the details about the reorganisation energies between dimers and the bandgap energies of molecules data.
	crystal_name : str.
		This is the name of the crystal.
	functional_and_basis_set : str.
		This is the name of the functional and basis set used in calculations to obtain EET and ATC information. 
	molecules_path : str.
		This is the path to all the molecules in the crystal.
	original_molecule_names : list of ints
		These are the names of the molecules in the crystal. 

	Returns
	------- 
	reorganisation_energy_data : dict.
		These are all the reorgansation energies for the dimers in the crystal.
	bandgap_energy_data : dict.
		These are all the bandgap energies for the molecules in the crystal.
	conformationally_equivalent_molecules : dict. 
		This dictionary contains informatino about which conformationally equivalent molecules are assigned to which conformatinoally unique molecules (and theirfore have the same reorganisation energy data assigned to them)
	"""

	# First, remove solvent data from molecule names in the original_molecule_names list.
	molecule_names = [int(str(molecule_name).replace('S','')) for molecule_name in original_molecule_names]

	# Second , obtain conformationally_equivalent_molecules_reverse, where unique molecules are the keys and equivalent molecule are the values. 
	conformationally_equivalent_molecules_reverse = {}
	for key, value in conformationally_equivalent_molecules.items():
		conformationally_equivalent_molecules_reverse.setdefault(value,[]).append(key)

	# Third, set up the reorganisation_energy_data and the bandgap_energy_data dictionaries. 
	reorganisation_energy_data = {(mol_1, mol2): None for mol_1, mol2 in zip(molecule_names, molecule_names)}
	bandgap_energy_data        = {mol_name:      None for mol_name    in molecule_names}

	# Fourth, if reorganisation_and_bandgap_energy_details['path_to_RE_folder'] exists, gather reorganisation energy and bandgap energy from reorganisation_and_bandgap_energy_details
	if 'path_to_RE_folder' in reorganisation_and_bandgap_energy_details:

		# 3.1: Get the path to the reorganisation energy folder.
		reorganisation_energy_path = reorganisation_and_bandgap_energy_details['path_to_RE_folder']

		# 3.2: Check if the relavant Egaps.txt file exists in the Individual_RE_Data folder of either the Unique_RE_Gaussian_Jobs or All_RE_Gaussian_Jobs folder. 
		if not os.path.exists(reorganisation_energy_path):
			raise Exception('Error: Can not find the reorganisation energy data folder: '+str(reorganisation_energy_path))
		if not os.path.exists(reorganisation_energy_path+'/Individual_RE_Data'):
			raise Exception('Error: Can not find the reorganisation energy data folder: '+reorganisation_energy_path+'/Individual_RE_Data')
		energies_TXT_filename = crystal_name+'_'+functional_and_basis_set+'_energies.txt'
		if not os.path.exists(reorganisation_energy_path+'/Individual_RE_Data/'+energies_TXT_filename):
			raise Exception('Error: Can not find the energies.txt file (contains energy gap data for molecule in crystal): '+reorganisation_energy_path+'/Individual_RE_Data/'+energies_TXT_filename)
		if not os.path.isfile(reorganisation_energy_path+'/Individual_RE_Data/'+energies_TXT_filename):
			raise Exception('Error: energies.txt (contains energy gap data for molecule in crystal) is not a file (maybe a folder?): '+reorganisation_energy_path+'/Individual_RE_Data/'+energies_TXT_filename)

		# 3.3: Record the energy for ground and excited structures for each molecule in the crystal.
		reorganisation_energy_molecule_data = {}
		with open(reorganisation_energy_path+'/Individual_RE_Data/'+energies_TXT_filename,'r') as EgapTXT:
			EgapTXT.readline()
			for line in EgapTXT:
				molecule_name, eGS_gGS_energy, eGS_gES_energy, eES_gGS_energy, eES_gES_energy, band_gap_energy = line.rstrip().split()
				molecule_no = int(molecule_name.replace('molecule_','').replace('S',''))
				if molecule_no in reorganisation_energy_molecule_data:
					raise Exception('Error: Molecule'+str(molecule_no)+' seems to appear twice in '+reorganisation_energy_path+'/Individual_RE_Data/'+energies_TXT_filename)
				reorganisation_energy_molecule_data[molecule_no] = (float(eGS_gGS_energy), float(eGS_gES_energy), float(eES_gGS_energy), float(eES_gES_energy))

		# 3.4: Add corresponding reorganisation energies and band-gap energies to reorganisation_energy_data and bandgap_energy_data.
		for index1, mol1_name in enumerate(molecule_names):

			# 3.4.1: Obtain the conformationally unique molecule name for molecule 1.
			unique_mol1_name = conformationally_equivalent_molecules[mol1_name]

			# 3.4.2: Get the energies of the excited and ground state of molecule 1.
			mol1_eGS_gGS, mol1_eGS_gES, mol1_eES_gGS, mol1_eES_gES = reorganisation_energy_molecule_data[unique_mol1_name]

			# 3.4.3: Go through the second molecule to get the reorganisation energies for all the dimer pairs. 
			for index2, mol2_name in enumerate(molecule_names, start=index1):

				# 3.4.4: Obtain the conformationally unique molecule name for molecule 2.
				unique_mol2_name = conformationally_equivalent_molecules[mol2_name]

				# 3.4.5: Get the energies of the excited and ground state of molecule 2.
				mol2_eGS_gGS, mol2_eGS_gES, mol2_eES_gGS, mol2_eES_gES = reorganisation_energy_molecule_data[unique_mol2_name]

				# 3.4.6: Get the dimer with arrangement (mol1_name, mol2_name), check it is not already in reorganisation_energy_data, and record the reorganisation energy for a molecule moving from mol1 --> mol2.
				dimer_pair1 = (mol1_name, mol2_name)
				reorganisation_energy_data[dimer_pair1] = convert_hartree_to_eV(get_reorganisation_energy(mol1_eGS_gGS, mol1_eGS_gES, mol2_eES_gGS, mol2_eES_gES))

				# 3.4.7: If mol1_name == mol2_name, dont need to do the next part of this for loop.
				if mol1_name == mol2_name:
					continue

				# 3.4.8: Get the dimer with arrangement (mol2_name, mol1_name), check it is not already in reorganisation_energy_data, and record the reorganisation energy for a molecule moving from mol2 --> mol1.
				dimer_pair2 = (mol2_name, mol1_name)
				reorganisation_energy_data[dimer_pair2] = convert_hartree_to_eV(get_reorganisation_energy(mol2_eGS_gGS, mol2_eGS_gES, mol1_eES_gGS, mol1_eES_gES))

		# 3.5: Record the band-gap energies of the molecules in the crystal.
		for mol_name in molecule_names:

			# 3.5.1: Obtain the conformationally unique molecule name for the molecule.
			unique_mol_name = conformationally_equivalent_molecules[mol_name]

			# 3.5.2: Get the energies of the excited and ground state of the molecule.
			mol_eGS_gGS, mol_eGS_gES, mol_eES_gGS, mol_eES_gES = reorganisation_energy_molecule_data[unique_mol_name]

			# 3.5.3: Obtain the band-gap for the molecule and add it to bandgap_energy_data.
			bandgap_energy_data[mol_name] = convert_hartree_to_eV(mol_eES_gES - mol_eGS_gGS)

	# ==================================================================================================================================

	# Fifth, if reorganisation energy data is given manually, record this and do consistancy checks. 
	if 'manual_RE_data' in reorganisation_and_bandgap_energy_details:

		# 5.1: Record the manually entered reorganisation energy data.
		for key, manual_reorganisation_energy in reorganisation_and_bandgap_energy_details['manual_RE_data'].items():

			# 5.2: If an int is given, make this a tuple of length two with the same molecule name in both:
			if isinstance(key, int) or isinstance(key, str):
				key = (key, key)

			# 5.3: Check the key has been given correctly. 
			if (not isinstance(key, tuple)) or (not len(key) == 2):
				toString  = "Error, your keys for reorganisation_and_bandgap_energy_details['manual_RE_data'] needs to be a tuple of length 2.\n"
				toString += "reorganisation_and_bandgap_energy_details['manual_RE_data'] = "+str(reorganisation_and_bandgap_energy_details['manual_RE_data'])+"\n"
				toString += "Check this."
				raise Exception(toString)

			# 5.4: Convert the molecule names into ints without solvent indicators.
			unique_molecule1_name = int(str(key[0]).replace('S',''))
			unique_molecule2_name = int(str(key[1]).replace('S',''))

			# 5.5: The inputs in manual_RE_data should be for unique molecules. 
			if not ((unique_molecule1_name in conformationally_unique_molecule_names) and (unique_molecule2_name in conformationally_unique_molecule_names)):
				toString  = 'Error: You have some pairs of molecules that are not unique. You can only include unique molecules here:\n'
				toString += "Molecule pairs in reorganisation_and_bandgap_energy_details['manual_RE_data']: "+str(reorganisation_and_bandgap_energy_details['manual_RE_data'].keys())+'\n'
				toString += 'ConformationallyUnique Molecules:      '+str(conformationally_unique_molecule_names)+'\n'
				toString += 'Conformationally Equivalent Molecules: '+str({key: value for key, value in conformationally_equivalent_molecules.items() if not (key == value)})+'\n'
				toString += 'Check this'
				raise Ecception(toString)

			# 5.6: Copy reorganisation energies for equivalent and unique molecules into reorganisation_energy_data
			for mol_1_name in conformationally_equivalent_molecules_reverse[unique_molecule1_name]:
				for mol_2_name in conformationally_equivalent_molecules_reverse[unique_molecule2_name]:

					# 5.6.1: Get the name of the dimer pair
					dimer_pair = (mol_1_name, mol_2_name)

					# 5.6.2: Given manually given reorganisation energy to dimer pair. 
					reorganisation_energy_data[dimer_pair] = manual_reorganisation_energy

	# ==================================================================================================================================

	# Sixth, if bandgap energy data is given manually, record this and do consistancy checks. 
	if 'manual_bandgap_data' in reorganisation_and_bandgap_energy_details:

		# 6.1: Record the manually entered reorganisation energy data.
		for key, manual_bandgap_energy in reorganisation_and_bandgap_energy_details['manual_bandgap_data'].items():

			# 6.2: Check the key has been given correctly. 
			if not (isinstance(key, int) or isinstance(key, str)):
				toString  = "Error, your key for reorganisation_and_bandgap_energy_details['manual_bandgap_data'] must be ints or strings.\n"
				toString += "reorganisation_and_bandgap_energy_details['manual_bandgap_data'] = "+str(reorganisation_and_bandgap_energy_details['manual_bandgap_data'])+"\n"
				toString += "Check this."
				raise Exception(toString)

			# 6.3: Convert the molecule names into ints without solvent indicators.
			unique_molecule_name = int(str(key).replace('S',''))

			# 6.4: The inputs in manual_RE_data should be for unique molecules. 
			if not (molecule1_name in conformationally_unique_molecule_names):
				toString  = "Error: You have some molecules give in reorganisation_and_bandgap_energy_details['manual_bandgap_data'] that are not unique. You can only include unique molecules here:\n"
				toString += "Molecule pairs in reorganisation_and_bandgap_energy_details['manual_bandgap_data']: "+str(reorganisation_and_bandgap_energy_details['manual_bandgap_data'].keys())+'\n'
				toString += 'Conformationally Unique Molecules:      '+str(conformationally_unique_molecule_names)+'\n'
				toString += 'Conformationally Equivalent Molecules: '+str({key: value for key, value in conformationally_equivalent_molecules.items() if not (key == value)})+'\n'
				toString += 'Check this'
				raise Ecception(toString)

			# 6.5: Copy reorganisation energies for equivalent and unique molecules into reorganisation_energy_data.
			for mol_name in conformationally_equivalent_molecules_reverse[unique_molecule_name]:
				bandgap_energy_data[mol_name] = manual_bandgap_energy

	# ==================================================================================================================================

	# Eighth, return dictionaries of the energy gaps for ground and excited structures for each molecule in the crystal. 
	return bandgap_energy_data, reorganisation_energy_data

	# ==================================================================================================================================




