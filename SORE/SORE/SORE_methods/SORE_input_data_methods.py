"""
SORE_input_data_methods.py, Geoffrey Weal, 26/3/22

SORE_input_data_methods contains information for obtaining repetitively obtained initial setup data.
"""

import os
import numpy as np
from ase.io import read, write

SORE_input_data_crystal_filename = 'SORE_input_data.txt'
def save_SORE_input_data(SORE_input_data_crystal_folderpath, molecule_names, molecules, crystal_cell_lattice, all_coupling_values):
	"""
	This method is designed to save the input data obtaind from the initial setup to file.

	Parameters
	----------
	SORE_input_data_crystal_folderpath : str.
		This is the folder path to save the SORE input data to.
	input_EKMC_settings : tuple
		This contains all the EKMC settings.

	"""

	# First, check that the SORE_input_data_filepath does not already exist, we dont want to accidently override it. 
	if os.path.exists(SORE_input_data_crystal_folderpath):
		raise Exception('Error: '+str(SORE_input_data_crystal_folderpath)+' already exists. Will fail program to prevent overriding this folder')

	# Second, make the SORE_input_data folder to store initial setup data into. 
	os.makedirs(SORE_input_data_crystal_folderpath)

	# Fifth, create the input_data file and add data to it. 
	with open(SORE_input_data_crystal_folderpath+'/'+SORE_input_data_crystal_filename, 'w') as SORE_input_dataTXT:
		SORE_input_dataTXT.write('molecule_names: '       + str(molecule_names)+'\n')
		SORE_input_dataTXT.write('crystal_cell_lattice: ' + str(np.array(crystal_cell_lattice).tolist())+'\n')
		SORE_input_dataTXT.write('all_coupling_values: '  + str(all_coupling_values)+'\n')

	# Sixth, write the xyz files for the molecules required for the SORE calculation. 
	for molecule_name, molecule in zip(molecule_names, molecules):
		write(SORE_input_data_crystal_folderpath+'/'+molecule_name+'.xyz', molecule)

def read_SORE_input_data(SORE_input_data_crystal_folderpath):
	"""
	This method is designed to read the input data from file, regarding the information from the initial setup.

	Parameters
	----------
	SORE_input_data_crystal_folderpath : str.
		This is the folder path to save the SORE input data to.

	Returns
	-------

	"""

	# First, check that the folders and files that we need to read exist. 
	if not os.path.exists(SORE_input_data_crystal_folderpath):
		raise Exception('Error: '+str(SORE_input_data_crystal_folderpath)+' does not exist. Check this')
	if not os.path.exists(SORE_input_data_crystal_folderpath+'/'+SORE_input_data_crystal_filename):
		raise Exception('Error: '+str(SORE_input_data_crystal_folderpath+'/'+SORE_input_data_crystal_filename)+' does not exist. Check this')

	# Second, read the input files
	input_data_variables = ('molecule_names', 'crystal_cell_lattice', 'all_coupling_values')
	input_data = []
	with open(SORE_input_data_crystal_folderpath+'/'+SORE_input_data_crystal_filename, 'r') as SORE_input_dataTXT:
		index = 0
		for variable_name in input_data_variables:
			line = SORE_input_dataTXT.readline().rstrip()
			if variable_name not in line:
				raise Exception('Error: Can not find '+str(variable_name)+' in line '+str(index+1)+' of '+str(SORE_input_data_crystal_folderpath+'/'+SORE_input_data_crystal_filename))
			data = eval(line.replace(variable_name+': ',''))
			input_data.append(data)
			index += 1

	# Third, 
	molecule_names, crystal_cell_lattice, all_coupling_values = input_data

	# Fourth, obtain the molecules files
	molecules  = []
	for molecule_name in molecule_names:
		if not os.path.exists(SORE_input_data_crystal_folderpath+'/'+str(molecule_name)+'.xyz'):
			raise Exception('Error: '+SORE_input_data_crystal_folderpath+'/'+str(molecule_name)+'.xyz'+' does not exist. Check this')
		molecule = read(SORE_input_data_crystal_folderpath+'/'+str(molecule_name)+'.xyz')
		molecules.append(molecule)

	# Fourth, return 
	return molecule_names, molecules, crystal_cell_lattice, all_coupling_values









