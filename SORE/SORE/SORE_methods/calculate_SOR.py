"""
check_max_jobs_in_queue_after_next_submission.py, Geoffrey Weal, 4/1/2023

This method will check to make sure that the next job can be submitted and not go over the allocated number of maximum jobs.
"""
import numpy as np
from tqdm import tqdm
from itertools import product
from tqdm.contrib.concurrent import process_map
from SORE.SORE.SORE_methods.calculate_SOR_methods.calculate_single_SOR_standard_method                      import calculate_single_SOR_standard_method
from SORE.SORE.SORE_methods.calculate_SOR_methods.auxillary_methods                                         import single_SOR_standard_method_input_generator, single_SOR_analytical_energetic_disorder_method_input_generator, k_b
from SORE.SORE.SORE_methods.calculate_SOR_methods.calculate_single_SOR_analytical_energetic_disorder_method import calculate_single_SOR_analytical_energetic_disorder_method

def calculate_SOR(SORE_mode, molecule_names, molecules, crystal_cell_lattice, kinetic_model, kinetics_details, molecule_bandgap_energy_data, dimer_reorganisation_energy_data, conformationally_equivalent_molecules, coupling_value_data, constant_rate_data, no_of_cpus):
    """
    This program is designed to simulate the movement of an exciton through a OPV crystal system using the sum-over-rates method. 

    Parameters
    ----------
    SORE_mode : str.
        This indicates the type of SOR calculation you want to perform. Options are:
            * standard:                      Calculate SOR using Marcus Theory.
            * numeric_energetic_disorder:    Calculate SOR using Marcus Theory, where site energy disorder in is included numerically (a number of repeated simulations are required).
            * analytical_energetic_disorder: Calculate SOR using Marcus Theory, where site energy disorder in is included analytically.
    molecule_names : list of strings
        This list contains all the names of the molecules in the crystal.
    molecules : list of ase.Atoms
        This is a list of the molecules in the crystal as ase objects.
    crystal_cell_lattice
        This is the cell object that contains the unit cell lattice matrix.
    kinetic_model : str.
        This is the type of kinetic model you would like to use.
    kinetics_details : dict.
        This contains all the information required to obtain the rate constants (except for electronic coupling and ATC charge data).
    molecule_bandgap_energy_data : dict.
        These are the bandgap energies of the unique molecules in the crystal.
    dimer_reorganisation_energy_data : dict.
        These are the reorganisation energies of the unique dimers in the crystal.
    conformationally_equivalent_molecules : dict.
        This dictionary contains information about which molecules are conformationally equivalent to each other. 
    coupling_value_data : dict. 
        This information contains all the coupling values between dimers in the crystals. These dimers are between molecule 1 (in the unit cell) and molecule 2 (in the unit cell or another cell).
    constant_rate_data : list
        This contains all the data required for the rate constant that are the same for all dimers. 
    calculate_diffusion_tensor : bool.
        This indicates if you want to calculate the diffusion tensor (True), or calculate the diffusion coefficient from the sum of scalar diffusion coefficients (False).
    no_of_cpus : int
        This is the number of cpus available to run the SOR method multiple times if you have set include_disorder_in_exciton_acceptor=True.
    """

    # Preamble, there seems to be a weird seeding problem when using multiprocessing. 
    # We will reseed every time this method is used to make sure that everything is random. 
    np.random.seed()

    # First, convert molecule_name from a string to an int.
    molecule_names_ints = [int(str(molecule_name).replace('S','')) for molecule_name in molecule_names]

    # Second, check that all the name of molecules in molecule_names are found in coupling_value_data keys
    if not (sorted(molecule_names_ints) == sorted(coupling_value_data.keys())):
        raise Exception('huh?')
    for mol1, mol1_data in coupling_value_data.items():
        if not (sorted(molecule_names_ints) == sorted(mol1_data.keys())):
            raise Exception('huh?')

    # Third, obtain the constant factors in the kinetic theory you want to use.
    if kinetic_model.lower() == 'marcus':
        M_constant, X_constant = constant_rate_data
    else:
        raise Exception('Huh?')

    # Fourth, obtan the centre of masses for each of the molecules in the crystal. 
    centre_of_masses = {molecule_name: np.array(molecule.get_center_of_mass()) for molecule_name, molecule in sorted(zip(molecule_names_ints, molecules), key=lambda x: x)}

    # Fifth, convert crystal_cell_lattice into a numpy array.
    crystal_cell_lattice_matrix = np.array(crystal_cell_lattice)

    # Sixth, calculate the energetic disorder component to include in the SOR if you have included energetic disorder in the kinetics_details dictionary. 
    if 'energetic_disorder' in kinetics_details:

        # 6.1: Obtain the energetic disorder value. 
        if not ('energetic_disorder' in kinetics_details):
            raise Exception('Error: You have not given a energetic_disorder value in your kinetic_details dictionary, but you have set include_disorder_in_exciton_acceptor to True. For this to be so, you need to include an energetic_disorder')
        energetic_disorder = kinetics_details['energetic_disorder']
        
        # 6.2: Obtain the temperature (if the energetic_disorder has not been set to 0.0).
        if not (energetic_disorder == 0.0):
            if not ('temperature' in kinetics_details):
                raise Exception('Error: If you want to include energetic_disorder in your simulation, you also need to give a Temperature. Add this to your kinetic_details dictionary.')
            temperature = kinetics_details['temperature']

    # Seventh, calculate the energetic_disorder_component for this system.
    if 'energetic_disorder' in kinetics_details:
        energetic_disorder_component = ((energetic_disorder ** 2.0) / (k_b * temperature)) if (not energetic_disorder == 0.0) else 0.0
    else:
        energetic_disorder_component = 0.0

    # Seventh, initialise the dictionary to record diffusion coefficients in.
    diffusion_coefficient_molecules = {}
    diffusion_tensor_molecules      = {}
    all_molecule_k_ijs_details      = {}
    all_average_hopping_distances   = {}

    # Eighth, if you want to include disorder in the exciton acceptor, run this calculation multiple times
    no_of_simulations = 100000 if ('numeric' in SORE_mode) else 1

    # Ninth, initialise a dictionary to store the total chance for an exciton to move between all the structurally unique molecules in the crystal. 
    exciton_donor_to_acceptor_probabilities = {(mol1,mol2): None for mol1, mol2 in product(coupling_value_data.keys(), coupling_value_data.keys())}

    # Tenth, get the method that we will used to calculate the SOR for exciton hopping.
    include_disorder_in_exciton_acceptor = False
    if SORE_mode in ['standard', 'numeric_energetic_disorder']:
        calculate_single_SOR_method = calculate_single_SOR_standard_method
        if SORE_mode == 'numeric_energetic_disorder':
            include_disorder_in_exciton_acceptor = True
    elif SORE_mode == 'analytical_energetic_disorder':
        calculate_single_SOR_method = calculate_single_SOR_analytical_energetic_disorder_method
    else:
        raise Exception('Error: SORE_mode must be set to either "standard", "numeric_energetic_disorder", or "analytical_energetic_disorder". SORE_mode = '+str(SORE_mode)+'. Check this.')

    # Ninth, obtain the diffusion coefficients for molecules in the crystal.
    # 9.1: For the first molecule in the dimer:
    pbar = tqdm(coupling_value_data.items(), leave=False)
    for mol1, mol1_data in pbar:
        pbar.set_description('Collecting Sum-Over-Rates Exciton Diffusion Coefficient for molecule '+str(mol1))

        # 9.2: Obtain the centre of mass for molecule 1.
        mol1_com = centre_of_masses[mol1]

        # 9.3: Get the unique molecule identifier for molecule 1.
        unique_molname_1 = conformationally_equivalent_molecules[mol1] if (mol1 in conformationally_equivalent_molecules) else mol1

        # 9.4: Obtain the bandgap energy for molecule 1.
        bandgap_energy_mol1 = molecule_bandgap_energy_data[unique_molname_1]

        # 9.5: Include disorder in the bandgap for molecule 1 if desired. 
        #      Only used for 
        #        * SORE_mode == 'numeric_energetic_disorder', and 
        #        * bandgap_energy_mol1_with_disorder == bandgap_energy_mol1 for SORE_mode = 'standard'
        bandgap_energy_mol1_with_disorder = bandgap_energy_mol1 - energetic_disorder_component

        # 9.6: Obtain the input data (via generator) for the type of calculation you would like to perform. 
        if SORE_mode in ['standard', 'numeric_energetic_disorder']:
            input_generator = single_SOR_standard_method_input_generator(no_of_simulations, mol1, mol1_data, mol1_com, unique_molname_1, bandgap_energy_mol1_with_disorder, centre_of_masses, conformationally_equivalent_molecules, molecule_bandgap_energy_data, energetic_disorder, dimer_reorganisation_energy_data, crystal_cell_lattice_matrix, kinetic_model, M_constant, X_constant, include_disorder_in_exciton_acceptor)
        elif SORE_mode == 'analytical_energetic_disorder':
            temperature = kinetics_details['temperature']
            input_generator = single_SOR_analytical_energetic_disorder_method_input_generator(no_of_simulations, mol1, mol1_data, mol1_com, unique_molname_1, bandgap_energy_mol1, centre_of_masses, conformationally_equivalent_molecules, molecule_bandgap_energy_data, energetic_disorder, dimer_reorganisation_energy_data, crystal_cell_lattice_matrix, kinetic_model, temperature)
        else:
            raise Exception('Error: SORE_mode must be set to either "standard", "numeric_energetic_disorder", or "analytical_energetic_disorder". SORE_mode = '+str(SORE_mode)+'. Check this.')

        # 9.7: Run the SOR for all simulation.
        if no_of_simulations > 1:
            mp_data = process_map(calculate_single_SOR_method, input_generator, total=no_of_simulations, unit='sims', leave=False)
        else:
            mp_data = [calculate_single_SOR_method(next(input_generator))]

        # 9.8: Convert mp_data into two lists containing the k_ijs and all the diffusion_coefficients
        all_probabilities                                   = [probabilities                  for probabilities, diffusion_coefficient, components_of_diffusion_tensor, k_ijs_details, hopping_distances in mp_data]
        all_diffusion_coefficient_simulations               = [diffusion_coefficient          for probabilities, diffusion_coefficient, components_of_diffusion_tensor, k_ijs_details, hopping_distances in mp_data]
        all_components_of_diffusion_tensors_simulations     = [components_of_diffusion_tensor for probabilities, diffusion_coefficient, components_of_diffusion_tensor, k_ijs_details, hopping_distances in mp_data]
        molecule_k_ijs_details                              = [k_ijs_details                  for probabilities, diffusion_coefficient, components_of_diffusion_tensor, k_ijs_details, hopping_distances in mp_data]
        molecule_average_hopping_distances                  = [hopping_distances              for probabilities, diffusion_coefficient, components_of_diffusion_tensor, k_ijs_details, hopping_distances in mp_data]

        # 9.9: Determine the average porobability for moving between two particular molecules in the crystal across all the simulations you performed.
        average_probs = get_average_probs(all_probabilities)
        for key, all_values in average_probs.items():
            exciton_donor_to_acceptor_probabilities[key] = sum(average_probs[key])

        # 9.10: Obtain the average diffusion coefficient across all the simulations run.
        average_diffusion_coefficient = sum(all_diffusion_coefficient_simulations)/float(len(all_diffusion_coefficient_simulations))

        # 9.11: Save the diffusion coefficient into the diffusion_coefficient_molecules dictionary.
        diffusion_coefficient_molecules[mol1] = average_diffusion_coefficient

        # 9.12: Get the average diffusion tensor all the simulations run.
        average_components_of_diffusion_tensor = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        for components_diffusion_tensors in all_components_of_diffusion_tensors_simulations:
            for index in range(6):
                average_components_of_diffusion_tensor[index] += components_diffusion_tensors[index]
        for index in range(6):
            average_components_of_diffusion_tensor[index] /= len(all_components_of_diffusion_tensors_simulations)
        diffusion_tensor_molecules[mol1] = average_components_of_diffusion_tensor

        ################################################################################################################
        # This part indicates if you are running more than one simulation. 
        # This part of the code still needs to be written. Here are safe stops the code will not run until 
        # You sort out what happens if simulations > 1

        # Here for doing checks, remove and figure out how to fix this afterwards, GRW. 
        if len(molecule_k_ijs_details) > 1:
            raise Exception('Error')
        all_molecule_k_ijs_details[mol1] = molecule_k_ijs_details[0]

        if len(molecule_average_hopping_distances) > 1:
            raise Exception('Error')
        all_average_hopping_distances[mol1] = molecule_average_hopping_distances[0]
        ################################################################################################################

    # ---------------------------------------------------------------------------------------------------------------------------------------
    # Tenth, use the probabilities of moving between structurally unique molecules in exciton_donor_to_acceptor_probabilities to 
    # calculate the steady state probabilities for an exciton to be in on each structurally unique molecule in the crystal. 
    
    # 10.1: Create a dictionary for converting the names of the molecules to indices in an array
    mol1_to_index_converter = {mol1: index for mol1, index in zip(coupling_value_data.keys(), range(len(coupling_value_data)))}
    index_to_mol1_converter = list(coupling_value_data.keys())

    # 10.2: Create a matrix for holding the probabilities for an exciton to transfer between structurally unique molecule in the crystal. 
    exciton_donor_to_acceptor_probability_matrix = np.empty((len(mol1_to_index_converter),len(mol1_to_index_converter)))
    for (mol1, mol2), probability in exciton_donor_to_acceptor_probabilities.items():
        mol1_index = mol1_to_index_converter[mol1]
        mol2_index = mol1_to_index_converter[mol2]
        exciton_donor_to_acceptor_probability_matrix[mol2_index][mol1_index] = probability

    # 10.3: Solve the steady state solution. This will give the steady state population probabilities for an exciton being on each 
    #       structurally unique molecule at any particular time. 
    #       This involves:
    #           1: Solve the eigenequation for exciton_donor_to_acceptor_probability_matrix
    #           2: Check that one of the eigenvalues is 1.0.
    #           3: Select the eivenvector corresponding with the eigenvalue of 1.0.
    eigenvalues, eigenvectors = np.linalg.eig(exciton_donor_to_acceptor_probability_matrix)
    eigenvalues = [round(eigenvalue,8) for eigenvalue in eigenvalues]
    try:
        eigenvalue_1_index = eigenvalues.index(1.0)
    except ValueError as error:
        raise Exception('Error: Did not obtain an eigenvalue of 1.0 for your probability eigenfunction calculation.')
    probabilities = eigenvectors[:,eigenvalue_1_index].tolist()
    probabilities = [probability/sum(probabilities) for probability in probabilities]
    #print('probabilities: '+str(probabilities))

    # 10.4: Assign probabilities to molecules.
    probabilities = {index_to_mol1_converter[index]: probabilities[index] for index in range(len(probabilities))}
    # ---------------------------------------------------------------------------------------------------------------------------------------

    # Eleventh, obtain the overall exciton diffusion coefficient based on 
    overall_diffusion_coefficient = 0.0
    for key, probability in probabilities.items():
        overall_diffusion_coefficient += probability * diffusion_coefficient_molecules[key]

    # Twelfth, obtain the overall exciton diffusion tensor based on 
    overall_diffusion_tensor_molecules = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    for index in range(len(overall_diffusion_tensor_molecules)):
        for key, probability in probabilities.items():
            overall_diffusion_tensor_molecules[index] += probability * diffusion_tensor_molecules[key][index]

    # Thirteenth, use the components of overall_diffusion_tensor_molecules to obtain the eigenvalues and eigenvectors
    # of the diffusion tensor, and use the eigenvalues to obtain the diffusion coefficient from the diffusion tensor
    # that includes anisotropy effects. 
    diffusion_tensor_row_1 = (overall_diffusion_tensor_molecules[0], overall_diffusion_tensor_molecules[3], overall_diffusion_tensor_molecules[4])
    diffusion_tensor_row_2 = (overall_diffusion_tensor_molecules[3], overall_diffusion_tensor_molecules[1], overall_diffusion_tensor_molecules[5])
    diffusion_tensor_row_3 = (overall_diffusion_tensor_molecules[4], overall_diffusion_tensor_molecules[5], overall_diffusion_tensor_molecules[2])
    diffusion_tensor = np.array((diffusion_tensor_row_1, diffusion_tensor_row_2, diffusion_tensor_row_3))
    diffusion_tensor_eigenvalues, diffusion_tensor_eigenvectors = np.linalg.eig(diffusion_tensor)
    #overall_diffusion_coefficient_from_diffusion_tensor = (diffusion_tensor_eigenvalues[0] * diffusion_tensor_eigenvalues[1] * diffusion_tensor_eigenvalues[2]) ** (1.0/3.0)
    overall_diffusion_coefficient_from_diffusion_tensor = (diffusion_tensor_eigenvalues[0] + diffusion_tensor_eigenvalues[1] + diffusion_tensor_eigenvalues[2]) / 3.0

    # Fourteenth, return diffusion_coefficient_molecules
    return diffusion_coefficient_molecules, overall_diffusion_coefficient, overall_diffusion_coefficient_from_diffusion_tensor, all_molecule_k_ijs_details, all_average_hopping_distances

# ================================================================================================================================================================

def get_average_probs(all_probabilities):
    """
    This method is designed to take all he probabilities of moving between exciton donor and acceptors for all the simulations performed, and 
    give the total probability for moving between structurally unique molecules.

    Parameters
    ----------
    all_probabilities : dict
        This dictionary contains all the probabilities for an exciton to move from a donor to an acceptor for all simulations performed. 

    Returns
    -------
    average_probs : dict
        These are all the average probabilities for an exciton to move from a donor to an acceptor for all simulations performed. 
    """

    # First, 
    average_probs = {key: [[] for _ in range(len(value))] for key, value in all_probabilities[0].items()}
    for probabilities in all_probabilities:
        for key, value in probabilities.items():
            for index in range(len(value)):
                average_probs[key][index].append(value[index])
    for key, all_values in average_probs.items():
        for index in range(len(all_values)):
            values = all_values[index]
            average_probs[key][index] = sum(values)/float(len(values))

    return average_probs

# ================================================================================================================================================================
