'''
Run_SORE.py, Geoffrey Weal, 3/11/2023

This program will calculate the exciton diffusion coefficients for OPV materials using the sum-over-rates method. 
'''
import xlsxwriter
import os, tqdm, math
from copy import deepcopy
from SORE.SORE.SORE_methods.get_SORE_version                          import get_SORE_version
from EKMC.EKMC_Setup.EKMC_Only_Setup.EKMC_Only_Setup                  import initial_setup, get_non_changing_lattice_kinetics_details
from EKMC.EKMC_Setup.EKMC_Only_Setup.get_constant_rate_law_data       import get_constant_rate_law_data
from SORE.SORE.SORE_methods.SORE_input_data_methods                   import save_SORE_input_data, read_SORE_input_data
from SORE.SORE.SORE_methods.get_conformationally_equivalent_molecules import get_conformationally_equivalent_molecules
from SORE.SORE.SORE_methods.get_RE_and_bandgap_data                   import get_RE_and_bandgap_data
from SORE.SORE.SORE_methods.calculate_SOR                             import calculate_SOR

def Run_SORE(EKMC_settings, SORE_settings={}, save_initial_data=False, no_of_cpus_for_setup=1):
    """
    This method is designed to setup the exciton kMC algorithm for performing multiple repeats of the same simulation in slurm.

    Parameters
    ----------
    EKMC_settings : dict.
        This contains the information about the crystal you want to simulate using this kMC program.
    no_of_cpus_for_setup : int.
        This is the number of cpus used to setup the EKMC simulations.
    """

    include_solvents = False # to deal withat some point

    # First, get the file names and paths of data to obtain local neighbourhood data for
    folder_name      = EKMC_settings['folder_name']
    molecules_path   = EKMC_settings['molecules_path']
    crystal_name     = os.path.basename(molecules_path)

    # Second, obtain the functional and basis set settings used for this simulation
    functional_and_basis_set = EKMC_settings['functional_and_basis_set']

    # Third, obtain the kinetic model that you would like to use for performing the kinetic monte carlo algorithm.
    kinetic_model = EKMC_settings['kinetic_model']

    # Fourth, obtain information about the short-range and long-range couplings from EKMC_settings
    short_range_couplings = EKMC_settings['short_range_couplings']
    long_range_couplings  = EKMC_settings['long_range_couplings']

    # Fifth, obtain information about rCut for short-range and long-range coupling models. 
    short_range_rCut          = EKMC_settings['short_range_rCut']
    long_range_rCut           = EKMC_settings.get('long_range_rCut',None)
    rCut_mol_dist_description = EKMC_settings['rCut_mol_dist_description']

    # Sixth, obtain information for obtaining reorganisation energies and bandgap energies 
    # for an excton moving from one molecule to another in the crystal.
    reorganisation_and_bandgap_energy_details = EKMC_settings['reorganisation_and_bandgap_energy_details']

    # Seventh, obtain information about the kinetic details and neighbourhood cut-off from EKMC_settings
    kinetics_details = EKMC_settings['kinetics_details']

    # Eighth, determine the type of SOR calculation you want to perform
    if not ('mode' in SORE_settings):
        print("Error: SORE_settings['mode'] must be set to either:")
        print('  * standard:                      Calculate SOR using Marcus Theory.')
        print('  * numeric_energetic_disorder:    Calculate SOR using Marcus Theory, where site energy disorder in is included numerically (a number of repeated simulations are required).')
        print('  * analytical_energetic_disorder: Calculate SOR using Marcus Theory, where site energy disorder in is included analytically.')
        print("Include SORE_settings['mode'] in your Run_SORE.py script and rerun.")
        exit('This program will exit without beginning.')
    if   SORE_settings['mode'].lower() == 'standard':
        SORE_mode = 'standard'
    elif SORE_settings['mode'].lower() == 'numeric_energetic_disorder':
        SORE_mode = 'numeric_energetic_disorder'
    elif SORE_settings['mode'].lower() == 'analytical_energetic_disorder':
        SORE_mode = 'analytical_energetic_disorder'
    else:
        print("Error: SORE_settings['mode'] must be set to either:")
        print('  * standard:                      Calculate SOR using Marcus Theory.')
        print('  * numeric_energetic_disorder:    Calculate SOR using Marcus Theory, where site energy disorder in is included numerically (a number of repeated simulations are required).')
        print('  * analytical_energetic_disorder: Calculate SOR using Marcus Theory, where site energy disorder in is included analytically.')
        print("Current SORE_settings['mode'] = "+str(SORE_settings['mode']))
        print("Include SORE_settings['mode'] in your Run_SORE.py script and rerun.")
        exit('This program will exit without beginning.')

    '''
    # Ninth, determine if you want to calculate the Sum-over-rates diffusion coefficient by either summing all scalar diffusion coefficient, or instead calculating the diffusion tensor. 
    if 'calculate_diffusion_tensor' in SORE_settings:
        calculate_diffusion_tensor = SORE_settings['calculate_diffusion_tensor']
        if not isinstance(calculate_diffusion_tensor,bool):
            raise Exception("Error: SORE_settings['calculate_diffusion_tensor'] must be a boolean.\nCheck this.\nSORE_settings['calculate_diffusion_tensor'] = "+str(SORE_settings['calculate_diffusion_tensor']))
    else:
        calculate_diffusion_tensor = False
    '''

    # Tenth, print information about the SORE program.
    print('######################################################################')
    print('######################################################################')
    print('######################################################################')
    print('---------------------------- SORE Program ----------------------------')
    version_no = get_SORE_version()
    version_string = 'Version '+str(version_no)
    space_no = math.floor((70 - len(version_string))/2.0)
    print(' '*space_no+str(version_string))
    title_name = 'Getting Sum-over-Rates for: '+str(crystal_name)
    space_no = math.floor((70 - len(title_name))/2.0)
    print(' '*space_no+str(title_name))
    title_name = 'No of CPUs: '+str(no_of_cpus_for_setup)
    space_no = math.floor((70 - len(title_name))/2.0)
    print(' '*space_no+str(title_name))
    print('######################################################################')
    print('######################################################################')
    print('######################################################################')

    # Eleventh, get the path to store this simulation in.
    SORE_Data_name = 'SORE_Data'
    if ('overall_folder_suffix_name' in EKMC_settings):
        overall_folder_suffix_name = EKMC_settings['overall_folder_suffix_name']
        if not (overall_folder_suffix_name == '' or overall_folder_suffix_name is None):
            SORE_Data_name += '_'+str(overall_folder_suffix_name)
    path_to_save_data_to = os.getcwd()+'/'+SORE_Data_name # +'/'+folder_name+'/'+functional_and_basis_set

    # ===========================================================================================================================================================================

    if save_initial_data:

        # Twelfth, check if the folder for saving initial setup data to exists, and if not make it. 
        SORE_input_data_folderpath = 'SORE_input_data'
        if not os.path.exists(SORE_input_data_folderpath):
            os.makedirs(SORE_input_data_folderpath)
        SORE_input_data_crystal_folderpath = SORE_input_data_folderpath+'/'+folder_name

        # Thirteenth, setup the kMC simulations and get the exciton kMC setup files for performing simulations.
        if not os.path.exists(SORE_input_data_crystal_folderpath):
            molecule_names, molecules, crystal_cell_lattice, all_coupling_values = initial_setup(molecules_path, functional_and_basis_set, kinetic_model, short_range_couplings, long_range_couplings, short_range_rCut, long_range_rCut, rCut_mol_dist_description, kinetics_details, include_solvents=include_solvents, no_of_cpus_for_setup=no_of_cpus_for_setup)
            save_SORE_input_data(SORE_input_data_crystal_folderpath, molecule_names, molecules, crystal_cell_lattice, all_coupling_values)
        else:
            molecule_names, molecules, crystal_cell_lattice, all_coupling_values = read_SORE_input_data(SORE_input_data_crystal_folderpath)

    else:

        molecule_names, molecules, crystal_cell_lattice, all_coupling_values = initial_setup(molecules_path, functional_and_basis_set, kinetic_model, short_range_couplings, long_range_couplings, short_range_rCut, long_range_rCut, rCut_mol_dist_description, kinetics_details, include_solvents=include_solvents, no_of_cpus_for_setup=no_of_cpus_for_setup)

    # ===========================================================================================================================================================================

    # Second, extract the non0change componenets of kinetics_details
    non_changing_lattice_kinetics_details = get_non_changing_lattice_kinetics_details(kinetic_model, kinetics_details)

    # Third, obtain the constant rate law data for the rate law you want to use. 
    constant_rate_data = get_constant_rate_law_data(kinetic_model, non_changing_lattice_kinetics_details)

    # Fourteenth, obtain the centre of masses for the molecules in the crystal unit cell. 
    molecule_centres_of_mass = {int(molecule_name): molecule.get_center_of_mass() for molecule_name, molecule in zip(molecule_names, molecules)}

    # Fiftheenth, obtain all the conformationally equivalent molecules in the crystal.
    conformationally_equivalent_molecules, conformationally_unique_molecule_names = get_conformationally_equivalent_molecules(molecules_path, molecule_names)

    # Sixteenth, obtain the bandgap for molecules (for dimers) and reorganisation energies (for dimers) for this crystal. 
    molecule_bandgap_energy_data, dimer_reorganisation_energy_data = get_RE_and_bandgap_data(reorganisation_and_bandgap_energy_details, crystal_name, functional_and_basis_set, conformationally_equivalent_molecules, molecule_names, conformationally_unique_molecule_names)

    # Seventeenth, calculate the Sum-Over-Rates values for your crystal with the desired settings. 
    diffusion_coefficient_molecules, overall_diffusion_coefficient, overall_diffusion_coefficient_from_diffusion_tensor, k_ijs_details, average_hopping_distances = calculate_SOR(SORE_mode, molecule_names, molecules, crystal_cell_lattice, kinetic_model, non_changing_lattice_kinetics_details, molecule_bandgap_energy_data, dimer_reorganisation_energy_data, conformationally_equivalent_molecules, all_coupling_values, constant_rate_data, no_of_cpus=no_of_cpus_for_setup)

    # Eighteenth, return diffusion_coefficient_molecules, overall_diffusion_coefficient, k_ijs_details, and average_hopping_distances
    return diffusion_coefficient_molecules, overall_diffusion_coefficient, overall_diffusion_coefficient_from_diffusion_tensor, k_ijs_details, average_hopping_distances
