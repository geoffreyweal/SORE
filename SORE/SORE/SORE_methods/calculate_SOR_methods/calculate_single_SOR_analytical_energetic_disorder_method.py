"""
calculate_single_SOR_analytical_energetic_disorder_method.py, Geoffrey Weal, 6/11/2023

This method is designed to obtain the Sum-Over_Rate equation, where energetic disroder has been included analytically. 
"""

import math
import numpy as np
from SORE.SORE.SORE_methods.calculate_SOR_methods.auxillary_methods import h_bar, k_b, Angstrom_to_cm

diffusion_tensor_multiplication = ((0,0), (1,1), (2,2), (0,1), (0,2), (1,2))
def calculate_single_SOR_analytical_energetic_disorder_method(data):
    """
    This program is designed to simulate the movement of an exciton through a OPV crystal system using the sum-over-rates method. 

    This method will obtain a single sum-over-rates calculation. 

    Parameters
    ----------
    mol1 : int
        This is the name of molecule 1.
    mol1_data : dict
        This dictionary contains the kinetic information of dimers that involve molecule 1 in the crystal (within some distance r_cut).
    mol1_com : tuple of (int, int, int)
        This tuple is the centre of mass coordinates for molecule 1. 
    unique_molname_1 : int
        This is the name of the molecule that molecule 1 is unique to. This will be different to mol1 if molecule 1 is structurally different to unique_molname_1. 
    bandgap_energy_mol1 : float
        This is the bandgap energy for molecule 1, given in eV.
    centre_of_masses : dict.
        This dictionary contains the centre of masses for the unique molecule in the crystal. 
    conformationally_equivalent_molecules : dict.
        This dictionary contains information about which molecules are conformationally equivalent to each other. 
    molecule_bandgap_energy_data : dict.
        These are the bandgap energies of the unique molecules in the crystal.
    energetic_disorder : float
        This is the energetic disorder value for the crystal, given in eV. 
    dimer_reorganisation_energy_data : dict.
        These are the reorganisation energies of the unique dimers in the crystal.
    crystal_cell_lattice_matrix : np.array
        This array describes the crystal lattice parameters.
    M_constant : float
        This is a constant given for the Marcus rate equation, which is specific for each dimer.
    X_constant : float
        This is a constant given for the Marcus rate equation, which is specific for each dimer.
    include_disorder_in_exciton_acceptor : bool.
        This boolean indicates if you would like to calculate the sum-over-rates method for multiple runs where the exciton acceptors (molecule 2s) contain disorder. The output given will be the average exciton diffusion as calculated by the SOR method. 
    """

    # Preamble, get all the information needed from data input.
    mol1, mol1_data, mol1_com, unique_molname_1, bandgap_energy_mol1, centre_of_masses, conformationally_equivalent_molecules, molecule_bandgap_energy_data, energetic_disorder, dimer_reorganisation_energy_data, crystal_cell_lattice_matrix, kinetic_model, temperature = data
    kT = k_b * temperature

    # First, initialise the diffusion coefficient for mol1 in the unit cell, as well as variables for storing k_ij vaules (which will be ued to calculate probabilities of movements).
    k_ijs = {}
    dhs   = []
    total_k_ij = 0.0
    diffusion_coefficient = 0.0
    components_of_diffusion_tensor = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    k_ijs_details = {}

    # Second, for the second molecule in the dimer:
    for mol2, m12_data in mol1_data.items():

        # 2.1: Obtain the centre of mass for molecule 2.
        mol2_com = centre_of_masses[mol2]

        # 2.2: Get the unique molecule identifier for molecule 2.
        unique_molname_2 = conformationally_equivalent_molecules[mol2] if (mol2 in conformationally_equivalent_molecules) else mol2

        # 2.3: Obtain the bandgap energy for molecule 2.
        bandgap_energy_mol2 = molecule_bandgap_energy_data[unique_molname_2]

        # 2.5: Obtain the difference in bamdgap energy between molecule 2 and molecule 1 in this dimer pair.
        difference_in_bandgap_energy = bandgap_energy_mol2 - bandgap_energy_mol1

        # 2.6: Get the dimer pair key.
        dimer_key = (unique_molname_1, unique_molname_2) # if (unique_molname_1 <= unique_molname_2) else (unique_molname_2, unique_molname_1)

        # 2.7: Get the reorganisation energy for this dimer pair
        reorganisation_energy = dimer_reorganisation_energy_data[dimer_key]

        # 2.8: Get the I term, that is obtained from reorganisation_energy, energetic_disorder, and kT
        I_term = get_I_term(reorganisation_energy, energetic_disorder, kT)

        # 2.9: For each displacement of molecule 2 by cell_displacement_ijk:
        for cell_displacement_ijk, coupling_value in m12_data.items():

            # 2.9.1: Obtain the displacement of molecule 2 into unit cell ijk
            cell_displacement = np.dot(np.array(cell_displacement_ijk), crystal_cell_lattice_matrix)

            # 2.9.2: Get the displacement between molecules in the dimer from molecule 1 to molecule 2
            dimer_displacement = (mol2_com + cell_displacement) - mol1_com

            # 2.9.3: Get the squared distance between molecule 1 and molecule 2
            distance_squared = np.dot(dimer_displacement,dimer_displacement)

            # 2.9.4: Get the rate constant based on the type of theory you want to use. 
            if kinetic_model.lower() == 'marcus':
                k_ij = (abs(coupling_value) ** 2.0) * I_term

            # 2.9.5: Save the k_ij data for each exciton acceptor environment 
            k_ijs.setdefault((mol1, mol2), []).append(k_ij)
            dhs  .append((distance_squared ** 0.5, k_ij))
            k_ijs_details[(mol1, mol2, cell_displacement_ijk[0], cell_displacement_ijk[1], cell_displacement_ijk[2])] = k_ij

            # 2.9.6: Add k_ij to the total k_ij for an exciton coming from the mol1 exciton donor.
            total_k_ij += k_ij

            # 2.9.7: Add this dimer's contribution to the diffusion coefficient.
            diffusion_coefficient += k_ij * distance_squared 

            for index, (index1, index2) in enumerate(diffusion_tensor_multiplication):
                components_of_diffusion_tensor[index] += k_ij * dimer_displacement[index1] * dimer_displacement[index2]

    if (total_k_ij == 0.0):
        raise Exception('Error: All rate constants were 0.0 s-1?. k_ijs: '+str(k_ijs))

    # Third, calculate the probabilities for take a certain step to. 
    probabilities = {key: [(value/total_k_ij) for value in values] for key, values in k_ijs.items()}

    # Fourth, obtain the average hopping distance from each exciton donor in the crystal.
    average_hopping_distance = 0
    for hopping_distance, k_ij in dhs:
        average_hopping_distance += hopping_distance * (k_ij/total_k_ij)

    # Fourth, multiply the diffusion_coefficient sum by 1/6 to convert it into a diffusion coefficient in 3D, and convert distance used from Angstroms to centimeters.
    diffusion_coefficient *= (1.0/6.0) * (Angstrom_to_cm ** 2.0)
    #diffusion_coefficient *= (1.0/6.0) * (1/total_k_ij) * (Angstrom_to_cm ** 2.0)

    # Fifth, multiply the components of diffusion_tensor by 1/2 to convert it into a diffusion coefficient in 1D, and convert distance used from Angstroms to centimeters.
    for index in range(len(components_of_diffusion_tensor)):
        components_of_diffusion_tensor[index] *= (1.0/2.0) * (Angstrom_to_cm ** 2.0)

    # Sixth, return the diffusion coefficient.
    return probabilities, diffusion_coefficient, components_of_diffusion_tensor, k_ijs_details, average_hopping_distance

# ================================================================================================================================================================

def get_I_term(reorganisation_energy, energetic_disorder, kT):
	"""
	This method is designed to obtain the I term in the SOR equation with disorder included analytically. 

	Parameters
	----------
	reorganisation_energy : float
		This is the reorganisation energy between molecule 1 and molecule 2 in the dimer. Given in eV
	energetic_disorder : float
		This is the energetic disorder value for the crystal, given in eV. 
	kT : float
		This is k_b * temperature (boltzmann constant times by temeprature). given in eV. 

	Returns
	-------
	k_ij (float): This is the rate constant for this dimer.
	"""

    # First, determine the denominator term that arrises in several parts of the equation.
	denominator_term      = 2.0*(energetic_disorder**2.0) + (4.0*reorganisation_energy*kT)

	# Second, calculate the prefactor term.
	prefactor_numerator   = 4 * math.pi
	prefactor_term        = math.sqrt(prefactor_numerator/denominator_term)

	# Third, calculation the exponential term. 
	exponential_numerator = ( reorganisation_energy + ((energetic_disorder ** 2.0)/kT) ) ** 2.0
	exponential_term      = math.exp(-exponential_numerator/denominator_term)

	# Fourth, determine the I term.
	I_term                = (1/h_bar) * prefactor_term * exponential_term

	# Fifth, return I_term 
	return I_term

# ================================================================================================================================================================
