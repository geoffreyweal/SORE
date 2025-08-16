"""
calculate_single_SOR_standard_method.py, Geoffrey Weal, 6/11/23

This method is designed to calculate the SOR for your simulation using Marcus theory. 

You can include disorder in your simulations if desired. 
"""

import numpy as np
from math import exp
try:
    from numpy.random.Generator import normal
except:
    from numpy.random import normal
from SORE.SORE.SORE_methods.calculate_SOR_methods.auxillary_methods import Angstrom_to_cm

def calculate_single_SOR_standard_method(data):
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
    bandgap_energy_mol1_with_disorder : float
        This is the bandgap energy for molecule 1. This will be the thermalised energy if you have given an energetic disorder in your python running script. 
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

    # Preamble 1, there seems to be a weird seeding problem when using multiprocessing. 
    # We will reseed every time this method is used to make sure that everything is random. 
    np.random.seed()

    # Preamble 2, get all the information needed from data input.
    mol1, mol1_data, mol1_com, unique_molname_1, bandgap_energy_mol1_with_disorder, centre_of_masses, conformationally_equivalent_molecules, molecule_bandgap_energy_data, energetic_disorder, dimer_reorganisation_energy_data, crystal_cell_lattice_matrix, kinetic_model, M_constant, X_constant, include_disorder_in_exciton_acceptor = data

    # First, initialise the diffusion coefficient for mol1 in the unit cell, as well as variables for storing k_ij vaules (which will be ued to calculate probabilities of movements).
    k_ijs = {}
    total_k_ij = 0.0
    diffusion_coefficient = 0.0

    # Second, for the second molecule in the dimer:
    for mol2, m12_data in mol1_data.items():

        # 2.1: Obtain the centre of mass for molecule 2.
        mol2_com = centre_of_masses[mol2]

        # 2.2: Get the unique molecule identifier for molecule 2.
        unique_molname_2 = conformationally_equivalent_molecules[mol2] if (mol2 in conformationally_equivalent_molecules) else mol2

        # 2.3: Obtain the bandgap energy for molecule 2.
        bandgap_energy_mol2 = molecule_bandgap_energy_data[unique_molname_2]

        # 2.4: If you want to include disorder in the exciton acceptor, do this here.
        if include_disorder_in_exciton_acceptor:
            bandgap_energy_mol2_with_disorder = normal(bandgap_energy_mol2, energetic_disorder)
        else:
            bandgap_energy_mol2_with_disorder = bandgap_energy_mol2

        # 2.5: Obtain the difference in bamdgap energy between molecule 2 and molecule 1 in this dimer pair.
        difference_in_bandgap_energy = bandgap_energy_mol2_with_disorder - bandgap_energy_mol1_with_disorder

        # 2.6: Get the dimer pair key.
        dimer_key = (unique_molname_1, unique_molname_2) if (unique_molname_1 <= unique_molname_2) else (unique_molname_2, unique_molname_1)

        # 2.7: Get the reorganisation energy for this dimer pair
        reorganisation_energy = dimer_reorganisation_energy_data[dimer_key]

        # 2.8: For each displacement of molecule 2 by cell_displacement_ijk:
        for cell_displacement_ijk, coupling_value in m12_data.items():

            # 2.8.1: Obtain the displacement of molecule 2 into unit cell ijk
            cell_displacement = np.dot(np.array(cell_displacement_ijk), crystal_cell_lattice_matrix)

            # 2.8.2: Get the displacement between molecules in the dimer from molecule 1 to molecule 2
            dimer_displacement = (mol2_com + cell_displacement) - mol1_com

            # 2.8.3: Get the squared distance between molecule 1 and molecule 2
            distance_squared = np.dot(dimer_displacement,dimer_displacement)

            # 2.8.4: Get the rate constant based on the type of theory you want to use. 
            if kinetic_model.lower() == 'marcus':
                k_ij = get_rate_constant_Marcus(M_constant, X_constant, difference_in_bandgap_energy, coupling_value, reorganisation_energy)

            # 2.8.5: Save the k_ij data for each exciton acceptor environment 
            k_ijs.setdefault((mol1, mol2), []).append(k_ij)

            # 2.8.6: Add k_ij to the total k_ij for an exciton coming from the mol1 exciton donor.
            total_k_ij += k_ij

            # 2.8.7: Add this dimer's contribution to the diffusion coefficient.
            diffusion_coefficient += k_ij * distance_squared
            #diffusion_coefficient += k_ij * k_ij * distance_squared

    # Third, calculate the probabilities for take a certain step to. 
    probabilities = {key: [(value/total_k_ij) for value in values] for key, values in k_ijs.items()}

    # Fourth, multiply the diffusion_coefficient sum by 1/6 to convert it into a diffusion coefficient in 3D, and convert distance used from Angstroms to centimeters.
    diffusion_coefficient *= (1.0/6.0) * (Angstrom_to_cm ** 2.0)
    #diffusion_coefficient *= (1.0/6.0) * (1/total_k_ij) * (Angstrom_to_cm ** 2.0)

    # Fifth, return the diffusion coefficient.
    return probabilities, diffusion_coefficient

# ================================================================================================================================================================

def get_rate_constant_Marcus(M_constant, X_constant, difference_in_bandgap_energy, coupling_value, reorganisation_energy):
    """
    This method is designed to gather the rate constants between the current molecule in the current cell point and every other neighbouring molecule in surrounding unit cells within the rCut_neighbourhood radius. 

    This method is based on Marcus theory.

    Parameters
    ----------
    M_constant : float
        This is a constant containing 2pi/h_bar * (1/sqrt(4*pi*k*T)).
    X_constant : float
        This is a constant containing 1/(4*k*T).
    difference_in_bandgap_energy : float
        This is the difference in the bandgap energy between molecule 1 and molecule 2 in the dimer. Given in eV
    coupling_value : float
        This is the coupling value between molecule 1 and molecule 2 in the dimer. Given in eV
    reorganisation_energy : float
        This is the reorganisation energy between molecule 1 and molecule 2 in the dimer. Given in eV
    Returns
    -------
    k_ij (float): This is the rate constant for this dimer.
    """

    # First, obtain the prefix value for the Marcus theory rate constant equation.
    prefix_value = (abs(coupling_value) ** 2.0) / (reorganisation_energy ** 0.5)

    # First, obtain the exponential value for the Marcus theory rate constant equation.
    exp_value = ((difference_in_bandgap_energy + reorganisation_energy) ** 2.0) / reorganisation_energy

    # Third, obtain the rate constant for this dimer as based on Marocus theory.
    k_ij = prefix_value * M_constant * exp( -X_constant * exp_value )

    # Fourth, return the rate constant, k_ij.
    return k_ij

# ================================================================================================================================================================

