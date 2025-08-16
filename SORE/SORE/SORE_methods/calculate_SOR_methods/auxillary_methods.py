

Angstrom_to_meter = 10.0 ** -10.0
meter_to_cm       = 10.0 **   2.0
#Angstrom_to_cm    = Angstrom_to_meter * meter_to_cm
Angstrom_to_cm    = 1.0 * (10.0 ** -8.0)
h_bar             = 6.582119569    * (10.0 ** -16.0) # eV s
k_b               = 8.617333262145 * (10.0 ** -5.0) # eV K-1

def single_SOR_standard_method_input_generator(no_of_simulations, mol1, mol1_data, mol1_com, unique_molname_1, bandgap_energy_mol1_with_disorder, centre_of_masses, conformationally_equivalent_molecules, molecule_bandgap_energy_data, energetic_disorder, dimer_reorganisation_energy_data, crystal_cell_lattice_matrix, kinetic_model, M_constant, X_constant, include_disorder_in_exciton_acceptor):
    '''
    Generator for the calculate_single_SOR_simulation method. 

    See calculate_single_SOR_simulation for parameters meanings. 
    '''
    for counter in range(no_of_simulations):
        yield (mol1, mol1_data, mol1_com, unique_molname_1, bandgap_energy_mol1_with_disorder, centre_of_masses, conformationally_equivalent_molecules, molecule_bandgap_energy_data, energetic_disorder, dimer_reorganisation_energy_data, crystal_cell_lattice_matrix, kinetic_model, M_constant, X_constant, include_disorder_in_exciton_acceptor)

def single_SOR_analytical_energetic_disorder_method_input_generator(no_of_simulations, mol1, mol1_data, mol1_com, unique_molname_1, bandgap_energy_mol1_with_disorder, centre_of_masses, conformationally_equivalent_molecules, molecule_bandgap_energy_data, energetic_disorder, dimer_reorganisation_energy_data, crystal_cell_lattice_matrix, kinetic_model, temperature):
    '''
    Generator for the calculate_single_SOR_simulation method. 

    See calculate_single_SOR_simulation for parameters meanings. 
    '''
    for counter in range(no_of_simulations):
        yield (mol1, mol1_data, mol1_com, unique_molname_1, bandgap_energy_mol1_with_disorder, centre_of_masses, conformationally_equivalent_molecules, molecule_bandgap_energy_data, energetic_disorder, dimer_reorganisation_energy_data, crystal_cell_lattice_matrix, kinetic_model, temperature)
