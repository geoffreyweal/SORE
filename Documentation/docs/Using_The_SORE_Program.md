# How To Use The SORE Program

SORE is a python library rather than a terminal command: you drive it from a ``Run_SORE.py`` script. Unlike EKMC, SORE does not write slurm submission files — it performs the calculation in the process you run it in.

## A basic ``Run_SORE.py`` script

```python title="Run_SORE.py"
from SORE import Run_SORE

# EKMC_settings is the same settings dictionary that the EKMC program uses.
# See the EKMC documentation for a description of each of its entries.
EKMC_settings = {
	'folder_name':                               'my_crystal',
	'molecules_path':                            'path/to/ECCP/output/for/my_crystal',
	'functional_and_basis_set':                  'F_wB97XD_B_6_31plusGd_p',
	'kinetic_model':                             'Marcus',
	'short_range_couplings':                     ...,
	'long_range_couplings':                      ...,
	'kinetics_details':                          {...},
	'reorganisation_and_bandgap_energy_details': {...},
}

# SORE_settings describes how you want the sum-over-rates calculation performed.
SORE_settings = {'mode': 'standard'}

Run_SORE(EKMC_settings, SORE_settings=SORE_settings, no_of_cpus_for_setup=1)
```

## The arguments of ``Run_SORE``

* ``EKMC_settings`` (*dict.*): The crystal, coupling and kinetics settings, in exactly the same format that the [EKMC program](https://geoffreyweal.github.io/EKMC/Using_The_EKMC_Program) uses. SORE reads the same ECCP output that EKMC does, so you can run both on a crystal without rewriting your settings.
* ``SORE_settings`` (*dict.*): The settings specific to the sum-over-rates calculation. See below.
* ``save_initial_data`` (*bool.*): If ``True``, save the crystal neighbourhood and coupling data obtained during setup so that a later run can reuse it rather than recomputing it. Useful when you want to try several ``SORE_settings`` on the same crystal. Default: ``False``.
* ``no_of_cpus_for_setup`` (*int.*): The number of cpus to use while building the crystal neighbourhood. This is the slow part of the calculation, so raise it for large crystals. Default: ``1``.

## The ``SORE_settings`` dictionary

### ``mode`` (required)

``SORE_settings['mode']`` must be set. It chooses how the sum-over-rates equation is evaluated:

* ``'standard'``: Sum over the rate constants directly, with no energetic disorder included.
* ``'analytical_energetic_disorder'``: Include energetic disorder between molecules, treated analytically.
* ``'numeric_energetic_disorder'``: Include energetic disorder between molecules, treated numerically.

!!! note

	If ``mode`` is missing, or is not one of the three values above, SORE will report the problem and stop before performing any calculation.

### ``calculate_diffusion_tensor``

* ``calculate_diffusion_tensor`` (*bool.*): If ``True``, also obtain the exciton diffusion coefficient from the diffusion tensor, as well as the overall diffusion coefficient.

## Output from the SORE Program

``Run_SORE`` obtains:

* the exciton diffusion coefficient for each molecule in the crystal,
* the overall exciton diffusion coefficient for the crystal,
* the overall exciton diffusion coefficient obtained from the diffusion tensor (where requested),
* the rate constants (``k_ij``) between neighbouring molecules, and
* the average hopping distances.

These are written to an excel file so you can inspect the per-molecule and per-pair values that the overall diffusion coefficient was obtained from.

## Choosing between SORE and EKMC

Both programs obtain an exciton diffusion coefficient from the same ECCP data.

| | SORE | [EKMC](https://geoffreyweal.github.io/EKMC) |
| --- | --- | --- |
| Method | Sum-over-rates equation | Kinetic Monte Carlo simulation |
| Cost | Fast; runs in one process | Slow; many repeat simulations on slurm |
| Gives you | The diffusion coefficient | The diffusion coefficient, plus explicit exciton trajectories |

SORE is a good first look at a crystal, and a good way to scan many crystals or many parameter sets. Use EKMC when you need the trajectories themselves, or want to check the sum-over-rates result against an explicit simulation.
