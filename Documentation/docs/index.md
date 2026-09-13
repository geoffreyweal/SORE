# The Sum-Over-Rates for Excitons (SORE) Program

Authors: Dr. Geoffrey Weal<sup>\*,†</sup>, Dr. Chayanit Wechwithayakhlung<sup>†</sup>, Assoc. Prof. Daniel Packwood<sup>†</sup>, Dr. Paul Hume<sup>\*</sup>, Prof. Justin Hodgkiss<sup>\*</sup>

<sup>\*</sup> Victoria University of Wellington, Wellington, New Zealand; The MacDiarmid Institute for Advanced Materials and Nanotechnology, Wellington, New Zealand. 

<sup>†</sup> Institute for Integrated Cell-Material Sciences (iCeMS), Kyoto University, Kyoto, Japan.

Group pages: https://people.wgtn.ac.nz/paul.hume/grants, https://www.packwood.icems.kyoto-u.ac.jp/, https://people.wgtn.ac.nz/justin.hodgkiss/grants


## What is the Sum-Over-Rates for Excitons (SORE) Program

The Sum-Over-Rates for Excitons (SORE) program calculates the exciton diffusion coefficient of a molecular crystal using a sum-over-rates equation.

SORE answers the same question as the [EKMC program](https://geoffreyweal.github.io/EKMC) — how fast does an exciton diffuse through this crystal — but does so analytically rather than by simulating individual exciton trajectories. Where EKMC runs many kinetic Monte Carlo simulations and averages them, SORE sums the hopping rate constants directly. This makes SORE much faster, at the cost of the detail that an explicit trajectory gives you.

SORE uses the same electronic data as EKMC: excited-state energies, reorganisation energies, and exciton (EET) couplings obtained from DFT by the [ECCP program](https://geoffreyweal.github.io/ECCP). It also takes the same ``EKMC_settings`` dictionary that you would give to EKMC, so you can run both on the same crystal without rewriting your settings.

SORE can be run in three modes, set by ``SORE_settings['mode']``:

* ``standard``: Sum over the rate constants directly, with no energetic disorder.
* ``analytical_energetic_disorder``: Include energetic disorder, treated analytically.
* ``numeric_energetic_disorder``: Include energetic disorder, treated numerically.

## Installation

It is recommended to read the installation page before using the SORE program. See [Installation: Setting Up SORE and Pre-Requisites Packages](Installation.md) for more information.

## Guide To Using SORE

The SORE program is one in a series of programs that are designed to be used in the workflow shown below. After you have installed SORE, see [How To Use The SORE Program](Using_The_SORE_Program.md) to learn about how to use this program.

## The Grand Scheme

The SORE program is used as part of a grand scheme for calculating the excited-state electronic properties of molecules in a crystal. This includes simulations of exciton and charge diffusion through crystal structures, in particular for organic molecules (but not limited to them). This scheme is shown below, along with where the SORE program is used in this scheme. 

<img alt="Schematic of Grand Scheme" src="Shared_Images/Grand_Scheme/Grand_Scheme.png?raw=true#only-light" />
<img alt="Schematic of Grand Scheme" src="Shared_Images/Grand_Scheme/Grand_Scheme_Dark.png?raw=true#only-dark" />

## Websites and Github Repositories for All Associated Programs

### Instructional Websites

* ACSD: https://geoffreyweal.github.io/ACSD
* ReCrystals: https://geoffreyweal.github.io/ReCrystals
* RSGC: https://geoffreyweal.github.io/RSGC
* ReJig: https://geoffreyweal.github.io/ReJig
* ECCP: https://geoffreyweal.github.io/ECCP
* EKMC: https://geoffreyweal.github.io/EKMC
* SORE: https://geoffreyweal.github.io/SORE
* SUMELF: https://geoffreyweal.github.io/SUMELF

### Github Repositories

* ACSD: https://github.com/geoffreyweal/ACSD
* ReCrystals: https://github.com/geoffreyweal/ReCrystals
* RSGC: https://github.com/geoffreyweal/RSGC
* ReJig: https://github.com/geoffreyweal/ReJig
* ECCP: https://github.com/geoffreyweal/ECCP
* EKMC: https://github.com/geoffreyweal/EKMC
* SORE: https://github.com/geoffreyweal/SORE
* SUMELF: https://github.com/geoffreyweal/SUMELF
