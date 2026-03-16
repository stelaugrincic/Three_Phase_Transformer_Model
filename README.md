DYNAMIC THREE-PHASE THREE-LIMB TRANSFORMER MODEL

Author: Stela Ugrincic

Affiliation: University of Maribor, Faculty of Electrical Engineering and Computer Science, Maribor, Slovenia

Email: stela.ugrincic1@um.si

Last modified: 16.03.2026

This repository contains the numerical model and source code of a three-phase three-limb transformer, for steady-state and low-frequency electromagnetic transient analyses.

## Overview

The dynamic transformer model is derived as a lumped-parameter model in the form of nonlinear ordinary differential equations.
It accounts for:
- Nonlinear magnetic properties of the transformer core through nonlinear characteristic;
- Different winding connections (Y, D, Yn, etc.);
- It can be modified for different core topologies.

The model implemented in the provided python scripts is especially suitable for steady-state, transformer inrush and short-circuit calculations.

## Repository Structure

- `transformer_data.py`: Rated parameters and construction data.
- `model_configuration.py`: Configuration of model parameters depending on the transformer topology and provided data.
- `transformer_model.py`: ODE system definitions.
- `main_transformer_model.py`: Main simulation and plotting script.
- `run_simulation.py`: Script for running simulations and deciding on saving and exporting configurations.

## How to Run - Requirements
1. Ensure you have Python 3.x installed.
2. Necessary libraries:
   ```bash
   pip install numpy scipy matplotlib pandas openpyxl

## How to Run - simulation
Run the simulation through the run_simulation.py script. To customize the run, you can uncomment the following lines in run_simulation.py:

    "--simple":        Uncomment to use simplified model, otherwise the extended (full) model is being used by default; 
                       the simplified model can be useful for inrush calculations (is=0), in case of no-load energization both models have the same results (validated)
    "--export-mat":    Uncomment to export simulation results to a MATLAB .mat file, other options: .csv and .xlsx
    "--save-plots":    Uncomment to save generated figures to the output directory
    "--no-show":       Uncomment to run the simulation without displaying plots on the screen

## Output
All results (saved plots, results) are saved in the directory specified by the "--output-dir" argument (default: ./results).   

## Default simulation:
No-load transformer energization at the zero crossing of the primary voltage of the first phase (at the time t=0.095 s, towards the negative maximum)

## Language
Please note that the source code was developed in a Slovenian research environment, therefore:
  - Variable names are mostly named in Slovenian,
  - Script names were translated into English for easier understanding,
  - Key sections of the code include comments in both English and Slovenian. 