#!/usr/bin/env python

import os
from pathlib import Path
os.environ.update(
    OMP_NUM_THREADS = '1',
    OPENBLAS_NUM_THREADS = '1',
    MKL_NUM_THREADS = '1',
    VECLIB_MAXIMUM_THREADS = '1',
    NUMEXPR_NUM_THREADS = '1'
)


########################################################
# Initialise settings                                  #
########################################################

from MDANSE.Chemistry import ATOMS_DATABASE
from MDANSE.Core.Platform import PLATFORM
from MDANSE.Core.Settings import Settings
from MDANSE.Framework.Units import UNITS_MANAGER

# uncomment this line to use MDANSE settings from the current path
# instead of global settings from the home directory
# you can also introduce a specific path instead of Path.cwd()
# PLATFORM._application_directory = Path.cwd()

ATOMS_DATABASE._load()
UNITS_MANAGER._load()
Settings.init(PLATFORM.main_settings, save=False)

########################################################
# This is an automatically generated MDANSE run script #
########################################################

from MDANSE.Framework.Jobs.IJob import IJob

########################################################
# Job parameters                                       #
########################################################

parameters_atm1_atm2 = {
    'atom_transmutation': '{}',                         # atom_transmutation
    'dcsf_input_file': '../mdanse_outputs/dcsf_atm1_atm2.mda',  # MDANSE Coherent Structure Factor
    'disf_input_file': '../mdanse_outputs/disf_atm1_atm2.mda',  # MDANSE Incoherent Structure Factor
    'output_files': ('../mdanse_outputs/ndtsf_atm1_atm2.mda', ['MDAFormat'], 'INFO'),  # output_files
    'running_mode': ('single-core',),                   # running_mode
    'trajectory': '../mdanse_outputs/converted_atm1_atm2.mdt',  # trajectory
}
parameters_atm1atm2 = {
    'atom_transmutation': '{}',                         # atom_transmutation
    'dcsf_input_file': '../mdanse_outputs/dcsf_atm1atm2.mda',  # MDANSE Coherent Structure Factor
    'disf_input_file': '../mdanse_outputs/disf_atm1atm2.mda',  # MDANSE Incoherent Structure Factor
    'output_files': ('../mdanse_outputs/ndtsf_atm1atm2.mda', ['MDAFormat'], 'INFO'),  # output_files
    'running_mode': ('single-core',),                   # running_mode
    'trajectory': '../mdanse_outputs/converted_atm1atm2.mdt',  # trajectory
}

########################################################
# Setup and run the analysis                           #
########################################################

if __name__ == "__main__":
    neutrondynamictotalstructurefactor = IJob.create('NeutronDynamicTotalStructureFactor')
    neutrondynamictotalstructurefactor.run(parameters_atm1_atm2, status=True)
    neutrondynamictotalstructurefactor = IJob.create('NeutronDynamicTotalStructureFactor')
    neutrondynamictotalstructurefactor.run(parameters_atm1atm2, status=True)
