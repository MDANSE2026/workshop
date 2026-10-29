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
PLATFORM._application_directory = Path.cwd()

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

import random

# Create the custom atom types
if "Atm1" not in ATOMS_DATABASE.atoms:
    ATOMS_DATABASE.add_atom("Atm1")
if "Atm2" not in ATOMS_DATABASE.atoms:
    ATOMS_DATABASE.add_atom("Atm2")
ATOMS_DATABASE.set_value("Atm1", "b_coherent", 1.0)
ATOMS_DATABASE.set_value("Atm1", "b_incoherent", 0.0)
ATOMS_DATABASE.set_value("Atm1", "color", "255;0;0")
ATOMS_DATABASE.set_value("Atm1", "vdw_radius", 0.215)
ATOMS_DATABASE.set_value("Atm2", "b_coherent", 10.0)
ATOMS_DATABASE.set_value("Atm2", "b_incoherent", 0.0)
ATOMS_DATABASE.set_value("Atm2", "color", "247;0;255")
ATOMS_DATABASE.set_value("Atm2", "vdw_radius", 0.15)
# uncomment the line below if you wish to save Atm1 and Atm2 to your
# local atom database
ATOMS_DATABASE.save()


N_ATOMS = 2048  # this should match the number of atoms in your trajectory.

isotopes = ["Atm1"] * int(N_ATOMS/2) + ["Atm2"] * int(N_ATOMS/2)

mapping = {}
for i in range(N_ATOMS):
    random_element = random.choice(isotopes)
    isotopes.remove(random_element)
    mapping[str(i)] = random_element


# setup to parameters to run the trajectory editor.
parameters = {
    'atom_charges': '{}',                               # atom_charges
    'atom_selection': '{"0": {"function_name": "select_all", "operation_type": "union"}}',                             # atom_selection
    'atom_transmutation': str(mapping).replace("'", '"'),  # atom_transmutation
    'frames': [0, 1001, 1],                             # frames
    'molecule_tolerance': [False, 0.04],                # molecule_tolerance
    'output_files': ('../mdanse_outputs/converted_atm1_atm2', 64, (1, 256), 'gzip', 'INFO', 16384),  # MDANSE trajectory (filename, format)
    'trajectory': ('../../tutorial_QENS/inputs/ar36_120fs_8_f64.mdt',
        'default',
        103809024,
        8191,
        1.0),  # Input trajectory file,  # trajectory
}


if __name__ == "__main__":
    trajectoryeditor = IJob.create('TrajectoryEditor')
    trajectoryeditor.run(parameters, status=True, prog_bar=True)
