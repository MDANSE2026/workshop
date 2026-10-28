#!/Users/maciej.bartkowiak/python/builder/bin/python3.14

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

parameters = {
    'atom_selection': """{
        "0": {
            "function_name": "select_all",
            "operation_type": "union"
        },
        "1": {
            "function_name": "select_dummy",
            "operation_type": "difference"
        }
    }""",  # Atom selection. The analysis can be run on a subset of atoms in the trajectory.
    'atom_transmutation': '{}',  # Atom transmutation. Atom types (and properties) can be temporarily changed for the current analysis run.
    'frames': [0, 1001, 1, 501],  # Frame sampling and correlation frames. The values are: first frame, last frame, step size.
    'grouping_level': 'atom',  # Grouping level. The partial results can be separated between different molecule types.
    'instrument_resolution': ('ideal', {}),  # Instrument resolution function
    'output_files': ('results/argon-disf',
        ['MDAFormat', 'TextFormat'],
        'no logs'),  # Analysis result files (name and format). It is possible to use MDAFormat and TextFormat simultaneously.
    'projection': ('NullProjector', []),  # Coordinate projection. The coordinates from the trajectory can be projected on a line or plane.
    'q_vectors': ('SphericalQVectors',
        {'n_vectors': 100,
        'seed': 0,
        'shells': [10.0, 44.0, 2.0],
        'width': 0.1}),  # Q vector generator. Determines how the q vectors will be generated and grouped into shells.
    'running_mode': ('single-core',),  # Parallelisation options. Run the job on a single core or parallelise over multiple cores.
    'trajectory': ('inputs/ar36_120fs_8_f64.mdt',
        'default',
        103809024,
        8191,
        1.0),  # Input trajectory file
    'weights': 'equal',  # Weights. Atom property selected here will be used for calculating the scaling factors of the results.
}

########################################################
# Setup and run the analysis                           #
########################################################

if __name__ == "__main__":
    dynamicincoherentstructurefactor = IJob.create("DynamicIncoherentStructureFactor")
    # Progress bars only available if tqdm available.
    # Install with `cli` optional dependency.
    dynamicincoherentstructurefactor.run(parameters, status=True, prog_bar=True)
