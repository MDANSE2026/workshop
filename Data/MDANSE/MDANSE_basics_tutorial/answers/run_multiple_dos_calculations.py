#!/Users/maciej.bartkowiak/python/builder/bin/python3.14

import os
os.environ.update(
    OMP_NUM_THREADS = '1',
    OPENBLAS_NUM_THREADS = '1',
    MKL_NUM_THREADS = '1',
    VECLIB_MAXIMUM_THREADS = '1',
    NUMEXPR_NUM_THREADS = '1'
)

########################################################
# This is an automatically generated MDANSE run script #
########################################################

from MDANSE.Framework.Jobs.IJob import IJob

########################################################
# Job parameters                                       #
########################################################

def run_dos(trajectory: str, projection_vector: tuple[str]):
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
        'frames': [0, 501, 1, 251],  # Frame sampling and correlation frames. The values are: first frame, last frame, step size.
        'grouping_level': 'atom',  # Grouping level. The partial results can be separated between different molecule types.
        'instrument_resolution': ('ideal', {}),  # Instrument resolution function
        'interpolation_order': 0,  # Interpolation order for velocity determination
        'output_files': (f'dos_{trajectory}_projection{",".join(projection_vector)}',
            ['MDAFormat'],
            'no logs'),  # Analysis result files (name and format). It is possible to use MDAFormat and TextFormat simultaneously.
        'projection': ('AxialProjector', projection_vector),  # Coordinate projection. The coordinates from the trajectory can be projected on a line or plane.
        'running_mode': ('single-core',),  # Parallelisation options. Run the job on a single core or parallelise over multiple cores.
        'trajectory': (f'{trajectory}.mdt',
            'default',
            103809024,
            8191,
            1.0),  # Input trajectory file
        'weights': 'atomic_weight',  # Weights. Atom property selected here will be used for calculating the scaling factors of the results.
    }
    densityofstates = IJob.create("DensityOfStates")
    # Progress bars only available if tqdm available.
    # Install with `cli` optional dependency.
    densityofstates.run(parameters, status=True, prog_bar=True)

########################################################
# Setup and run the analysis                           #
########################################################

if __name__ == "__main__":
    for traj_name in ["isotropic", "anisotropic"]:
        for projection in [("1","0","0"),
                           ("0","1","0"),
                           ("0","0","1"),]:
            run_dos(traj_name, projection)
