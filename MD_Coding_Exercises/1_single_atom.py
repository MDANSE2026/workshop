import numpy as np


KB = 8.314462618e-3  # AMU nm^2 ps^-2 K^-1 = kJ mol^-1 K^-1
rng = np.random.default_rng()


def next_positions(pos: np.ndarray, vels: np.ndarray, ds: float) -> np.ndarray:
    """Returns the next set of positions after an amount of time ds has
    elapsed.

    Parameters
    ----------
    pos : np.ndarray
        Positions of the atoms in the simulation.
    vels : np.ndarray
        Velocities of the atoms in the simulation.
    ds : float
        The time step.

    Returns
    -------
    np.ndarray
        Next position after an amount of time ds has elapsed.
    """
    pass


def write_positions(pos: np.ndarray, b_len: float, file_name: str):
    """Append positions to xyz file.

    Parameters
    ----------
    pos : np.ndarray
        Positions of the atoms in the simulation.
    b_len : float
        The simulation box length.
    file_name : str
        File append positions to.

    """
    with open(file_name, "a") as f:
        f.write(f"1\n")
        f.write(f'Lattice="{10 * b_len} 0.0 0.0 0.0 {10 * b_len} 0.0 0.0 0.0 {10 * b_len}" pbc="T T T"\n')
        f.write(f"Ar {10 * pos[0]} {10 * pos[1]} {10 * pos[2]}\n")


# numerical simulation of system of a single atom
box_len = 4.585  # nm
positions = box_len * np.array([0.5, 0.5, 0.5])
velocities = np.random.rand(3)  # nm ps^-1
time_step = 0.01  # 0.01 ps = 10 fs
num_time_steps = 100


open("single_atom.xyz", "w").close()  # clear the file
for _ in range(num_time_steps):
    write_positions(positions, box_len, "single_atom.xyz")
    positions = next_positions(positions, velocities, time_step)


# to check your results you can convert this xyz file using the ASE
# converter in MDANSE.
