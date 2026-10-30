import numpy as np


KB = 8.314462618e-3  # AMU nm^2 ps^-2 K^-1 = kJ mol^-1 K^-1
rng = np.random.default_rng()


def random_velocities_from_maxwell_boltzmann(n_atms: int, temp: float, mass: float) -> np.ndarray:
    """Generate velocities whose speeds are drawn from a maxwell boltzmann
    distribution.

    Parameters
    ----------
    n_atms : int
        Number of atoms to generate velocities for.
    temp : float
        Temperature in kelvin.
    mass : float
        Mass of the particles.

    Returns
    -------
    np.ndarray
        Numpy array of particle velocities drawn from maxwell boltzmann
        distribution.
    """
    pass


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
    # copy from single atom
    pass


def wrap_coordinates(pos: np.ndarray, b_len: float) -> np.ndarray:
    """Folds the atom positions back into the simulation box.

    Parameters
    ----------
    pos : np.ndarray
        Positions of the atoms in the simulation.
    b_len : float
        The simulation box length.

    Returns
    -------
    np.ndarray
        Atom positions folded into the simulation box.
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
    n_atms = pos.shape[0]
    with open(file_name, "a") as f:
        f.write(f"{n_atms}\n")
        f.write(f'Lattice="{10 * b_len} 0.0 0.0 0.0 {10 * b_len} 0.0 0.0 0.0 {10 * b_len}" pbc="T T T"\n')
        for p in pos:
            f.write(f"Ar {10 * p[0]} {10 * p[1]} {10 * p[2]}\n")


# numerical simulation of system of free particles
mass = 39.948  # amu
box_len = 4.585  # nm
num_atoms = 2048
temperature = 85  # K
positions = box_len * np.random.rand(num_atoms, 3)
velocities = random_velocities_from_maxwell_boltzmann(num_atoms, temperature, mass)
time_step = 0.1  # 0.1 ps = 100 fs
num_time_steps = 1000


open("free_particles.xyz", "w").close()  # clear the file
for _ in range(num_time_steps):
    write_positions(positions, box_len, "free_particles.xyz")
    positions = next_positions(positions, velocities, time_step)
    positions = wrap_coordinates(positions, box_len)


# to check your results you can convert this xyz file using the ASE
# converter in MDANSE. Then try running the temperature calculation
# it should be near the temperature specified above
