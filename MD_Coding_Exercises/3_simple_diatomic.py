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
    # copy from free particles
    pass


def next_positions(pos: np.ndarray, vels: np.ndarray, acc: np.ndarray, ds: float) -> np.ndarray:
    """Returns the next positions.

    Parameters
    ----------
    pos : np.ndarray
        Positions of the atoms in the simulation.
    vels : np.ndarray
        Velocities of the atoms in the simulation.
    acc : np.ndarray
        Acceleration of the atoms in the simulation.
    ds : float
        The time step.

    Returns
    -------
    np.ndarray
        The next position.
    """
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
    # copy from free particles
    pass


def calculate_acceleration(pos: np.ndarray, mass: float, bond_eq: float, f_k: float):
    """Returns the acceleration on each atom.

    Parameters
    ----------
    pos : np.ndarray
        Positions of the atoms in the simulation.
    mass : float
        Mass of the atoms in the simulation.
    bond_eq : float
        The equilibrium bond length.
    f_k : float
        The force constant between the atoms in the diatomic.

    Returns
    -------
    np.ndarray
        Acceleration on each atom.
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
            f.write(f"N {10 * p[0]} {10 * p[1]} {10 * p[2]}\n")


# numerical simulation of a simple diatomic
mass = 14.007  # amu
box_len = 1  # nm
num_atoms = 2
temperature = 85  # K
force_constant = 1.35e+06  # AMU ps^-2
eq_bond_len = 0.1098  # nm
positions = np.array([
    [0.5, 0.5, 0.5 + 0.95 * eq_bond_len / 2],
    [0.5, 0.5, 0.5 - 0.95 * eq_bond_len / 2]
])
velocities = random_velocities_from_maxwell_boltzmann(2, temperature, mass)
time_step = 0.001  # 0.001 ps = 1 fs
num_time_steps = 10000


open("simple_diatomic.xyz", "w").close()  # clear the file
acceleration = calculate_acceleration(positions, mass, eq_bond_len, force_constant)
positions_0 = positions
write_positions(positions_0, box_len, "simple_diatomic.xyz")
positions_1 = positions + velocities * time_step + 0.5 * acceleration * time_step**2

for i in range(num_time_steps):
    wrapped = wrap_coordinates(positions_1, box_len)
    write_positions(wrapped, box_len, "simple_diatomic.xyz")

    # verlet
    acceleration = calculate_acceleration(positions_1, mass, eq_bond_len, force_constant)
    positions_2 = next_positions(positions_0, positions_1, acceleration, time_step)
    positions_0 = positions_1
    positions_1 = positions_2
