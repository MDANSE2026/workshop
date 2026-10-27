import numpy as np


KB = 8.314462618e-3  # AMU nm^2 ps^-2 K^-1 = kJ mol^-1 K^-1
rng = np.random.default_rng()


def generate_fcc_lattice(n_cells: int):
    """Generates an FCC supercell with a box which is duplicated by
    n_cells for each side.

    Parameters
    ----------
    n_cells : int
        Number of cell to duplicate along each side.

    Returns
    -------
    np.ndarray
        Array of atom positions.
    """
    fcc_pts = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ])
    r = np.arange(n_cells, dtype=float)
    cells = np.stack(np.meshgrid(r, r, r, indexing="ij"), axis=-1).reshape(-1, 3)
    return (cells[:, None, :] + fcc_pts[None, :, :]).reshape(-1, 3)


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


def next_velocities(vels: np.ndarray, acc_0: np.ndarray, acc_1: np.ndarray, ds: float):
    """Returns the next velocities.

    Parameters
    ----------
    vels : np.ndarray
        Velocities of the atoms in the simulation.
    acc_0 : np.ndarray
        Acceleration of the atoms in the simulation at the current time.
    acc_1 : np.ndarray
        Acceleration of the atoms in the simulation at the next time step.
    ds : float
        The time step.

    Returns
    -------
    np.ndarray
        The next velocities.
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


def calculate_acceleration(
        pos: np.ndarray, mass: float, sigma: float, eps: float, b_len: float):
    """Returns the acceleration on each atom. Forces are calculated
    between atom up to a cutoff of one half the box length.

    Parameters
    ----------
    pos : np.ndarray
        Positions of the atoms in the simulation.
    mass : float
        Mass of the atoms in the simulation.
    sigma : float
        Lennard-Jones sigma parameter.
    eps : float
        Lennard-Jones epsilon parameter.
    b_len : float
        The simulation box length.

    Returns
    -------
    np.ndarray
        Acceleration on each atom.
    """
    # diff = pos[:, None, :] - pos[None, :, :]  # big array of separation vectors
    # minimum image convention on diff
    # calculate distance matrix (r2)
    # np.fill_diagonal(r2, np.inf)  # avoids self interaction
    # r2[r2 > (b_len / 2) ** 2] = np.inf  # removes interactions beyond cutoff
    # calculate the forces
    pass


def write_xyz(pos: np.ndarray, vel: np.ndarray, mass, b_len: float, file_name: str):
    """Append positions and velocities to xyz file.

    Parameters
    ----------
    pos : np.ndarray
        Positions of the atoms in the simulation.
    vel : np.ndarray
        Positions of the atoms in the simulation.
    mass : float
        Mass of the atoms in the simulation.
    b_len : float
        The simulation box length.
    file_name : str
        File append positions to.

    """
    pos = 10 * pos
    mom = mass * vel / 100
    n_atms = pos.shape[0]
    with open(file_name, "a") as f:
        f.write(f"{n_atms}\n")
        f.write(f'Lattice="{10 * b_len} 0.0 0.0 0.0 {10 * b_len} 0.0 0.0 0.0 {10 * b_len}" Properties=species:S:1:pos:R:3:momenta:R:3 pbc="T T T"\n')
        for q, p in zip(pos, mom):
            f.write(f"Ar {q[0]} {q[1]} {q[2]}   {p[0]} {p[1]} {p[2]}\n")


# numerical simulation of a simple LJ system
mass = 39.948  # amu
n_cells = 4
cell_len = 0.57311 # nm
box_len = n_cells * cell_len  # nm
temperature = 85  # K
positions = cell_len * generate_fcc_lattice(n_cells)
num_atoms = positions.shape[0]
velocities = random_velocities_from_maxwell_boltzmann(num_atoms, temperature, mass)
time_step = 0.002  # 0.002 ps = 2 fs
num_time_steps = 100
lj_sigma = 0.3405  # nm
lj_eps = 0.99542  # AMU nm^2 ps^-2 = kJ mol^-1


acceleration_0 = calculate_acceleration(positions, mass, lj_sigma, lj_eps, box_len)
open("../lennard_jones.xyz", "w").close()  # clear the file
for i in range(num_time_steps):
    if i % 5 == 0:
        # dump coordinates every 10 fs
        wrapped = wrap_coordinates(positions, box_len)
        write_xyz(wrapped, velocities, mass, box_len, "../lennard_jones.xyz")

    # velocity verlet
    positions = next_positions(positions, velocities, acceleration_0, time_step)
    acceleration_1 = calculate_acceleration(positions, mass, lj_sigma, lj_eps, box_len)
    velocities = next_velocities(velocities, acceleration_0, acceleration_1, time_step)
    acceleration_0 = acceleration_1

