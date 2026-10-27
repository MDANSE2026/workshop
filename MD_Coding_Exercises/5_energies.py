import numpy as np
import matplotlib.pyplot as plt


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
    # copy from simple lennard jones
    pass


def potential_energy(
        pos: np.ndarray, sigma: float, eps: float, b_len: float):
    """Returns the potential energy per atom of the system.

    Parameters
    ----------
    pos : np.ndarray
        Positions of the atoms in the simulation.
    sigma : float
        Lennard-Jones sigma parameter.
    eps : float
        Lennard-Jones epsilon parameter.
    b_len : float
        The simulation box length.

    Returns
    -------
    float
        The potential energy per atom.
    """
    # similar to the force calculation
    pass


def kinetic_energy(vels: np.ndarray, mass: float):
    """Returns the kinetic energy per atom of the system.

    Parameters
    ----------
    vels : np.ndarray
        Velocities of the atoms in the simulation.
    mass : float
        Mass of the atoms in the simulation.

    Returns
    -------
    float
        The kinetic energy per atom.
    """
    pass


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
num_time_steps = 5000
lj_sigma = 0.3405  # nm
lj_eps = 0.99542  # AMU nm^2 ps^-2 = kJ mol^-1


pes = []
kes = []
tot_es = []
times = []

acceleration_0 = calculate_acceleration(positions, mass, lj_sigma, lj_eps, box_len)
open("../lennard_jones.xyz", "w").close()  # clear the file
for i in range(num_time_steps):

    if i % 5 == 0:
        pe = potential_energy(positions, lj_sigma, lj_eps, box_len)
        ke = kinetic_energy(velocities, mass)
        pes.append(pe)
        kes.append(ke)
        tot_es.append(pe + ke)
        times.append(i * time_step)

    # velocity verlet
    positions = next_positions(positions, velocities, acceleration_0, time_step)
    acceleration_1 = calculate_acceleration(positions, mass, lj_sigma, lj_eps, box_len)
    velocities = next_velocities(velocities, acceleration_0, acceleration_1, time_step)
    acceleration_0 = acceleration_1


plt.plot(times, pes, label="Potential energy")
plt.plot(times, kes, label="Kinetic energy")
plt.plot(times, tot_es, label="Total energy")

plt.xlabel("Time (ps)")
plt.ylabel("Energy (kJ/mol)")
plt.legend()
plt.grid(True)
plt.show()
