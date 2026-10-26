"""Adapted from ASE tutorials, generates an isotropic trajectory."""
from asap3 import EMT
from ase import Atoms, units
from ase.io.trajectory import Trajectory
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.nptberendsen import NPTBerendsen
from ase.md.verlet import VelocityVerlet

size = 6

T = 900  # Kelvin

# Set up a crystal
small_atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          symbol="Cu",
                          pbc=True)
small_atoms.set_chemical_symbols(["Cu", "Au", "Cu", "Au"])

atoms = small_atoms*[size, size, size]

# Describe the interatomic interactions with the Effective Medium Theory
atoms.calc = EMT()

# Running MD with constant temperature and pressure
# this will allow the unit cell volume to change
dyn = NPTBerendsen(atoms, timestep=5.0 * units.fs, temperature_K=T,
                   taut=100 * units.fs, pressure_au=1.01325 * units.bar,
                   taup=1000 * units.fs, compressibility_au=4.57e-7 / units.bar)

def printenergy(a:Atoms=atoms):  # store a reference to atoms in the definition.
    """Function to print the potential, kinetic and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    volume = a.get_volume()
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK V=%.3f)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), volume, epot + ekin))


dyn.attach(printenergy, interval=50)

print("NPT run is about to start")
printenergy()
dyn.run(10000)

# Next step is an NVE run
dyn2 = VelocityVerlet(atoms, timestep=5.0*units.fs)
traj = Trajectory('isotropic_CuAu_T900K.traj', 'w', atoms)
dyn2.attach(traj.write, interval=1)
dyn2.attach(printenergy, interval=10)

print("NVE run is about to start")
dyn2.run(500)
print("Finished! All done!")
