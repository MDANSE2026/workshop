import itertools as it
import sys

import openmm as mm
from openmm import app, unit


OUTPUT_DCD_NAME = "argon.dcd"
OUTPUT_PDB_NAME = "argon.pdb"

TEMPERATURE = 85.0 * unit.kelvin
FRICTION = 1.0 / unit.picosecond
TIMESTEP = 2.0 * unit.femtosecond

NVT_STEPS = 12000
NVE_EQUIL_STEPS = 12000
NVE_PROD_STEPS = 60060
REPORT_INTERVAL = 60

MASS = 39.948 * unit.amu

N_CELLS = 4
CELL_LEN = 0.57311  # nm
BOX_LEN = N_CELLS * CELL_LEN  # nm
LJ_SIGMA = 3.405 * unit.angstrom
LJ_EPS = 0.99542 * unit.kilojoules_per_mole


# create topology and set positions
fcc_basis = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5)]
topology = app.Topology()
chain = topology.addChain()
argon = app.Element.getBySymbol("Ar")

positions = []
for (i, j, k), (bx, by, bz) in it.product(
    it.product(range(N_CELLS), repeat=3), fcc_basis
):
    residue = topology.addResidue("AR", chain)
    topology.addAtom("Ar", argon, residue)
    positions.append(mm.Vec3(i + bx, j + by, k + bz) * CELL_LEN)
positions = unit.Quantity(positions, unit.nanometer)

topology.setPeriodicBoxVectors(
    (
        mm.Vec3(BOX_LEN, 0, 0),
        mm.Vec3(0, BOX_LEN, 0),
        mm.Vec3(0, 0, BOX_LEN)
    ) * unit.nanometer
)

# create openmm system, and create and add forces
system = mm.System()
system.setDefaultPeriodicBoxVectors(*topology.getPeriodicBoxVectors())

force = mm.NonbondedForce()
force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
force.setCutoffDistance(1.0 * unit.nanometer)
force.setUseSwitchingFunction(True)
force.setSwitchingDistance(0.9 * unit.nanometer)
force.setUseDispersionCorrection(True)

for _ in range(topology.getNumAtoms()):
    system.addParticle(MASS)
    force.addParticle(0.0, LJ_SIGMA, LJ_EPS)

system.addForce(force)

# create the langevin and verlet integrators
nvt_integrator = mm.LangevinMiddleIntegrator(TEMPERATURE, FRICTION, TIMESTEP)
nve_integrator = mm.VerletIntegrator(TIMESTEP)
integrator = mm.CompoundIntegrator()
integrator.addIntegrator(nvt_integrator)
integrator.addIntegrator(nve_integrator)
integrator.setCurrentIntegrator(0)

# create the openmm simulation and set the positions and velocities
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(TEMPERATURE)

# create some reporters which outputs some details of the simulation
simulation.reporters.append(app.StateDataReporter(
    sys.stdout, REPORT_INTERVAL, step=True, temperature=True,
    density=True, speed=True))

# minimisation
print(f"Running OpenMM under the platform: {simulation.context.getPlatform().getName()}")
simulation.minimizeEnergy()

# NVT equilibration
print(f"Running NVT equilibration for: {NVT_STEPS} steps")
simulation.step(NVT_STEPS)
simulation.reporters.clear()

# NVE equilibration
print(f"Running NVE equilibration for: {NVE_EQUIL_STEPS} steps")
integrator.setCurrentIntegrator(1)
simulation.step(NVE_EQUIL_STEPS)

# NVE production
print(f"Running NVE production for: {NVE_PROD_STEPS} steps")
# clear and create new reporters for the production stage
simulation.reporters.clear()
simulation.reporters.append(app.DCDReporter(OUTPUT_DCD_NAME, REPORT_INTERVAL))
simulation.reporters.append(app.StateDataReporter(
    sys.stdout, REPORT_INTERVAL, step=True, temperature=True,
    totalEnergy=True, speed=True, totalSteps=NVE_PROD_STEPS,
    progress=True))
# reset the time back to zero and run the production stage
simulation.context.setStepCount(0)
simulation.context.setTime(0 * unit.picosecond)
simulation.step(NVE_PROD_STEPS)

# Save the final structure. This is the topology file needed to the trajectory
# conversion with the DCD converter.
state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
with open(OUTPUT_PDB_NAME, "w") as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)