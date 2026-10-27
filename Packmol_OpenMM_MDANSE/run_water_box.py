import sys

import openmm as mm
import openmm.app as app
import openmm.unit as unit


PACKED_PDB = "water_box.pdb"

TEMPERATURE = 300 * unit.kelvin
PRESSURE = ...
TIMESTEP = ...

NPT_STEPS = int(round(100.0 / TIMESTEP.value_in_unit(unit.picosecond)))
NVT_STEPS = int(round(50.0 / TIMESTEP.value_in_unit(unit.picosecond)))
NVE_EQUIL_STEPS = int(round(25.0 / TIMESTEP.value_in_unit(unit.picosecond)))
NVE_PROD_STEPS = int(round(100.0 / TIMESTEP.value_in_unit(unit.picosecond)))
REPORT_INTERVAL = int(round(1.0 / TIMESTEP.value_in_unit(unit.picosecond)))
DCD_INTERVAl = int(round(20.0 / TIMESTEP.value_in_unit(unit.femtoseconds)))

# load the pdf file and set the box vectors
pdb = app.PDBFile(PACKED_PDB)
topology = pdb.topology
topology.setPeriodicBoxVectors(
    [mm.Vec3(3.105, 0, 0),
     mm.Vec3(0, 3.105, 0),
     mm.Vec3(0, 0, 3.105)] * unit.nanometer
)

# create the force field and add extra particles
forcefield = app.ForceField(...)
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addExtraParticles(forcefield)
topology = modeller.topology

# create the openmm system
system = forcefield.createSystem(
    topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=0.9 * unit.nanometer,
    rigidWater=...,
)

# create and add the barostat
barostat_index = system.addForce(mm.MonteCarloBarostat(PRESSURE, TEMPERATURE, 25))

# create the integrators for our NPT, NVT, and NVT runs
nvt_integrator = ...
nve_integrator = mm.VerletIntegrator(TIMESTEP)
nve_integrator.setConstraintTolerance(1e-6)
integrator = mm.CompoundIntegrator()
integrator.addIntegrator(nvt_integrator)
integrator.addIntegrator(nve_integrator)
integrator.setCurrentIntegrator(0)  # set to the NVT integrator

# create the openmm simulation and set the positions and velocities
simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(modeller.positions)
simulation.context.setVelocitiesToTemperature(TEMPERATURE)

# create some reporters which outputs some details of the simulation
simulation.reporters.append(app.StateDataReporter(
    "equilibration.csv", REPORT_INTERVAL, step=True, time=True,
    potentialEnergy=True, totalEnergy=True, temperature=True, volume=True,
    density=True, separator=","))
simulation.reporters.append(app.StateDataReporter(
    sys.stdout, REPORT_INTERVAL, step=True, temperature=True,
    density=True, speed=True))

# minimisation
print(f"Running OpenMM under the platform: {simulation.context.getPlatform().getName()}")
simulation.minimizeEnergy(maxIterations=5000)

# NPT equilibration
print(f"Running NPT equilibration for: {NPT_STEPS} steps")
simulation.step(NPT_STEPS)

# NVT equilibration
print(f"Running NVT equilibration for: {NVT_STEPS} steps")
system.getForce(barostat_index).setFrequency(0)  # 0 disables the barostat
simulation.context.reinitialize(preserveState=True)  # we modified the force, we have to reinitialize
simulation.step(NVT_STEPS)

# NVE equilibration
print(f"Running NVE equilibration for: {NVE_PROD_STEPS} steps")
integrator.setCurrentIntegrator(1)  # switch to the velocity verlet
simulation.step(NVE_EQUIL_STEPS)

# NVE production
print(f"Running NVE production for: {NVE_PROD_STEPS} steps")
# clear and create new reporters for the production stage
simulation.reporters.clear()
simulation.reporters.append(app.DCDReporter("production.dcd", DCD_INTERVAl))
simulation.reporters.append(app.StateDataReporter(
    "production.csv", REPORT_INTERVAL, step=True, time=True,
    potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
    temperature=True, volume=True, density=True, separator=","))
simulation.reporters.append(app.StateDataReporter(
    sys.stdout, REPORT_INTERVAL, step=True, temperature=True,
    totalEnergy=True, speed=True, totalSteps=NVE_PROD_STEPS,
    progress=True))
# reset the time back to zero and run the production stage
simulation.context.setStepCount(0)
simulation.context.setTime(0 * unit.picosecond)
simulation.step(NVE_PROD_STEPS)

# Save the final structure. This is the topology file needed to the trajectory
# conversion with the DCD converter. We need to use this one and not the
# input PDB file since we might have added some virtual particles.
state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
with open("water_box_final.pdb", "w") as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
