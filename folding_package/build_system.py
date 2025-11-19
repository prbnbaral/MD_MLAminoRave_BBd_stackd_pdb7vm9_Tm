## MD related libraries (these should be in a folder called files)
from files.omm_readinputs import *
from files.omm_readparams import *
from files.omm_vfswitch import *
from files.omm_barostat import *
from files.omm_restraints import *
from files.omm_rewrap import *
import MDAnalysis as mda

# Openn MM libraries
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

# ParmEd libraries
from parmed import load_file
from parmed.openmm.reporters import NetCDFReporter
from parmed import unit as u


def build_system(top, gro, params, inputs, system):
    print('Building the system')
    gen_box(top, gro)
    system = top.createSystem(params, nonbondedMethod=inputs.coulomb,
                              nonbondedCutoff=inputs.r_off*nanometers,
                              constraints=inputs.cons,
                              ewaldErrorTolerance=inputs.ewald_Tol,
                             )

    if inputs.vdw == 'Force-switch': system = vfswitch(system, top, inputs)
    if inputs.pcouple == 'yes':      system = barostat(system, inputs)
    if inputs.rest == 'yes':         system = restraints(system, gro, inputs)
    print('System has been built')
    
    return system
    

def build_context(top, system, integrator, platform):
    print("Building simulation context")
    simulation = Simulation(top.topology, system, integrator, platform)
    simulation.context.setPositions(gro.positions)
    
    if 'irst' in locals() or 'irst' in globals():
        with open(irst, 'r') as f:
            simulation.context.setState(XmlSerializer.deserialize(f.read()))
    if 'ichk' in locals() or 'ichk' in globals():
        with open(ichk, 'rb') as f:
            simulation.context.loadCheckpoint(f.read())

    # Re-wrap
    if 'rewrap' in locals() or 'orst' in globals():
        print("\nRe-wrapping the simulation")
        simulation = rewrap(simulation)
    return simulation
