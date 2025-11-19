from __future__ import print_function
import sys
import os
sys.executable


import numpy
numpy.set_printoptions(threshold=sys.maxsize)
import pandas as pd


## MD related libraries (these should be in a folder called files)
from input.openmm.omm_readinputs import *
from input.openmm.omm_readparams import *
from input.openmm.omm_vfswitch import *
from input.openmm.omm_barostat import *
from input.openmm.omm_restraints import *
from input.openmm.omm_rewrap import *
import MDAnalysis as mda

# Openn MM libraries
from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

# ParmEd libraries
from parmed import load_file
from parmed.openmm.reporters import NetCDFReporter
from parmed import unit as u

# Plumed plugin libraries
Platform.loadPluginsFromDirectory('/usr/local/openmm/lib')
Platform.loadPluginsFromDirectory('/usr/local/openmm/lib/plugins')
from openmmplumed import PlumedForce
print(Platform.getPluginLoadFailures())

# Own libraries
## path to the python script for data generation
Path2script = "/home/mert/Data/Scripting"
sys.path.append(Path2script)
from folding_package.traj_gen import traj_gen
from folding_package.distance_cal import distance_cal
from folding_package.col_ids import col_ids

# ML libraries

from keras import backend as K
from keras.models import Model
from keras.layers import Input, Dense, Lambda
from keras.initializers import RandomUniform, Constant
from keras.optimizers import RMSprop
from keras.constraints import unit_norm
from keras import regularizers
from keras.callbacks import Callback
from keras.losses import mean_squared_error

from folding_package.COLVAR2npy import *
from folding_package.Analyze_prave import *
from folding_package.rave import *
from folding_package.amino import *


num_iter=50

for iteration in range(0, num_iter):
    #Defining input files
        
    psffile='input/md/'+'step3_charmm2omm.psf'
    crdfile='input/md/'+'step3_charmm2omm.crd'
    topparfile='input/md/'+'toppar.str'
    if iteration == 0:
        dcdfile='twister_MgCl2_0.000M_KCl_0.100M_at_303K.200ns.3.dcd'
    else:
        dcdfile='output/md/''output_'+str(iteration-1)+'.dcd'
    # Load parameters
    print("Loading parameters")
    #inputs = read_inputs(inpfile)
    params = read_params(topparfile)
    top = CharmmPsfFile(psffile)
    gro = CharmmCrdFile(crdfile)
    print('Parameters loaded')
    
    #Define universe for MDAnalysis
    rna = mda.Universe(psffile, dcdfile)
    traj_len = len(rna.trajectory)
    print("Length of trajectory is: "+str(traj_len))
    
    ## Atoms selection subroutine
    ura_wc_n3=rna.select_atoms('resname URA and name N3')
    ura_wc_o4=rna.select_atoms('resname URA and name O4')
    ade_wc_n1=rna.select_atoms('resname ADE and name N1')
    ade_wc_n6=rna.select_atoms('resname ADE and name N6')
    gua_wc_o6=rna.select_atoms('resname GUA and name O6')
    gua_wc_n1=rna.select_atoms('resname GUA and name N1')
    gua_wc_n2=rna.select_atoms('resname GUA and name N2')
    cyt_wc_n4=rna.select_atoms('resname CYT and name N4')
    cyt_wc_n3=rna.select_atoms('resname CYT and name N3')
    cyt_wc_o2=rna.select_atoms('resname CYT and name O2')
    
    num_ura_sel = ura_wc_n3.positions.shape[0]
    num_ade_sel = ade_wc_n1.positions.shape[0]
    num_gua_sel = gua_wc_o6.positions.shape[0]
    num_cyt_sel = cyt_wc_n4.positions.shape[0]
    
    print("Number of URACIL: "+str(num_ura_sel))
    print("Number of ADENINE: "+str(num_ade_sel))
    print("Number of GUANINE: "+str(num_gua_sel))
    print("Number of CYTOSINE: "+str(num_cyt_sel))
    
    #Calculating the distances between selected atoms
    dist_n3n1=traj_gen(ura_wc_n3, ade_wc_n1, 1, rna)
    dist_o4n6=traj_gen(ura_wc_o4, ade_wc_n6, 1, rna)
    
    #Merge all order parameters in one dataframe
    OPs=pd.concat([dist_n3n1, dist_o4n6], axis=1, sort=False)
    OPs.to_csv('output/ml/'+'data.csv', sep=' ')
    
    #initialize AMINO package
    if iteration==0:
        
        all_ops = []
        for i in OPs.columns:
            all_ops.append(OrderParameter(i, OPs[i]))
        final_ops = find_ops(all_ops, 10, 100)
        
        print("\nAMINO order parameters:")
        cols = []
        for i in final_ops:
            print(i)
            cols.append(str(i))
    
    OPs_reduced=OPs.filter(cols)
    OPs_reduced.to_csv('input/data_0', sep=' ')
    print("OPs trajectory length: "+str(OPs_reduced.shape[0]))
    with open('input/data_0', "r+") as myfile:
        myfile.seek(0)
        myfile.write("#")
        myfile.close()
    
    
    # Defining inputs for RAVE
    system_name = 'data'
    n_trajs = 1   
    save_path = 'output/rave/'  
    input_path = 'input/'
    T = 300     
   
    
    if iteration == 0: bias = False
    if iteration != 0: bias = True
    
    
    
    #predictive time delay #
    time_delay= list(range(0, 100, 10)) #predictive time delay
    trials = range(4) 
    
    #network variables
    batch_size = 10
    training_size = 100
    #training_size = int(OPs_reduced.shape[0]/batch_size)*batch_size-batch_size
    print("training size: "+str(training_size))
    
    
    op_dim = OPs_reduced.shape[1] 
    rc_dim = 1  
    int_dim = 128  
    s_vari = 0.005 
    learning_rate = 0.0002 
    decay = 0.0       
    epochs = 20 
    
    # End of input definitions for RAVE
    
    rave(system_name, n_trajs, save_path, input_path, T, bias, time_delay, training_size, batch_size, op_dim, rc_dim, int_dim, s_vari, learning_rate, decay, epochs, trials, iteration)
    
    network_info = '_int_dim'+str(int_dim)+'_lr'+str(learning_rate)+'_decay'+str(decay)+'_batch_size'+str(batch_size)
    if not bias:
        system_name2 = 'unbiased_' + system_name
    else:
        system_name2 = system_name
        
    save_result(system_name2, op_dim, time_delay, trials, s_vari, training_size, network_info, save_path)  
    
    # Extracting reactions coordinate weights from RAVE output
    RCs = pd.read_csv('output/rave/final_result_'+system_name2+'_svar'+str(s_vari)+'_train_size'+str(training_size)+network_info+'.txt', delimiter= '\s+')

    # MD part
    if iteration == 0:
        print("EQUILIBRATION STARTS")
        print('Building the system')
        
        inpfile='input/md/'+'step4_equilibration.inp'
        inputs = read_inputs(inpfile)
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
        
        # Defining the metadynamics inputs
        from folding_package.plumed_script import *
        bias_factor=10
        label='RC_1'
        sigma=0.092
        height=10
        pace=100
        plumed_script(OPs_reduced, RCs, iteration, T, bias_factor, label, sigma, height, pace, 'output/meta/')
        
        file = open('output/meta/'+'plumed_s_'+str(iteration)+'.dat',"r+")
        script = file.read()
        file.close()
        
        system.addForce(PlumedForce(script))
        
        integrator = LangevinIntegrator(inputs.temp*kelvin, inputs.fric_coeff/picosecond, inputs.dt*picoseconds)
        
        
        print("Setting platform")
        platform = Platform.getPlatformByName('CUDA')
        #prop = dict(CudaPrecision='single')
        print("Platform has been set")
        
        print("Building simulation context")
        simulation = Simulation(top.topology, system, integrator, platform)
        simulation.context.setPositions(gro.positions)
        
        
        
        
        
        # Energy minimization
        
        if inputs.mini_nstep > 0:
            print("\nEnergy minimization: %s steps" % inputs.mini_nstep)
            simulation.minimizeEnergy(tolerance=inputs.mini_Tol*kilojoule/mole, maxIterations=inputs.mini_nstep)
            print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
        
        
        # Generate initial velocities
        if inputs.gen_vel == 'yes':
            print("\nGenerate initial velocities")
            if inputs.gen_seed:
                simulation.context.setVelocitiesToTemperature(inputs.gen_temp, inputs.gen_seed)
            else:
                simulation.context.setVelocitiesToTemperature(inputs.gen_temp) 
      
        # Production
        if inputs.nstep > 0:
            print("\nMD run: %s steps" % inputs.nstep)
            if inputs.nstdcd > 0:
                odcd = 'output/md/'+'equilibration_'+str(iteration)+'.dcd'
                simulation.reporters.append(DCDReporter(odcd, inputs.nstdcd))
            simulation.reporters.append(
                StateDataReporter(sys.stdout, inputs.nstout, step=True, time=True, potentialEnergy=True, temperature=True, volume=True, progress=True,
                                  remainingTime=True, speed=True, totalSteps=inputs.nstep, separator='\t')
            )
            simulation.reporters.append(PDBReporter('output/md/'+'equilibration_'+str(iteration)+'.pdb', inputs.nstep))
            simulation.step(inputs.nstep) 
        # Write restart file
        orst = 'output/md/'+'equilibration_'+str(iteration)+'.rst'
        
        if 'orst' in locals() or 'orst' in globals():
            state = simulation.context.getState( getPositions=True, getVelocities=True )
            with open(orst, 'w') as f:
                f.write(XmlSerializer.serialize(state))
        if 'ochk' in locals() or 'ochk' in globals():
            with open(ochk, 'wb') as f:
                f.write(simulation.context.createCheckpoint())
                
        print("EQUILIBRATION ENDS")  
    
    
    print("PRODUCTION STARTS")
    print('Building the system')
    inpfile='input/md/'+'step5_production.inp'
    inputs = read_inputs(inpfile)
    
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
    
    # Defining the metadynamics inputs
    from folding_package.plumed_script import *
    bias_factor=10
    label='RC_1'
    sigma=0.092
    height=10
    pace=100
    plumed_script(OPs_reduced, RCs, iteration, T, bias_factor, label, sigma, height, pace, 'output/meta/')
    
    
    file = open('output/meta/'+'plumed_s_'+str(iteration)+'.dat',"r+")
    script = file.read()
    file.close()
    
    system.addForce(PlumedForce(script))

    integrator = LangevinIntegrator(inputs.temp*kelvin, inputs.fric_coeff/picosecond, inputs.dt*picoseconds)
    
    
    print("Setting platform")
    platform = Platform.getPlatformByName('CUDA')
    #prop = dict(CudaPrecision='single')
    print("Platform has been set")
    
    print("Building simulation context")
    simulation = Simulation(top.topology, system, integrator, platform)
    simulation.context.setPositions(gro.positions)
    
    
    if iteration == 0:
        irst = 'output/md/'+'equilibration_'+str(iteration)+'.rst'
    if iteration != 0:
        irst = 'output/md/'+'output_'+str(iteration-1)+'.rst'
    
    if 'irst' in locals() or 'irst' in globals():
        with open(irst, 'r') as f:
            simulation.context.setState(XmlSerializer.deserialize(f.read()))
    if 'ichk' in locals() or 'ichk' in globals():
        with open(ichk, 'rb') as f:
            simulation.context.loadCheckpoint(f.read())
    
    
    # Production
    if inputs.nstep > 0:
        print("\nMD run: %s steps" % inputs.nstep)
        if inputs.nstdcd > 0:
            odcd = 'output/md/'+'output_'+str(iteration)+'.dcd'
            simulation.reporters.append(DCDReporter(odcd, inputs.nstdcd))
        simulation.reporters.append(
            StateDataReporter(sys.stdout, inputs.nstout, step=True, time=True, potentialEnergy=True, temperature=True, volume=True, progress=True,
                              remainingTime=True, speed=True, totalSteps=inputs.nstep, separator='\t')
        )
        simulation.reporters.append(PDBReporter('output/md/'+'output_'+str(iteration)+'.pdb', inputs.nstep))
        simulation.step(inputs.nstep) 
    # Write restart file
    orst = 'output/md/'+'output_'+str(iteration)+'.rst'
    
    if 'orst' in locals() or 'orst' in globals():
        state = simulation.context.getState( getPositions=True, getVelocities=True )
        with open(orst, 'w') as f:
            f.write(XmlSerializer.serialize(state))
    if 'ochk' in locals() or 'ochk' in globals():
        with open(ochk, 'wb') as f:
            f.write(simulation.context.createCheckpoint())
            
    print("PRODUCTION ENDS")  
      
    print("End of iteration: "+str(iteration))
