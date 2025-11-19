# A template to run Machine Learning combined MD  

## Overview

This repository provides a **template workflow** for combining **Machine Learning** with **enhanced sampling Molecular Dynamics (Metadynamics)** simulations. The workflow utilizes machine learning to predict appropriate reaction coordinates for metadynamics simulations, with reaction coordinates spanned by combination of center of mass of atoms of BACKBONE DISTANCE and STACKING DISTANCE (different nucleic bases) of a DNA hairpin (PDB ID: 7vm9).

The computational framework integrates OpenMM, metadynamics, and deep learning for iterative simulations, originally developed for studying the twister ribozyme with Mg²⁺ ions using both additive and Drude polarizable force fields. 
For example, refer to 
*J. Chem. Phys.* **161**, 225102 (2024)  
DOI: [10.1063/5.0241246](https://doi.org/10.1063/5.0241246)

## Directory Structure
```
simulation_folder/
├── input/
│   ├── ml/
│   ├── md/
│   │   ├── inputs (psf, crd, inp...)
│   │   ├── toppar_drude.str
│   │   ├── drude_toppar_2020/
│   │   ├── toppar.str
│   │   └── toppar/
│   └── openmm/
│       ├── omm_vfswitch.py
│       ├── omm_restraints.py
│       ├── omm_barostat.py
│       ├── omm_readinputs.py
│       └── omm_readparams.py
├── output/
│   ├── meta/
│   ├── ml/
│   ├── md/
│   └── rave/
├── restraints/
│   ├── prot_ext.txt
│   └── prot_pos.txt
├── toppar/
├── drude_2024/
├── folding_package/
└── folding_drude_na_mixedops.py 
```

## Directory Descriptions

### Input Files

**`toppar/`** - Contains all additive force field related files including `.rtf`, `.str`, and `.prm` files.

**`drude_2024/`** - Contains all Drude force field related `.str` files.

**`input/`** - Contains subfolders with input files organized by workflow component:
- **`openmm/`** - Files for running unbiased molecular dynamics simulations
- **`md/`** - Force field files (`toppar.str` for additive, `toppar_drude.str` for Drude), equilibration/production settings, and initial simulation configurations
- **`ml/`** - Initially empty; populated by the code with machine learning input files during execution

### Output Files

**`output/`** - Contains subfolders for output files (should be created empty):
- **`meta/`** - Metadynamics PLUMED input files generated during each iteration
- **`ml/`** - Machine learning related output files
- **`md/`** - Unbiased MD simulation results
- **`rave/`** - RAVE code output files

### Additional Components

**`restraints/`** - Contains restraint definition files:
- `prot_pos.txt` - Point restraint definitions
- `prot_ext.txt` - External restraining forces for secondary structure

**`folding_package/`** - Python package containing the complete iterative workflow code (imported by main script)

**`folding_drude_na_mixedops.py`** - Main execution script that coordinates deep learning, OpenMM, and metadynamics components

## Required Setup and Configuration

Before running the workflow, modify the following settings in `main_code_drude.py`:

### Path Configuration
```python
Path2script = "/home/username/python_codes"
```
Change this to the path where you've placed the complete workflow code.

### OpenMM Plugin Path
```python
Platform.loadPluginsFromDirectory('/home/username/modules/local/openmm/lib/plugins')
```
Update this path to point to your OpenMM PLUMED plugin location.

### Additional Settings
Customize metadynamics and deep learning parameters in the corresponding sections of `main_code_drude.py`.

## Dependencies

The following Python packages were used:
```
keras==2.8.0
MDAnalysis==2.2.0
numpy==1.23.1
pandas==1.4.2
parmed==3.4.3
OpenMM==7.7.0
```

## Usage

To run the simulation workflow:
```bash
python folding_drude_na_mixedops.py -i 0 1 2 3 4
```

This command executes:
- Iteration #0: Equilibration
- Iterations 1-4: Production runs

## Workflow Overview

This simulation workflow iteratively combines:
1. **OpenMM** for molecular dynamics simulations
2. **Metadynamics** for enhanced sampling
3. **Deep Learning** for analysis and optimization
4. **RAVE** for additional analysis

The system is designed for protein folding studies using both additive and Drude polarizable force fields.
