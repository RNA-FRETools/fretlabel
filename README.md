# <img src="https://github.com/fdsteffen/fluordynamics/blob/master/docs/source/_static/fluordynamics_logo.png" width="50">Fluordynamics
[![Build Status](https://github.com/fdsteffen/fluordynamics/workflows/fluordynamics%20build/badge.svg)](https://github.com/fdsteffen/fluordynamics/actions)

## Introduction

### Interactive fragment generation

Fluordynamics is a package for setting up and running Molecular dynamics simulations with extrinsic organic fluorophores. The eponymous Python module allows to build new fragments (base, linker and dye) interactively. It integrates into PyMOL such that the molecular viewer can be controlled directly from a Jupyter notebook. The module also features a dedicated PyMOL plugin called *Fluorlabel* to attach the new fragments to the nucleic acid of interest. 

### Patching AMBER force fields

Together with established pipelines for topology generation such as *Antechamber* and *Acpype*, Fluordynamics patches existing Gromacs ports of the AMBER force fields (including AMBERDYES, Graen et al. JCTC, 2014) with the newly generated fragment parameters.

### MD simulation setup

The package further includes run files for setting up molecular dynamics simulations in Gromacs.


## Download and install

Install gromacs-tools from Github with pip
```
pip install --user git+https://github.com/fdsteffen/fluordynamics.git
```


## Dependencies

To generate your own fragments make sure you have two following Python dependencies installed:

- numpy
- pandas
- biopandas

You may also need the following programs and packages

- *PyMOL* https://pymol.org/2/#download
- *Antechamber*  https://ambermd.org/GetAmber.php#ambertools
```
conda install ambertools=19 -c ambermd
```

- *Acpype* https://alanwilter.github.io/acpype/
```
conda install -c conda-forge acpype
```

- A quantum chemistry package such as *Gaussian* https://gaussian.com/ or *GAMESS* https://www.msg.chem.iastate.edu/gamess/
