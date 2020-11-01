# <img src="https://github.com/fdsteffen/fluordynamics/blob/master/docs/source/_static/fluordynamics_logo.png" width="200">
[![Build Status](https://github.com/fdsteffen/fluordynamics/workflows/Fluordynamics%20build/badge.svg)](https://github.com/fdsteffen/fluordynamics/actions)
[![Docs Status](https://github.com/fdsteffen/fluordynamics/workflows/Fluordynamics%20docs/badge.svg)](https://github.com/fdsteffen/fluordynamics/actions)

*Fluordynamics* simplifies the workflow of setting up, running and evaluating Molecular dynamics simulations with extrinsic organic fluorophores for *in silico* FRET experiments.

### Interactive fragment generation

The Python module allows to build new fragments (base, linker and dye) **interactively**. It integrates into PyMOL such that the molecular viewer can be controlled directly from a **Jupyter notebook**. The module also features a dedicated **PyMOL plugin** called *Fluorlabel* to attach the new fragments to a nucleic acid of interest with the click of a button.

### Patching AMBER force fields

Using established pipelines for topology generation such as *Antechamber* and *Acpype*, Fluordynamics builds on top of the AMBERDYES package (Graen et al. *JCTC* 2014) and extends the force field with parameters of common **nucleic acid** linker chemistries.

### *in silico* FRET prediction

MD trajectories with all-atom organic dyes are used to as distance and orientation distributions (:math:`R_\text{DA}` and :math:`\kappa^2`) to compute photon bursts as part of an *in silico* FRET experiment.


### Documentation

*Fluordynamics* is documented [here](https://fdsteffen.github.io/fluordynamics/)


## Download and install

Install Fluordynamics from Github with pip
```
pip install --user git+https://github.com/fdsteffen/fluordynamics.git
```

To generate your own fragments you may further need:

- *PyMOL* https://pymol.org/2/#download
- *Antechamber*  https://ambermd.org/GetAmber.php#ambertools
```
conda install -c conda-forge ambertools=20
```

- *Acpype* https://alanwilter.github.io/acpype/
```
conda install -c conda-forge acpype
```

- A quantum chemistry package such as *Gaussian* https://gaussian.com/ or *GAMESS* https://www.msg.chem.iastate.edu/gamess/


## References

If you use Fluordynamics in your work please refer to the following paper:
 F.D. Steffen, R.K.O. Sigel, R. Börner, *Phys. Chem. Chem. Phys.* **2016**, *18*, 29045-29055. 
 |Steffen2016|

For further information see a list of [related projects](https://fdsteffen.github.io/fluordynamics/references)
