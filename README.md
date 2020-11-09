# <img src="https://github.com/fdsteffen/fluordynamics/blob/master/docs/source/_static/fluordynamics_logo.png" width="200">
[![Build Status](https://github.com/fdsteffen/fluordynamics/workflows/Fluordynamics%20build/badge.svg)](https://github.com/fdsteffen/fluordynamics/actions)
[![Docs Status](https://github.com/fdsteffen/fluordynamics/workflows/Fluordynamics%20docs/badge.svg)](https://github.com/fdsteffen/fluordynamics/actions)

*Fluordynamics* simplifies the workflow of setting up, running and evaluating Molecular dynamics simulations with extrinsic organic fluorophores for *in silico* FRET experiments. Documentation including tutorials can be found [here](https://fdsteffen.github.io/fluordynamics/).


<img src="docs/source/_static/graphical_abstract.png" width="550">

### Fluorescent labeling *in silico*

Label your nucleic acid of interest with the click of a button. A **PyMOL plugin**<sup>[1](#pymol)</sup> called *Fluorlabel* extends the AMBERDYES package (Graen et al. *JCTC* 2014) designed originally for proteins with geometries and force field parameters of common nucleic acid linker chemistries.

### Interactive fragment generation

The module further describes how to build new fragments (base, linker and dye) by integrating with established pipelines for topology generation such as *Antechamber* and *Acpype*.

### FRET prediction

Calculate FRET distributions on the basis of MD simulation with all-atom organic dyes. *Fluorburst* uses distance and orientation trajectories of the fluorophores to compute photon bursts as part of an *in silico* FRET experiment.


## Download and install

Clone or download the git repository from

```
git clone https://github.com/fdsteffen/fluordynamics.git
```

or from its mirror
```
git clone https://github.com/RNA-FRETools/fluordynamics.git
```

- Install the PyMOL plugin via the Plugin manager: `Plugin` &rarr; `Plugin manager` &rarr; `Install New Plugin` &rarr; `Choose file...` and select the *Fluorlabel* `__init__` file.

- To generate your own fragments you may further need:

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

F.D. Steffen, R.K.O. Sigel, R. Börner, *Phys. Chem. Chem. Phys.* **2016**, *18*, 29045-29055. [![](https://img.shields.io/badge/DOI-10.1039/C6CP04277E-blue.svg)](https://doi.org/10.1039/C6CP04277E)


For further information see a list of [related projects](https://fdsteffen.github.io/fluordynamics/references).

----

<sup><a name="pymol">1</a></sup> PyMOL is a trademark of Schrödinger, LLC.