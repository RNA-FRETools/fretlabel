# Welcome to FRETlabel
*FRETlabel* is a **PyMOL plugin**[^PyMOL] to label nucleic acids with explicit fluorescent dyes. It aims to facilitate the workflow of setting up, running and evaluating **molecular dynamics simulations with explicit organic fluorophores** for *in silico* FRET calculations.

## Features

- **PyMOL plugin for fluorescent labeling**: Label your nucleic acid of interest with the click of a button. The PyMOL plugin extends the AMBERDYES package (Graen et al. *JCTC* 2014) geometries and force field parameters of common nucleic acid linker chemistries.
- **Build new fragments interactively**: Tutorials guide you step-by-step through the process of creating new base, linker and dye fragments by integrating with established pipelines for topology generation such as *Antechamber* and *Acpype*.
- **FRET prediction**: Calculate FRET distributions from MD simulation with all-atom organic dyes. *FRETlabel* integrates with *FRETraj* to compute photon bursts based on distance $R_\text{DA}(t)$ and orientation trajectories $\kappa^2(t)$ of the fluorophores. 

```{figure} images/graphical_abstract.png
---
width: 100%
name: graphical_abstract
align: left
---
Schematic of *FRETlabel*: (i) A fluorescent dye (here Cy3) is coupled to a nucleic acid via a PyMOL plugin. (ii) An existing force field (e.g. AMBERDYES) is patched with parameters for linker fragments to enable MD simulations with explicit fluorophores (dots represent the spatial distribution of the dye).
```

```{admonition} FRETraj
*FRETlabel* uses explicit fluorophores to calculate FRET efficiencies. If you instead like to use an implicit, geometrical dye model that relies on **accessible-contact volumes (ACV)** then have a look our sister project [*FRETraj*](https://rna-fretools.github.io/fretraj/intro.html) (Steffen, *Bioinformatics*, 2021)
``` 

[^PyMOL]: PyMOL is a trademark of Schrödinger, LLC.
