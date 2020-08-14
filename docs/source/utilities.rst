Utilities
=========

gromacs-tools
-------------

*gromacs-tools* is a set of bash scripts to set up and run molecular dynamics simulation with Gromacs. It currently includes the following scripts:

- *resp_fit.sh* implements two-stage restrained electrostatic potential (RESP) fitting using Antechamber_
- *solvate.sh* generates a solvent box around the molecule of interest, neutralizes the charge and carries out an energy minimization
- *single_run.sh* performs initial temperature and pressure equilibration and starts a production MD run
- *continue_run.sh* restarts an interrupted or terminated MD run from the last checkpoint file


.. _Antechamber: https://ambermd.org/GetAmber.php#ambertools


.. toctree::
   :hidden:
