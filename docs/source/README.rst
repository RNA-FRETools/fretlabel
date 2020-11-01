Fluordynamics - Fluorescent labeling *in silico*
================================================

.. toctree::
   :maxdepth: 2


.. image:: _static/graphical_abstract.png
   :width: 100pc
   :align: center

What is in the box?
*******************

*Fluordynamics* simplifies the workflow of setting up, running and evaluating Molecular dynamics simulations with extrinsic organic fluorophores for *in silico* FRET calculations.

Interactive fragment generation
-------------------------------

The Python module allows to build new fragments (base, linker and dye) **interactively**. It integrates into PyMOL such that the molecular viewer can be controlled directly from a **Jupyter notebook**. The module also features a dedicated **PyMOL plugin** called *Fluorlabel* to attach the new fragments to a nucleic acid of interest with the click of a button. 

Patching AMBER force fields
---------------------------

Together with established pipelines for topology generation such as *Antechamber* and *Acpype*, Fluordynamics patches existing **Gromacs** ports of the AMBER force fields (including AMBERDYES) with new fragment parameters.

*in silico* FRET
----------------

On top of the MD trajectories, photon bursts are computed allowing to simulate a FRET experiment.

------

Download and install
********************

Install fluordynamics from Github with pip ::

    pip install --user git+https://github.com/fdsteffen/fluordynamics.git


To generate your own fragments you may further need:

- *PyMOL* https://pymol.org/2/#download
- *Antechamber*  https://ambermd.org/GetAmber.php#ambertools ::

    conda install -c conda-forge ambertools=20
- *Acpype* https://alanwilter.github.io/acpype/ ::

    conda install -c conda-forge acpype
- A quantum chemistry package such as *Gaussian* https://gaussian.com/ or *GAMESS* https://www.msg.chem.iastate.edu/gamess/

------

Getting started
***************

.. role::  raw-html(raw)
    :format: html

Want to label your nucleic acid *in silico*?
--------------------------------------------
:raw-html:`&rarr;` See the quick intro to the :doc:`PyMOL plugin <fluorlabel>`

Extend the labeling portfolio and build your own linker fragments?
------------------------------------------------------------------
:raw-html:`&rarr;` Visit one of the :doc:`tutorials <fragment_building>`

You have an MD trajectory and want to calculate FRET?
-----------------------------------------------------
:raw-html:`&rarr;` Check out the :doc:`Fluorburst tutorial <fluorburst>`