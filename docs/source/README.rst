Fluordynamics - Fluorescent labeling *in silico*
================================================

.. toctree::
   :maxdepth: 2

.. role:: raw-html(raw)
   :format: html

.. image:: _static/graphical_abstract.png
   :width: 100pc
   :align: center

What is in the box?
*******************

*Fluordynamics* simplifies the workflow of setting up, running and evaluating Molecular dynamics simulations with extrinsic organic fluorophores for *in silico* FRET calculations.

Fluorescent labeling *in silico*
--------------------------------

Label your nucleic acid of interest with the click of a button. A **PyMOL plugin** [#]_ called *Fluorlabel* extends the AMBERDYES package (Graen et al. *JCTC* 2014) designed originally for proteins with geometries and force field parameters of common nucleic acid linker chemistries.


Interactive fragment generation
-------------------------------

The module further describes how to build new fragments (base, linker and dye) by integrating with established pipelines for topology generation such as *Antechamber* and *Acpype*.


FRET prediction
---------------

Calculate FRET distributions on the basis of MD simulation with all-atom organic dyes. *Fluorburst* uses distance :math:`R_\text{DA}(t)` and orientation trajectories :math:`\kappa^2(t)` of the fluorophores to compute photon bursts as part of *in silico* FRET experiment.


------


Download and install
********************

- Clone or download the git repository from ::

    git clone https://github.com/fdsteffen/fluordynamics.git

    # or from the mirror
    git clone https://github.com/RNA-FRETools/fluordynamics.git


- Install the PyMOL plugin via the Plugin manager: ``Plugin`` :raw-html:`&rarr;` ``Plugin manager`` :raw-html:`&rarr;` ``Install New Plugin`` :raw-html:`&rarr;` ``Choose file...`` and select the ``__init__`` file from the *fluorlabel* subdirectory.

- To generate your own fragments you may further need:

    - *PyMOL* https://github.com/schrodinger/pymol-open-source
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


-----

References
**********

.. |Steffen2016| image:: https://img.shields.io/badge/DOI-10.1039/C6CP04277E-blue.svg
  :target: https://doi.org/10.1039/C6CP04277E

If you use Fluordynamics in your work please refer to the following paper:
 F.D. Steffen, R.K.O. Sigel, R. Börner, *Phys. Chem. Chem. Phys.* **2016**, *18*, 29045-29055. 
 |Steffen2016|

For further information see a list of :doc:`related projects <references>`

----

.. [#] PyMOL is a trademark of Schrödinger, LLC.