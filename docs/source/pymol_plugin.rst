.. _PyMOL_plugin:

PyMOL Plugin
============

.. toctree::
   :maxdepth: 2

The FluorLabel PyMOL Plugin allows to label an RNA or DNA molecule of interest (provided as a PDB file) with a variety of cyanine, Alexa and Atto dyes. Currently implemented coupling chemistries include:
    
    - **3'-end** via the phosphate or alternatively via a hydrazide integrated into a six-membered sugar ring (see Zhao, *Nucleic Acids Res.* **2018**)
    - **5'-end** via phosphate
    - **uracil** / **deoxythimidine** via C5/C7
    - **(deoxy)adenine** / **(deoxy)cytosine** via an ethenoadduct (see Zhao, *Nucleic Acids Res.* **2018**)


.. image:: _static/FluorLabel.gif
   :width: 400
   :align: center

Installation
------------

1. Clone or download the Fluordynamics repository from Github ::

    git clone https://github.com/fdsteffen/fluordynamics.git


2. Open PyMOL_ and go to :code:`Plugin -> Plugin Manager -> Install New Plugin -> Choose file...`, navigate to the Fluordynamics package and select the :code:`__init__` file in the *fluorlabel* subdirectory

3. Restart PyMOL and open *FluorLabel* from the Plugin menu

.. _PyMOL: https://pymol.org/2/#download
