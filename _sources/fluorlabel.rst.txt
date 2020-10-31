FluorLabel
==========

*FluorLabel* describes a workflow to label nucleic acids *in silico*. For this purpose, a number of labeling chemistries are implemented, including standard 3'/5'-end labeling as well as internal labeling via modified uracils/deoxythymines or ethenoadducts on adenines and cytosines (Zhao, *Nucleic Acids Res.* **2018**). The process of building these base-linker fragments is documented in a series of Jupyter notebooks in the section :ref:`fragment_building`

.. image:: tutorial/images/dye-linker-base.png
   :width: 500
   :align: center

*FluorLabel* already includes a number of cyanines, Alexa and Atto dyes which can be coupled to an RNA or DNA molecule of interest using a :ref:`PyMOL_plugin`. The Plugin integrates abovementioned coupling chemistries in a GUI such that the biomolecules can be labeled with just a few clicks. 

.. image:: _static/FluorLabel_static.jpg
   :width: 500
   :align: center

.. toctree::
   :hidden:

   fragment_building
   pymol_plugin
