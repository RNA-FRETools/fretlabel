FluorBurst
==========

.. toctree::
   :hidden:

   tutorial/fluorburst

*FluorBurst* uses a Markov chain to compute fluorescent bursts from an MD trajectory of a biomolecule with fluorophores explicitely included. Distances and orientation factors are calculated from the transition dipoles of the donor and acceptor dye using *gmx dyecouple* which is part of Gromacs. Specifically, the orientation factor :math:`\kappa^2` is calculated as

.. math:: 

   \kappa^2 = (\cos\theta_T - 3\cos\theta_D\cos\theta_A)^2

where :math:`\theta_D` and :math:`\theta_A` are the angles between the vector connecting the central carbons and the transition dipole of the donor or acceptor respectively, and :math:`\theta_T` is the angle between the two transition dipoles. [1]_ 


.. jupyter-execute::
    
   import os
   import fluordynamics as fd
   import matplotlib.pyplot as plt 
   path = os.path.dirname(fd.__file__)

   experiment = fd.fluorburst.Experiment.load(path+'/../docs/source/tutorial/trajectory_examples/dna')

   f, ax = plt.subplots(nrows=1, ncols=1, figsize=(2.5, 2), sharex=False, sharey=False, squeeze=False)
   plt.hist(experiment.FRETefficiencies,bins=20, range=(0,1), color='grey', edgecolor='black')
   ax[0,0].set_xlabel('FRET')
   ax[0,0].set_ylabel('occupancy')


.. [1] J.R. Lakowicz, Principles of fluorescence spectroscopy, 3rd ed., Springer, New York, 2006.
