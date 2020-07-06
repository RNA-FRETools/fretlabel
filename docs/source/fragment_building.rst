Fragment building
=================

.. toctree::
   :maxdepth: 2
   :hidden:

   tutorial/create_dye_fragment
   tutorial/base_linker_generation
   
.. jupyter-execute::
    
    import numpy as np
    from matplotlib import pyplot
    %matplotlib inline

    x = np.linspace(1E-3, 2 * np.pi)

    pyplot.plot(x, np.sin(x) / x)
    pyplot.plot(x, np.cos(x))
    pyplot.grid()

.. jupyter-execute::
    
    import numpy as np
    from matplotlib import pyplot
    %matplotlib inline

    x = np.linspace(1E-3, 2 * np.pi)

    pyplot.plot(x, np.sin(x) / x)
    pyplot.plot(x, np.cos(x))
    pyplot.grid()




