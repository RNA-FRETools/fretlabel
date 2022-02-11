# Installation

Depending on your operating system and preference there are multiple options to install PyMOL and *FRETlabel*
- install PyMOL and add *FRETlabel* as a Plugin manually &rarr; [Install PyMOL and FRETlabel](#install-manually)
- install a Docker image with PyMOL and *FRETlabel* preconfigured (see Docker image) &rarr; [Install with Docker](#install-docker)


## 1.1 Install PyMOL and FRETlabel
<a name="install-manually"></a>
````{tabbed} For Windows
- Get PyMOL from [Schrödinger](https://pymol.org/) or follow the instructions [here](https://pymolwiki.org/index.php/Windows_Install).
- Open the **Anaconda/Miniforge prompt** which comes bundled with PyMOL
````

````{tabbed} For Linux

- Get PyMOL either from [Schrödinger](https://pymol.org/) or from your package manager (e.g. on Ubuntu `apt-get install pymol`). Alternatively, you can compile PyMOL yourself from the source code on [Github](https://github.com/schrodinger/pymol-open-source). 
````
- Install *FRETlabel* by running the following commands

    ```
    git clone https://github.com/RNA-FRETools/fretlabel.git
    cd fretlabel
    pip install fretlabel
    ```

- Locate the installation directory by running

    ```
    fretlabel --path
    ```


## 1.2 Register the Plugin
- Start PyMOL and install the *FRETlabel* GUI with PyMOL's Plugin manager: `Plugin` &rarr; `Plugin manager` &rarr; `Install New Plugin` &rarr; `Choose file...` and select the `fretlabel_gui.py` file located in the directory that was issued by `fretlabel --path`. In the popup window select where you would like to install the plugin (default: `<PyMOL_path>/site-packages/pmg_tk/startup/`). Confirm with OK.


## 1.3 Install with Docker
<a name="install-docker"></a>
As an alternative to the native installation outlined above, you may also use a Docker image with PyMOL and *FRETlabel* preinstalled. Make sure you have [Docker](https://www.docker.com/products/docker-desktop) and an X11 server installed (e.g. [VcXsrv](https://sourceforge.net/projects/vcxsrv/) for Windows, configured with `Wgl="False"` and `DisableAC="True"`). Then pull and run the image from DockerHub (replace `hostdir` with a directory on your host system that you would like to mount into the image)

```
docker pull fdsteffen/fretlabel-pymol:latest
docker run -e DISPLAY=host.docker.internal:0 -v hostdir:/mnt fdsteffen/fretlabel-pymol
```

```{admonition} Incentive or open-source PyMOL
PyMOL was developed by Warren DeLano {cite}`DeLano.2002` and is currently maintained by [Schrödinger](https://pymol.org/). 
Binaries for Windows, Linux and macOS are distributed by Schrödinger under academic and commercial licensing options. The source-code is available on [Github](https://github.com/schrodinger/pymol-open-source).
```

```{tip} 
To generate your own fragments you further need:
- the preprocessing toolbox [*Antechamber*](https://ambermd.org/GetAmber.php#ambertools) packaged with Ambertools
    ```
    conda install -c conda-forge ambertools=20
    ```
- the AnteChamber PYthon Parser interfacE [*Acpype*](https://alanwilter.github.io/acpype/)
    ```
    conda install -c conda-forge acpype
    ```
- A quantum chemistry package such as [*Gaussian*](https://gaussian.com/) or [*GAMESS*](https://www.msg.chem.iastate.edu/gamess/)
```
