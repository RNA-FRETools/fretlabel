Parameter file
==============

*Fluorburst* uses a parameter file to define the settings of the Markov-chain Monte Carlo simulation.
The settings follow in large parts the ones introduced in Hoefling et al., *Comput. Phys. Commun.* (2013). The structure is as follows: ::

    {
    "dyes":{
        "tauD":0.75,
        "tauA":1.5,
        "QD":0.2,
        "QA":0.3
        },
    "sampling":{
        "nbursts":2000,
        "skipframesatstart":0,
        "skipframesatend":1000,
        "multiprocessing":true
        },
    "fret":{
        "R0":5.4,
        "kappasquare":0.666666,
        "no_gamma":false,
        "quenching_radius":1.0
        },
    "species":{
        "name":["all"],
        "unix_pattern_rkappa": ["*.dat"],
        "unix_pattern_don_coords": ["*Cy3*.xvg"],
        "unix_pattern_acc_coords": ["*Cy3*.xvg"],
        "probability": [1]
        },
    "bursts":{
        "lower_limit":20,
        "upper_limit":150,
        "lambda":-2.3,
        "QY_correction":false,
        "averaging": "all"
        }
    }

- **dyes** defines the fluorescence lifetimes (:code:`tauD` & :code:`tauA`) and quantum yields (:code:`QD` & :code:`QA`) of the donor and acceptor dyes
- **sampling** specifies the total number of bursts to calculate (:code:`nbursts`) as well as the number of frames to skip at the beginning and end of each MD trajectory (:code:`skipframesatstart` & :code:`skipframesatend`) and whether to perform the calculation on multiple cores in parallel
- **fret** defines the Förster radius (:code:`R0`), the dye distance below which contact quenching will occur (:code:`quenching_radius`)