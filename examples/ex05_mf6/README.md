## Transport near a low-permeability zone
This problem is presented in the MT3DMS documentation (Zheng & Wang, 1999), discussing different numerical schemes for advection. Problem illustrates solute injection near a low-permeability zone. Flow model is steady, solved with MODFLOW-6 and concentration at the injection well is stored in the flow model as an auxiliary variable.

**Note**: In order to run MODPATH-RW, it is first necessary to execute the MODFLOW-6 simulation. 

### Random walk 
Simulation is of type endpoint and forward tracking. The program employs the ``SRC`` package to extract concentration and from here defines the release of particles. An observation monitors the flux concentration at the extraction well. 
 
MODPATH-RW name file includes the following packages:

```
MPBAS      ex05a_mf6.mpbas
GRBDIS     ex05_mf6.dis.grb
TDIS       ex05_mf6.tdis
HEAD       ex05_mf6.hds
BUDGET     ex05_mf6.bud
DSP        ex05a_mf6.dsp
RWOPTS     ex05a_mf6.rwopts
SRC        ex05a_mf6.src
OBS        ex05a_mf6.obs
```

Simulation is run with the command:

```
mpathrw ex05a_mf6.mprw
```

Explore the program command line interface options for help and some basic instructions (``mpathrw -h``).

## References
Zheng, C. & Wang, P., 1999, MT3DMS: a modular three-dimensional multispecies transport model for simulation of advection, dispersion, and chemical reactions of contaminants in groundwater systems; documentation and user’s guide

[MT3DMS Problem 9 in modflow6-examples@readthedocs](https://modflow6-examples.readthedocs.io/en/master/_examples/ex-gwt-mt3dms-p09.html)
