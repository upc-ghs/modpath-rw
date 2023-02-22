# Solute transport near a low-permeability zone
This problem is presented in the MT3DMS documentation (Zheng & Wang, 1999), discussing different numerical schemes for advection. Problem illustrates solute injection near a low-permeability zone.

Flow model is steady (with MODFLOW-6) and MODPATH-RW model injects particles with the ``SRC`` package, using auxiliary variables given to the flow model. 

MODPATH-RW name file includes the following packages:

```
MPBAS      p09sim.mpbas
GRBDIS     gwf-p09-mf6.dis.grb
TDIS       gwf-p09-mf6.tdis
HEAD       gwf-p09-mf6.hds
BUDGET     gwf-p09-mf6.bud
DSP        p09sim.dsp
RWOPTS     p09sim.rwopts
OBS        p09sim.obs
SRC        p09sim.src
```

In order to run the program, it is first necessary to execute the MODFLOW-6 simulation. MODPATH-RW simulation file is configured to run a MODPATH-RW timeseries and can be executed with the command

```
mpathrw p09sim.mprw
```

The example considers observation cells at injection and extraction wells. Three different files are provided as alternative to illustrate the different cell specification methods. The user would need to modify the ``OBS`` file name in order to load a configuration. The alternatives are:

- ``p09sim.obs``: resident observation cells are specified as list of cells.

- ``p09sim.obs2``: resident observation cells are specified as an internal array.

- ``p09sim.obs3``: resident observation cells are specified from an external file (``obs1cellarea.csv``)


Explore the program command line interface options for help and some basic instructions (``mpathrw -h``).

## References
Zheng, C. & Wang, P., 1999, MT3DMS: a modular three-dimensional multispecies transport model for simulation of advection, dispersion, and chemical reactions of contaminants in groundwater systems; documentation and user’s guide

[MT3DMS Problem 9 in modflow6-examples@readthedocs](https://modflow6-examples.readthedocs.io/en/master/_examples/ex-gwt-mt3dms-p09.html)
