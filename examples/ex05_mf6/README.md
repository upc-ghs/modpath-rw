## Transport near a low-permeability zone
The problem is presented in the MT3DMS documentation (Zheng & Wang, 1999), discussing different numerical schemes for advection. Simulation illustrates a solute injection near a low-permeability zone.

Flow model is steady, solved with MODFLOW-6. Concentration at the injection well is stored in the flow model as an auxiliary variable.

**Note**: In order to run MODPATH-RW, it is first necessary to execute the MODFLOW-6 simulation. 

### Random walk: case a
Simulation is of type endpoint, with forward tracking. The program employs the ``SRC`` package to extract the ``CONCENTRATION`` auxiliary variable from the ``WEL-1`` package. With this information, and the user provided particles' mass and release template, the program internally configures a release of particles. 

A sink observation monitors the flux concentration at the extraction well.
 
Simulation is run with the command:

```
mpathrw ex05a_mf6.mprw
```

### Random walk: case a
Similar to the previous case, but now considering that two solutes of different inflow concentration (``CONCENTRATION``,``CONCENTRATION2``) and specific dispersion properties are being injected. Simulation illustrates the usage of the multispecies capabilities of the program. The ``SPC`` package links the dispersion parameters and solutes. The species id's are given in the ``SRC`` package for each auxiliary variable. 

The sink observation monitors flux concentration at the extraction well. Time vector for each solute is different breakthrough data is concatenated vertically, being the first data column the species id. 
 
Simulation is run with the command:

```
mpathrw ex05b_mf6.mprw
```

**Note**: Explore the command line interface options for help and some basic usage instructions (``mpathrw -h``).

## References
Zheng, C. & Wang, P., 1999, MT3DMS: a modular three-dimensional multispecies transport model for simulation of advection, dispersion, and chemical reactions of contaminants in groundwater systems; documentation and user’s guide

[MT3DMS Problem 9 in modflow6-examples@readthedocs](https://modflow6-examples.readthedocs.io/en/master/_examples/ex-gwt-mt3dms-p09.html)
