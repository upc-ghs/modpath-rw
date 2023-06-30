## Transport near a low-permeability zone, backward tracking. 
This is the same setup than ``ex05_mf6``, in this case, with backward tracking. 

A high-resolution series of flux concentration at the extraction well is given to the flow model as timeseries data. This series is given as auxiliary variable to the extraction well. 

Flow boundaries remain constant in time, but the simulation considers multiple time steps for the purposes of representing the changes in time of the flux concentration. 

**Note**: In order to run MODPATH-RW, it is first necessary to execute the MODFLOW-6 simulation. 

### Random walk
Simulation is of type endpoint, with backward tracking. The reference time is the end of the MODFLOW simulation. The ``SRC`` package makes use of the ``CONCENTRATION`` variable at the extraction well in order to configure a release of particles in backward direction. Notice that in this case the extraction well can inject particles due to the inversed flow boundaries.

A sink observation monitors the flux concentration at the injection well. In this case, this series should be interpreted as the probability density function of the injection. Time vector in the observation is the MODPATH tracking time. This means that while interpreting the series, it should be considered that the origin is the given reference time and time is reversed. 
 
```
mpathrw ex06a_mf6.mprw
```

**Note**: Explore the command line interface options for help and some basic usage instructions (``mpathrw -h``).

## References
Zheng, C. & Wang, P., 1999, MT3DMS: a modular three-dimensional multispecies transport model for simulation of advection, dispersion, and chemical reactions of contaminants in groundwater systems; documentation and user’s guide

[MT3DMS Problem 9 in modflow6-examples@readthedocs](https://modflow6-examples.readthedocs.io/en/master/_examples/ex-gwt-mt3dms-p09.html)
