## Transport in two-dimensional heterogeneous aquifer 
Simulation considers a rectangular initial condition of solute concentrations released near aquifer inlet. 

Flow model considers an heterogeneous hydraulic conductivity distribution and unit-mean head gradient. The latter is induced with prescribed heads at the aquifer inlet and outlet. 

**Note**: In order to run MODPATH-RW, it is first necessary to execute the MODFLOW-6 simulation. 

### Random walk
Simulation is of type timeseries with forward tracking. An initial distribution of concentrations is specifed with the ``IC`` package. Smoothed reconstruction of concentrations is done only once at the final tracking time.  

```
mpathrw ex07a_mf6.mprw
```

**Note**: Explore the command line interface options for help and some basic usage instructions (``mpathrw -h``).

