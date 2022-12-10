# modpath-rw
A Random Walk Particle Tracking Code for Solute Transport in Heterogeneous Aquifers

## Overview
MODPATH-RW is an extension of MODPATH version 7 implementing particles displacement based on the Random Walk Particle Tracking (RWPT) method. Code is built on top of MODPATH-OMP, so it can displace particles in parallel with the OpenMP library. 

## Compilation
Download source code and access the `make` folder. Here execute one of the makefiles with

```
make -f makefile-gfortran-pc
```

This will compile source files and will output a binary, by default with name `mpathrwpt.exe`. 
Compilation has been verified with ``gfortran@>=9.2.1`` and with ``ifort`` from Intel ``oneAPI@2021.3.0``.


## Running parallel simulations
``modpath`` interface has been extended to simplify execution of parallel runs via the command line. 
Some alternatives are:

- Use the ``-np`` argument to specify the number of processes:
```
mpath7omp example.mpsim -np 4 
```

- Use the ``-parallel`` argument to run in parallel employing the maximum number of available processors:

```
mpath7omp example.mpsim -parallel 
```

- Specify the ``OMP_NUM_THREADS`` environment variable:
```
OMP_NUM_THREADS=4 mpath7omp example.mpsim
```
In this last case, if the variable is defined at a system level, it will be employed when submitting a normal run instruction ``mpath7omp example.mpsim``


## Parallel output for timeseries
Three different output protocols for timeseries running in parallel have been implemented. The protocol can be selected via the command line argument ``-tsoutput``:

- ``-tsoutput 1``: is the default format, output is performed into a single output unit with OpenMP thread exclusive clause (critical). Only difference versus a serial run is that the output file contains non-sorted particle indexes.
- ``-tsoutput 2``: timeseries records are written into thread specific binary units and then consolidated into a single file after each timeseries output time. Timeseries file generated with this format does not contains a file header.
- ``-tsoutput 3``: timeseries records are written into thread specific output units. Timeseries file header is only written to output unit related to the first thread ``1_example.timeseries``. Initial particle positions are also written to the file of the first thread.

## Contributing
[Contribution guidelines for this project](readme/CONTRIBUTING.md)

## License
MIT License

## Resources

* [MODPATH](https://www.usgs.gov/software/modpath-particle-tracking-model-modflow)
* [Central modpath-v7 repository](https://github.com/MODFLOW-USGS/modpath-v7)
* [MODPATH-OMP](https://github.com/MARSoluT/modpath-omp)
* [gfortran](https://gcc.gnu.org/wiki/GFortran)
* [OpenMP](https://www.openmp.org/)
* [Intel oneApi HPC toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html)
* [MIT License](https://mit-license.org/)
