# modpath-rw
A Random Walk Particle Tracking Code for Solute Transport in Heterogeneous Aquifers

## Overview
MODPATH-RW is an extension of MODPATH-v7 implementing particles displacement based on the Random Walk Particle Tracking (RWPT) method. Code is built on top of MODPATH-OMP, so it can displace particles in parallel with the OpenMP library. 

## Get the code 
MODPATH-RW requires the module for smoothed reconstruction of concentrations [GPKDE](https://github.com/upc-ghs/gpkde). While cloning, do it recursively with the command

```
git clone --recursive https://github.com/upc-ghs/modpath-rw.git
```

Previous instruction should clone both the modified files of MODPATH-RW plus the ``gpkde`` submodule. The latter is required for compilation of the program. For users cloning via ``ssh`` it is convenient to configure ``git`` with the following

```
git config --global url."ssh://git@".insteadOf https://
```

If you decide to clone not-recursively, then be sure to clone the [GPKDE](https://github.com/upc-ghs/gpkde) repository into the root folder of MODPATH-RW as makefiles assume this submodule exists.


## Compilation
Download source code and access the `make` folder. Here execute one of the makefiles with

```
make -f makefile-gfortran-pc
```

This will compile source files and will output a binary, by default with name `mpathrwpt`. Compilation has been verified with ``gfortran@>=9.2.1`` and with ``ifort`` from Intel ``oneAPI@2021.3.0``.


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


## Input files
Users familiarized with [FloPy](https://github.com/modflowpy/flopy) are encouraged to write input files for MODPATH-RW with the extension [flopyrw](https://github.com/modflowpy/flopyrw). The latter is based on the classes for MODPATH-v7, adapted for specific MODPATH-RW requirements.


## Contributing
Follow the [contribution guidelines](readme/CONTRIBUTING.md) for this project.

## License
MIT License

## Resources

* [MODPATH](https://www.usgs.gov/software/modpath-particle-tracking-model-modflow)
* [MODPATH-v7](https://github.com/MODFLOW-USGS/modpath-v7)
* [MODPATH-OMP](https://github.com/upc-ghs/modpath-omp)
* [GPKDE](https://github.com/upc-ghs/gpkde)
* [flopy](https://github.com/modflowpy/flopy)
* [flopyrw](https://github.com/upc-ghs/flopyrw)
* [gfortran](https://gcc.gnu.org/wiki/GFortran)
* [Intel oneApi HPC toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html)
* [OpenMP](https://www.openmp.org/)
* [MIT License](https://mit-license.org/)
