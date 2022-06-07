# modpath-omp
MODPATH version 7 with parallel processing of particles implemented with the OpenMP library. 

## Compilation
Download the code and access the `make` folder. Here execute the `makefile` with

```
make -f makefile-gfortran-pc
```
This will compile source files and will output a binary, by default with name `mpath7omp.exe`. 
Compilation has been verified with ``gfortran@9.2.1`` and with ``ifort`` from Intel oneApi ``2021.3.0``.


## Running parallel simulations
``modpath`` interface has been extended to simplify execution of parallel runs via the command line. 
Some alternatives are:

- Use the ``np`` argument to specify the number of processes:
```
mpath7omp -np 4 example.mpsim
```

- Use the ``parallel`` argument to run in parallel employing the maximum number of available processors

```
mpath7omp -parallel example.mpsim
```

- Specify the ``OMP_NUM_THREADS`` environment variable.
```
OMP_NUM_THREADS=4 mpath7omp example.mpsim
```
In this last case, if the variable is defined at a system level, it will be employed when submitting a normal run instruction ``mpath7omp example.mpsim``


## Parallel output for timeseries
Three different output protocols for timeseries running in parallel have been implemented. The protocol can be selected via the command line argument ``-tsoutput``:

- ``-tsoutput 1``: is the default format, output is performed in a single output unit with OpenMP thread exclusive clause (critical)
- ``-tsoutput 2``: timeseries records are written into thread specific binary units and then consolidated in a single file after each timeseries output time.
- ``-tsoutput 3``: timeseries records are written into thread specific output units.


## Resources

* [MODPATH](https://www.usgs.gov/software/modpath-particle-tracking-model-modflow)
* [Central modpath-v7 repository](https://github.com/MODFLOW-USGS/modpath-v7)
* [gfortran](https://gcc.gnu.org/wiki/GFortran)
* [OpenMP](https://www.openmp.org/)
* [Intel oneApi HPC toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html)
