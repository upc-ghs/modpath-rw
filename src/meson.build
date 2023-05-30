external_sources = files(
  # gpkde
  'lib' / 'gpkde' / 'Precision.F90',
  'lib' / 'gpkde' / 'Constants.f90',
  'lib' / 'gpkde' / 'Histogram.f90',
  'lib' / 'gpkde' / 'KernelMultiGaussian.f90',
  'lib' / 'gpkde' / 'GridCell.f90', 
  'lib' / 'gpkde' / 'GridProjectedKDE.F90',
  # finterp
  'lib' / 'finterp' / 'linear_interpolation_module.F90',
  # ziggurat random number
  'lib' / 'rngzig' / 'rng_par_zig.f90',
)

sources = files(
  'Coordinate.f90',
  'TrackCellOptions.f90',
  'BudgetListItem.f90',
  'HeadRecordHeader.f90',
  'GridLocation.f90',
  'GlobalData.f90',
  'DispersionData.f90',
  'Solute.f90',
  'Observation.f90',
  'ParticleTrackingOptions.f90',
  'compilerversion.fpp',
  'openspec.f90',
  'Particle.f90',
  'utl7.f',
  'GridLocationList.f90',
  'GeoReference.f90',
  'ModpathUnstructuredBasicData.f90',
  'UtilMisc.f90',
  'ParticleLocation.f90',
  'BudgetRecordHeader.f90',
  'ParticleCoordinate.f90',
  'ParticleGroup.f90',
  'ParticleLocationList.f90',
  'ParticleCoordinateList.f90',
  'utl8.f90',
  'ParticleLocationPtr.f90',
  'ModflowGrid.f90',
  'BudgetReader.f90',
  'TimeDiscretizationData.f90',
  'ModpathSubCellData.f90',
  'ParticlePath.f90',
  'TrackSubCellResult.f90',
  'HeadReader.f90',
  'TrackPathResult.f90',
  'ModflowRectangularGrid.f90',
  'ModpathBasicData.f90',
  'RectangularGridDisMf6.f90',
  'ModpathCellData.f90',
  'TrackSubCell.F90',
  'RectangularGridDisvMf6.f90',
  'StartingLocationReader.f90',
  'RectangularGridDis.f90',
  'RectangularGridDisuMfusg.f90',
  'TrackCellResult.f90',
  'RectangularGridDisuMf6.f90',
  'ModpathCellDataContainer.f90',
  'FlowModelData.f90',
  'ModpathSimulationData.f90',
  'ParticleManager.F90',
  'TrackCell.f90',
  'TransportModelData.f90',
  'ParticleTrackingEngine.f90',
)

extlib = static_library(
           'extlib',
           external_sources
         )
core   = static_library(
           'core',
           sources,  
           dependencies: dependencies, 
           link_with   : [extlib]
         )
theexe = executable(
           'mpathrw', 
           'MPathRW.F90', 
           dependencies: dependencies, 
           link_with   : [extlib,core],
           install     : true
         )
