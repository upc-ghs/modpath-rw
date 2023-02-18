module GlobalDataModule
    implicit none
    integer, save :: narealsp, issflg, nper
    integer, save :: mpbasUnit, disUnit, tdisUnit, gridMetaUnit, headUnit,      & 
      headuUnit, budgetUnit, traceModeUnit, binPathlineUnit
    integer, save :: inUnit, pathlineUnit, endpointUnit, timeseriesUnit,        &
      mplistUnit, mpsimUnit, traceUnit, budchkUnit, aobsUnit, logUnit
    integer, save :: gpkdeUnit, obsUnit, dspUnit, rwoptsUnit   ! RWPT 
    integer, save :: spcUnit, icUnit, srcUnit, impUnit         ! RWPT 
    integer, save :: particleGroupCount
    integer, save :: gridFileType
    integer, save :: logType
    integer, parameter :: niunit = 100
    character(len=200) :: mpnamFile, mpsimFile, mplistFile, mpbasFile, disFile, &
      tdisFile, gridFile, headFile, budgetFile, traceFile, gridMetaFile,        &
      mplogFile, gpkdeFile, obsFile, dspFile, rwoptsFile, spcFile, icFile,      & 
      srcFile, impFile ! RWPT
end module
