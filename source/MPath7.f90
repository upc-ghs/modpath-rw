!  MPath7.f90 
!  PROGRAM: MPath7

!  This software is preliminary or provisional and is subject to revision.       ! kluge provisional
!  It is being provided to meet the need for timely best science. The software
!  has not received final approval by the U.S. Geological Survey (USGS).
!  No warranty, expressed or implied, is made by the USGS or the U.S. Government
!  as to the functionality of the software and related material nor shall the
!  fact of release constitute any such warranty. The software is provided on
!  the condition that neither the USGS nor the U.S. Government shall be held
!  liable for any damages resulting from the authorized or unauthorized use of
!  the software.

    program MPath7
!*********************************************************************************
! Main program code for USGS MODPATH particle tracking model - Version 7
!
!   Specifications:
!---------------------------------------------------------------------------------
    use GlobalDataModule,only : niunit, narealsp, issflg, nper, mpbasUnit,      &
        disUnit, tdisUnit, gridMetaUnit, headUnit, headuUnit, budgetUnit,       &
        inUnit, pathlineUnit, endpointUnit, timeseriesUnit, binPathlineUnit,    &
        mplistUnit, traceUnit, budchkUnit, aobsUnit, logUnit, mpsimUnit,        &
        dispersionUnit,                                                         & ! RWPT
        traceModeUnit, mpnamFile, mplistFile, mpbasFile, disFile, tdisFile,     &
        gridFile, headFile, budgetFile, mpsimFile, traceFile,  gridMetaFile,    &
        mplogFile, logType, particleGroupCount, gridFileType
    use UtilMiscModule,only : ulog
    use utl8module,only : freeunitnumber, ustop, ugetnode ! GPDKE
    use ModpathCellDataModule,only : ModpathCellDataType
    use ModpathBasicDataModule,only : ModpathBasicDataType
    use ModpathSimulationDataModule,only : ModpathSimulationDataType
    use BudgetReaderModule,only : BudgetReaderType
    use HeadReaderModule,only : HeadReaderType
    
    use ModflowRectangularGridModule,only : ModflowRectangularGridType
    use RectangularGridDisModule,only : RectangularGridDisType
    use RectangularGridDisMf6Module,only : RectangularGridDisMf6Type
    use RectangularGridDisvMf6Module,only : RectangularGridDisvMf6Type
    use RectangularGridDisuMfusgModule,only : RectangularGridDisuMfusgType    
    
    use TimeDiscretizationDataModule,only : TimeDiscretizationDataType
    use ParticleTrackingEngineModule,only : ParticleTrackingEngineType
    use FlowModelDataModule, only: FlowModelDataType
    use TransportModelDataModule, only: TransportModelDataType
    use TrackPathResultModule,only : TrackPathResultType
    use ParticleLocationModule,only : ParticleLocationType
    use ParticleCoordinateModule,only : ParticleCoordinateType
    use ParticleGroupModule,only : ParticleGroupType
    use ParticleModule,only : ParticleType
    use ParticleManagerModule
    use BudgetRecordHeaderModule,only : BudgetRecordHeaderType
    use GeoReferenceModule,only : GeoReferenceType
    use CompilerVersion,only : get_compiler
    use SoluteModule, only : SoluteType ! RWPT
    use ObservationModule, only : ObservationType ! OBS
    use GridProjectedKDEModule, only : GridProjectedKDEType ! GPKDE
    use omp_lib ! OpenMP
    !--------------------------------------------------------------------------
    implicit none
    
    ! Variables declarations
    type(HeadReaderType),allocatable :: headReader
    type(BudgetReaderType), allocatable :: budgetReader
    
    class(ModflowRectangularGridType), pointer :: modelGrid
    class(RectangularGridDisType), allocatable, target :: disGrid
    class(RectangularGridDisMf6Type), allocatable, target :: disMf6Grid
    class(RectangularGridDisvMf6Type), allocatable, target :: disvMf6Grid
    class(RectangularGridDisuMfusgType), allocatable, target :: disuMfusgGrid

    type(TimeDiscretizationDataType), allocatable :: tdisData
    type(ParticleTrackingEngineType), allocatable,target :: trackingEngine
    type(FlowModelDataType), allocatable :: flowModelData
    type(TransportModelDataType), allocatable, target:: transportModelData  ! RWPT
    type(ModpathBasicDataType), allocatable, target :: basicData
    type(ModpathSimulationDataType), allocatable, target :: simulationData
    type(ModpathCellDataType), allocatable, target :: cellData
    type(TrackPathResultType), target :: trackPathResult
    type(ParticleLocationType) :: pLoc
    type(ParticleCoordinateType),pointer :: pCoordFirst, pCoordLast, pCoordTP
    type(ParticleCoordinateType) :: pCoord
    type(ParticleGroupType),pointer :: pGroup
    type(ParticleType),pointer :: p
    type(BudgetRecordHeaderType) :: budgetRecordHeader
    type(GeoReferenceType) :: geoRef
    type(GridProjectedKDEType), allocatable:: gpkde                 ! GPKDE
    type(ModpathCellDataType) :: cellDataBuffer                     ! RWPT
    type(SoluteType), pointer :: solute                             ! RWPT
    doubleprecision,dimension(:),allocatable :: timePoints
    doubleprecision,dimension(:),allocatable :: tPoint
    integer,dimension(7) :: budgetIntervalBins
    doubleprecision,dimension(6) :: budgetIntervalBreaks
    logical :: traceModeOn, unitOpened
    integer :: budgetIntervalBreakCount, maxErrorCell
    doubleprecision :: maxError
    integer :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
    doubleprecision :: elapsedTime
    integer :: groupIndex, particleIndex, pendingCount,  &
      activeCount, timePointCount, tPointCount, pathlineRecordCount
    integer :: stressPeriodCount, recordHeaderCount
    integer :: n, m, ktime, kfirst, klast, kincr, period, step, nt, count,      &
      plCount, tsCount, status, itend, particleID, topActiveCellNumber,         &
      auxCount
    integer :: bufferSize, cellConnectionCount
    integer,dimension(:),allocatable :: buffer
    doubleprecision :: t, stoptime, maxTime, tsMax, time
    character(len=132) message
    character(len=20) version
    character(len=75) terminationMessage
    character(len=80) compilerVersionText
    logical :: isTimeSeriesPoint, timeseriesRecordWritten

    ! GPKDE
    doubleprecision, dimension(:,:), allocatable :: activeParticleCoordinates
    doubleprecision, dimension(:), allocatable   :: activeParticleMasses
    integer :: activeCounter, itcount, ns, npg, pgid
    doubleprecision, dimension(:,:), allocatable :: gpkdeDataCarrier
    doubleprecision, dimension(:), allocatable :: gpkdeWeightsCarrier

    ! OBSERVATIONS
    integer :: nlines, io, irow, krow, nobs, nit, countTS, nTimesHigher, cellNumber 
    integer :: baserow, lastrow, srow, timeIndex, solCount
    doubleprecision :: dTObsSeries
    doubleprecision :: initialTime, initialGlobalX, initialGlobalY, initialGlobalZ, QSinkCell 
    doubleprecision, dimension(3) :: sbuffer
    type( ObservationType ), pointer :: obs => null()
    doubleprecision, allocatable, dimension(:,:) :: obsSinkFlowInTime ! (ntimes,ncells)
    doubleprecision, allocatable, dimension(:)   :: obsAccumSinkFlowInTime ! (ntimes)
    doubleprecision, allocatable, dimension(:)   :: qSinkBuffer
    integer            :: idColFormat
    character(len=200) :: colFormat
    character(len=200) :: qSinkFormat
    doubleprecision    :: obsAccumPorousVolume, dX, dY, dZ, porosity
    doubleprecision    :: modelX, modelY
    integer            :: sequenceNumber, timePointIndex, timeStep
    doubleprecision, allocatable, dimension(:,:) :: BTCPerSolute
    doubleprecision, allocatable, dimension(:,:) :: BTCHistPerSolute
    logical :: anyFromThisSolute = .false.

    ! Parallel variables
    integer :: ompNumThreads
    integer :: ompThreadId
    integer :: baseTimeseriesUnit
    integer :: timeseriesBinUnit
    integer :: reclen, lastRecord
    integer, allocatable, dimension(:) :: timeseriesTempUnits
    integer, allocatable, dimension(:) :: timeseriesRecordCounts
    character(len=200), allocatable, dimension(:) :: timeseriesTempFiles
    character(len=200) :: tempChar
    logical :: parallel = .false.
    integer :: tsOutputType = 0

    ! Interface for timeseries output protocol
    procedure(TimeseriesWriter), pointer :: WriteTimeseries=>null()
!---------------------------------------------------------------------------------
    
    ! Set version
    version = '0.0.1'
    
    call get_compiler(compilerVersionText)
    write(*,'(1x/a,a)') 'MODPATH-RW Version ', version
    write(*,'(a)') compilerVersionText
    write(*,*)
    
    ! Set the default termination message
    terminationMessage = "Normal termination."
    
    ! Assign dedicated file unit numbers
    disUnit = 101
    endpointUnit = 102
    pathlineUnit = 103
    timeseriesUnit = 104
    mplistUnit = 105
    traceUnit = 106
    budchkUnit = 107
    aobsUnit = 108
    logUnit = 109
    mpsimUnit = 110
    tdisUnit = 111
    mpbasUnit = 112
    headUnit = 113
    budgetUnit = 114
    traceModeUnit = 115
    binPathlineUnit = 116
    gridMetaUnit = 117
    dispersionUnit = 118 ! RWPT
    baseTimeseriesUnit = 6600 ! OpenMP

    ! Parse the command line for simulation file name, log file name, and options
    call ParseCommandLine(mpsimFile, mplogFile, logType, parallel, tsOutputType)
    ! Open the log file (unless -nolog option)
    if (logType /= 0) then
        open(unit=logUnit, file=mplogFile, status='replace', form='formatted', access='sequential')
    else
        logUnit = -logUnit
    end if
    ! Get the number of threads for the parallel region
    if ( parallel ) then
        ompNumThreads = omp_get_max_threads()
    else
        ompNumThreads = 1
    end if
    call ulog('Command line parsed.', logUnit)

    ! If simulation file name not on command line, prompt user for name
    if (mpsimFile == "") then
        call ulog('Prompt for the name of the MODPATH simulation file.', logUnit)
        call PromptSimulationFile(mpsimFile)
    end if
    
    ! Read the first two records of the simulation file to get the names of the 
    ! name file (mpnamFile) and the listing file (mplistFile)
    allocate(simulationData)
    open(unit=mpsimUnit, file=mpsimFile, status='old', form='formatted', access='sequential')
    call ulog('Read the first two records of the simulation file to get mpnamFile and mplistFile ...', logUnit)
    call simulationData%ReadFileHeaders(mpsimUnit)
    mpnamFile = simulationData%NameFile
    mplistFile = simulationData%ListingFile
    
    ! Open the MODPATH output listing file
    open(unit=mplistUnit, file=mplistFile, status='replace', form='formatted', access='sequential')
    
    write(mplistUnit,'(1x/a,a)') 'MODPATH Version ', version
    write(mplistUnit,'(a)') compilerVersionText
    write(mplistUnit, *)
    write(mplistUnit, '(a)') 'This software has been approved for release by the U.S. Geological'
    write(mplistUnit, '(a)') 'Survey (USGS). Although the software has been subjected to rigorous'    
    write(mplistUnit, '(a)') 'review, the USGS reserves the right to update the software as needed'    
    write(mplistUnit, '(a)') 'pursuant to further analysis and review. No warranty, expressed or'    
    write(mplistUnit, '(a)') 'implied, is made by the USGS or the U.S. Government as to the'    
    write(mplistUnit, '(a)') 'functionality of the software and related material nor shall the'    
    write(mplistUnit, '(a)') 'fact of release constitute any such warranty. Furthermore, the'    
    write(mplistUnit, '(a)') 'software is released on condition that neither the USGS nor the U.S.'    
    write(mplistUnit, '(a)') 'Government shall be held liable for any damages resulting from its'    
    write(mplistUnit, '(a)') 'authorized or unauthorized use. Also refer to the USGS Water'    
    write(mplistUnit, '(a)') 'Resources Software User Rights Notice for complete use, copyright,'    
    write(mplistUnit, '(a)') 'and distribution information.'    
    write(mplistUnit, *)
   


    ! Read the MODPATH name file
    call ReadNameFile(mpnamFile, mplistUnit, gridFileType)
    
    ! Process spatial and time discretization data
    call ulog('Allocate rectangular unstructured grid component.', logUnit)
    allocate(modelGrid)
    call ulog('Allocate time discretization data component ...', logUnit)
    allocate(tdisData)
    
    write(mplistUnit, '(1x/a)') 'Grid data'
    write(mplistUnit, '(a)')    '---------'
      
    select case (gridFileType)
        case (1) 
            ! MODFLOW-2005 discretization file (DIS)
            ! Read spatial and time discretization. 
            write(mplistUnit, '(a,1x)') 'Grid file type: MODFLOW-2005 discretization file (DIS)'
            call ulog('Allocate disGrid.', logUnit)
            allocate(disGrid)
            call ulog('Read grid file.', logUnit)
            call disGrid%ReadData(disUnit, gridMetaUnit, mplistUnit, stressPeriodCount)
            modelGrid => disGrid
            call tdisData%ReadData(disUnit, mplistUnit, stressPeriodCount)
            
            ! Close discretization files
            close(disUnit)
            
            ! RWPT
            ! isUnstructured remains as false

        case (2)
            ! MODPATH spatial(MPUGRID) and time (TDIS) discretization files 
            ! Read spatial discretization
            write(mplistUnit, '(a,1x)') 'Grid file type: MODFLOW-USG unstructured grid file (DISU).'
            call ulog('Allocate disuMfusgGrid.', logUnit)
            allocate(disuMfusgGrid)
            call ulog('Read grid file.', logUnit)
            call disuMfusgGrid%ReadData(disUnit, gridMetaUnit, mplistUnit, stressPeriodCount)
            modelGrid => disuMfusgGrid
            call tdisData%ReadData(disUnit, mplistUnit, stressPeriodCount)
            ! Close discretization file
            close(disUnit)
        
            ! RWPT: Is there a way to check if is really Unstructured ? 
            ! In the meantime, assume it is at RectangularGridDisuMfusgInit1 

        case (3)
            ! MODFLOW-6 DIS binary grid file
            ! Read spatial discretization
            write(mplistUnit, '(a,1x)') 'Grid file type: MODFLOW-6 DIS binary grid file.'
            call ulog('Allocate disMf6Grid.', logUnit)
            allocate(disMf6Grid)
            call ulog('Read DIS binary grid file.', logUnit)
            call disMf6Grid%ReadData(disUnit, gridMetaUnit, mplistUnit)
            modelGrid => disMf6Grid
            
            ! Read time discretization file
            if(len_trim(tdisFile) .gt. 0) then
                write(mplistUnit, '(a,1x)') 'Time discretization file type: MODFLOW-6 time discretization file.'
                call ulog('R ead time discretization data component ...', logUnit)
                call tdisData%ReadData(tdisUnit, mplistUnit)
            else
                call ulog('The time discretization file was not specified.', logUnit)
                call ustop('The time discretization file was not specified.')
            end if
            
            ! Close discretization files
            close(disUnit)
            close(tdisUnit)
           
            ! RWPT
            ! isUnstructured remains as false

            ! Notice that this grid could have different (distributed)
            ! delr, delc. Meaning that although is structured,
            ! is not regular.

        case (4)
            ! MODFLOW-6 DISV binary grid file
            ! Read spatial discretization
            write(mplistUnit, '(a,1x)') 'Grid file type: MODFLOW-6 DISV binary grid file.'
            call ulog('Allocate disvMf6Grid.', logUnit)
            allocate(disvMf6Grid)
            call ulog('Read DISV binary grid file.', logUnit)
            call disvMf6Grid%ReadData(disUnit, gridMetaUnit, mplistUnit)
            modelGrid => disvMf6Grid
           
            ! Read time discretization file
            if(len_trim(tdisFile) .gt. 0) then
                write(mplistUnit, '(a,1x)') 'Time discretization file type: MODFLOW-6 time discretization file.'
                call ulog('Read time discretization data component ...', logUnit)
                call tdisData%ReadData(tdisUnit, mplistUnit)
            else
                call ulog('The time discretization file was not specified.', logUnit)
                call ustop('The time discretization file was not specified.')
            end if
            
            ! Close discretization files
            close(disUnit)
            close(tdisUnit)
            
            ! RWPT: In RectangularGridDisvMf6Module it is checked whether it is unstructured or not,
            ! while counting midPoints, CheckRectangular, checking if isSmoothed 

        case (5)
            ! MODFLOW-6 DISU binary grid file
            write(mplistUnit, '(1x,a)') 'MODFLOW-6 DISU binary grid files are not yet supported. Stop.' 
            stop
            
        case default
            write(mplistUnit, '(1x,a)') 'Unknown grid file type. Stop.'
            stop
            
        end select
    
    ! Write connection data
    if (logType == 2) then
        call ulog('Skip output of cell connection data.', logUnit)
    else if (logType == 1) then
        write(logUnit, *)
        write(logUnit, '(1x,a)') '----------------------------------------------------------'
        write(logUnit, '(1x,a)') 'Cell connection data:'
        write(logUnit, '(1x,a)') '----------------------------------------------------------'
        write(logUnit, '(1x,a)') 'Format has two lines for each cell, listed by cell number.'
        write(logUnit, '(1x,a)') 'Line 1: Cell connections'
        write(logUnit, '(1x,a)') 'Line 1: Face assignment codes'
        write(logUnit, '(1x,a)') '----------------------------------------------------------'
        bufferSize = 25
        allocate(buffer(bufferSize))
        do n = 1, modelGrid%CellCount
            call modelGrid%GetJaCellConnections(n, buffer, bufferSize, cellConnectionCount)
            if (cellConnectionCount == 0) cycle
            write(logUnit, '(1x,25i8)') (buffer(m), m = 1, cellConnectionCount)
            call modelGrid%GetCellConnectionFaces(n, buffer, bufferSize, cellConnectionCount)
            write(logUnit, '(1x,25i8)') (buffer(m), m = 1, cellConnectionCount)
            write(logUnit, *)        
        end do
    end if
    
    ! Initialize the georeference data
    call geoRef%SetData(modelGrid%OriginX, modelGrid%OriginY, modelGrid%RotationAngle)
    
    ! Initialize the budgetReader component
    call ulog('Allocate budget reader component.', logUnit)
    allocate(budgetReader)
    call ulog('Open budget file in budget reader.', logUnit)
    write(mplistUnit, *)
    call budgetReader%OpenBudgetFile(budgetFile, budgetUnit, mplistUnit)
    if(budgetReader%GetFileOpenStatus()) then
        write(mplistUnit, '(1x,a)') 'The budget file was opened successfully.'
    else
        call ustop('An error occurred processing the budget file. Stopping.')
    end if
    
    ! Initialize the headReader component
    call ulog('Allocate head reader component.', logUnit)
    allocate(headReader)
    call ulog('Open head file in head reader component.', logUnit)
    call headReader%OpenFile(headFile, headUnit, mplistUnit)
    
    ! Read the MODPATH basic data file
    call ulog('Allocate MODPATH basic data component.', logUnit)
    allocate(basicData)
    call ulog('Read MODPATH basic data component.', logUnit)   
    call basicData%ReadData(mpbasUnit, mplistUnit, modelGrid)

    ! Read the remainder of the MODPATH simulation file
    call ulog('Read the remainder of the MODPATH simulation data component.', logUnit)
    call simulationData%ReadData(mpsimUnit, mplistUnit, basicData%IBound, tdisData, modelGrid)
       
    ! Budget File Data Summary
    ! If budget output option = 2, then write a list of budget record headers.
    call WriteBudgetFileInfo(mplistUnit, budgetReader) 
    if(simulationData%BudgetOutputOption .eq. 2) call WriteBudgetRecordHeaders(mplistUnit, budgetReader)

        
    ! Initialize the particle tracking engine:
    call ulog('Allocate particle tracking engine component.', logUnit)
    allocate(trackingEngine)
    allocate(flowModelData)
    call flowModelData%Initialize(headReader, budgetReader, modelGrid,&
                                    basicData%HNoFlow, basicData%HDry )
    call flowModelData%SetIBound(basicData%IBound,modelGrid%CellCount)
    call flowModelData%SetPorosity(basicData%Porosity, modelGrid%CellCount)
    call flowModelData%SetZones(simulationData%Zones, modelGrid%CellCount)
    call flowModelData%SetRetardation(simulationData%Retardation, modelGrid%CellCount)
    call flowModelData%SetDefaultIface(basicData%DefaultIfaceLabels, &
            basicData%DefaultIfaceValues, basicData%DefaultIfaceCount)
    call ulog('Initialize particle tracking engine component.', logUnit)
    if ( simulationData%TrackingOptions%RandomWalkParticleTracking ) then 
        ! Initialize transportModelData
        allocate( transportModelData ) 
        call transportModelData%Initialize( modelGrid )
        call transportModelData%ReadData( dispersionUnit, simulationData%DispersionFile, mplistUnit, &
                        simulationData, flowModelData, basicData%IBound, modelGrid, simulationData%TrackingOptions )
        call trackingEngine%Initialize(modelGrid, simulationData%TrackingOptions, flowModelData, transportModelData)
    else 
        call trackingEngine%Initialize(modelGrid, simulationData%TrackingOptions, flowModelData)
    end if 
    ! The trackingEngine initialization is complete

    ! Prepare to stop if there are no particles to track
    if(simulationData%TotalParticleCount .eq. 0) then
        terminationMessage = 'The simulation was terminated because there are no particles to track.'
        goto 100
    end if

    ! Initialize GPKDE reconstruction 
    if ( simulationData%TrackingOptions%GPKDEReconstruction ) then
        allocate( gpkde )
        ! Initialization should be performed once grid properties are known.
        ! Moreover, reconstruction can employ a grid different than flow model grid.
        ! So for USG grids, reconstructed information could be obtained in a regular
        ! rectangular grid, given particles position.
        call ulog('Initialize GPKDE object and output unit', logUnit)
        call gpkde%Initialize(& 
            simulationData%TrackingOptions%gpkdeDomainSize,                          &
            simulationData%TrackingOptions%gpkdeBinSize,                             &
            domainOrigin=simulationData%TrackingOptions%gpkdeDomainOrigin,           &
            nOptimizationLoops=simulationData%TrackingOptions%gpkdeNOptLoops,        &
            databaseOptimization=simulationData%TrackingOptions%gpkdeKernelDatabase, &
            minHOverLambda=simulationData%TrackingOptions%gpkdeKDBParams(1),         &
            deltaHOverLambda=simulationData%TrackingOptions%gpkdeKDBParams(2),       &
            maxHOverLambda=simulationData%TrackingOptions%gpkdeKDBParams(3)          &
        )
        ! Initialize output unit/file
        open(unit=simulationData%TrackingOptions%gpkdeOutputUnit, &
             file=simulationData%TrackingOptions%gpkdeOutputFile, &
           status='replace', form='formatted', access='sequential')
        
        !! Is there a way to verify whether gpkde and the flow-model have the
        !! same cells structure ?
        !select case (gridFileType)
        !    case (1)
        !      ! MODFLOW-2005 discretization file (DIS)
        !      if(&
        !        ( gpkde%nBins(1) .eq. modelGrid%columnCount ) .and. &
        !        ( gpkde%nBins(2) .eq. modelGrid%rowCount    ) .and. &
        !        ( gpkde%nBins(3) .eq. modelGrid%layerCount  ) )
        !        ! Is the same grid
        !        print *, 'MF62005DIS: YES IS THE SAME !'
        !        continue
        !      end if
        !    case (2)
        !      ! MODPATH spatial(MPUGRID) and time (TDIS) discretization files 
        !      continue
        !    case (3) 
        !      ! MODFLOW-6 DIS binary grid file
        !      if(&
        !        ( gpkde%nBins(1) .eq. modelGrid%columnCount ) .and. &
        !        ( gpkde%nBins(2) .eq. modelGrid%rowCount    ) .and. &
        !        ( gpkde%nBins(3) .eq. modelGrid%layerCount  ) )
        !        ! Is the same grid
        !        print *, 'MF6DIS: YES IS THE SAME !'
        !        continue
        !      end if 
        !    case (4)
        !      ! MODFLOW-6 DISV binary grid file
        !      continue
        !    case (5)
        !      ! MODFLOW-6 DISU binary grid file
        !      continue
        !end select

    end if


    ! Compute range of time steps to use in the time step loop
    message ='Compute range of time steps. Prepare for time step loop'
    call ulog(message, logUnit)
    !
    kfirst = tdisData%FindContainingTimeStep(simulationData%ReferenceTime)
    if(simulationData%TrackingDirection .eq. 1) then
        klast = tdisData%CumulativeTimeStepCount
        kincr = 1
    else
        klast = 1
        kincr = -1
    end if 


    ! Set the appropriate value of stoptime. Start by setting stoptime to correspond to the start or the
    ! end of the simulation (depending on the tracking direction)
    if(simulationData%TrackingDirection .eq. 1) then
        stoptime = tdisData%TotalTimes(tdisData%CumulativeTimeStepCount) -      &
          simulationData%ReferenceTime
        call tdisData%GetPeriodAndStep(tdisData%CumulativeTimeStepCount, period, step)
        call flowModelData%LoadTimeStep(period, step)
        if ( simulationData%TrackingOptions%RandomWalkParticleTracking ) then 
            call transportModelData%LoadTimeStep(period, step)
        end if
    else
        stoptime = simulationData%ReferenceTime
        call tdisData%GetPeriodAndStep(tdisData%CumulativeTimeStepCount, period, step)
        call flowModelData%LoadTimeStep(1, 1)
        if ( simulationData%TrackingOptions%RandomWalkParticleTracking ) then 
            call transportModelData%LoadTimeStep(1, 1)
        end if
    end if
    !
    if(simulationData%StoppingTimeOption .eq. 2) then
        ! Set stoptime to 1.0d+30 if the EXTEND option is on and the boundary time step is steady state.
        ! If the boundary time step is transient, leave stoptime set to correspond to the beginning or 
        ! end of the simulation.
        if(flowModelData%SteadyState) stoptime = 1.0d+30
    else if(simulationData%StoppingTimeOption .eq. 3) then
        ! If a specific stoptime was specified, always apply it if there is a steady-state time step at the beginning
        ! or end of the time domain of the simulation.
        if(flowModelData%SteadyState) then
            stoptime = simulationData%StopTime
        else
        ! If the boundary time step is transient, do not set stoptime to the specified value if it would extend beyond
        ! the time domain of the simulation.
            if(simulationData%StopTime .lt. stoptime)                           &
              stoptime = simulationData%StopTime
        end if
    end if
    
    write(mplistUnit, '(1x/a,e15.7)')                                           &
      'The simulation will be run with stoptime = ', stoptime

    write(*,*)
    write(*,'(A)') 'Run particle tracking simulation ...'    
    write(mplistUnit, *)
    write(mplistUnit, *)
    write(mplistUnit,'(1X,A)') 'Run particle tracking simulation ...'
    
    ! Allocate tPoint array
    tPointCount = 0
    if(simulationData%SimulationType .eq. 2) then
        tPointCount = simulationData%TimePointCount
        if(tPointCount .gt. 0) tPointCount = 1
    ! RWPT
    else if( (simulationData%SimulationType .ge. 3) .and. &
             (simulationData%SimulationType .lt. 7) ) then
        tPointCount = 1
    end if
    if(allocated(tPoint)) deallocate(tPoint)
    allocate(tPoint(tPointCount))
    
    ! Open particle output files

    ! Endpoint
    open(unit=endpointUnit, file=simulationData%EndpointFile, status='replace', &
      form='formatted', access='sequential')

    ! Pathline
    if((simulationData%SimulationType .eq. 2) .or.                              &
       (simulationData%SimulationType .eq. 4) .or.                              & 
       (simulationData%SimulationType .eq. 6)) then
        open(unit=pathlineUnit, file=simulationData%PathlineFile,               &
          status='replace', form='formatted', access='sequential')
        open(unit=binPathlineUnit, status='scratch', form='unformatted',        &
          access='stream', action='readwrite')
        call WritePathlineHeader(pathlineUnit, simulationData%TrackingDirection,&
          simulationData%ReferenceTime, modelGrid%OriginX, modelGrid%OriginY,   &
          modelGrid%RotationAngle)
    end if

    ! Timeseries
    if((simulationData%SimulationType .eq. 3) .or.                              &
       (simulationData%SimulationType .eq. 4) .or.                              &
       (simulationData%SimulationType .eq. 5) .or.                              & ! RWPT
       (simulationData%SimulationType .eq. 6)) then                               ! RWPT

        ! Allocate arrays for parallel output
        ! In serial will be allocated with 
        ! ompNumThreads = 1, but is not used for writing outputs,
        ! it should be present only for consistency with TimeseriesWriter interface
        allocate( timeseriesTempUnits(ompNumThreads) )
        allocate( timeseriesRecordCounts(ompNumThreads) )

        ! If not parallel 
        if ( .not. parallel ) then 
            ! Default output 
            open(unit=timeseriesUnit, file=simulationData%TimeseriesFile,           &
              status='replace', form='formatted', access='sequential')
            call WriteTimeseriesHeader(timeseriesUnit,                              &
              simulationData%TrackingDirection, simulationData%ReferenceTime,       &
              modelGrid%OriginX, modelGrid%OriginY, modelGrid%RotationAngle)
            ! Assign writing interface
            WriteTimeseries => WriteTimeseriesRecordSerial
        else

            ! Select parallel output protocol
            select case( tsOutputType ) 
                case (1)
                    ! Default critical output 
                    open(unit=timeseriesUnit, file=simulationData%TimeseriesFile,           &
                      status='replace', form='formatted', access='sequential')
                    call WriteTimeseriesHeader(timeseriesUnit,                              &
                      simulationData%TrackingDirection, simulationData%ReferenceTime,       &
                      modelGrid%OriginX, modelGrid%OriginY, modelGrid%RotationAngle)
                    ! Assign writing interface
                    WriteTimeseries => WriteTimeseriesRecordCritical 
                case (2) 
                    ! Parallel consolidated

                    ! Open consolidated direct access output unit
                    ! This case does not allow file header
                    ! reclen = 2I8+es18.9e3+i10+i5+2i10+6es18.9e3+i10+jumpchar
                    reclen = 2*8+18+10+5+2*10+6*18+10+1 
                    open(unit=timeseriesUnit, file=simulationData%TimeseriesFile,     &
                      status='replace', form='formatted', access='direct', recl=reclen)
                    lastRecord = 0

                    ! Initialize binary temporal units 
                    do m = 1, ompNumThreads
                        timeseriesTempUnits( m ) = baseTimeseriesUnit + m
                        open(unit=timeseriesTempUnits( m ), status='scratch', form='unformatted', &
                                                            access='stream', action='readwrite'   )
                    end do
                    timeseriesRecordCounts = 0
                    ! Assign writing interface
                    WriteTimeseries => WriteTimeseriesRecordConsolidate 
                case (3) 
                    ! Parallel not consolidated
                    allocate( timeseriesTempFiles(ompNumThreads) )
                    
                    ! Initialize formatted parallel units 
                    do m = 1, ompNumThreads
                        timeseriesTempUnits( m ) = baseTimeseriesUnit + m
                        write( unit=tempChar, fmt=* )m 
                        write( unit=timeseriesTempFiles( m ), fmt='(a)')&
                            trim(adjustl(tempChar))//'_'//trim(adjustl(simulationData%TimeseriesFile))
                        open( unit=timeseriesTempUnits( m ),     &
                              file=timeseriesTempFiles( m ),     & 
                              status='replace', form='formatted', access='sequential')
                    end do
                    ! Write header to file of first thread 
                    call WriteTimeseriesHeader(timeseriesTempUnits(1),                      &
                      simulationData%TrackingDirection, simulationData%ReferenceTime,       &
                      modelGrid%OriginX, modelGrid%OriginY, modelGrid%RotationAngle)
                    ! Assign writing interface
                    WriteTimeseries => WriteTimeseriesRecordThread
                case default 
                    call ustop('Invalid timeseries output protocol on commmand line. Stop.')
                end select
        end if

    end if


    ! Trace
    if(simulationData%TraceMode .gt. 0) then
        open(unit=traceModeUnit, file=simulationData%TraceFile,                 &
          status='replace', form='formatted', access='sequential')
        write(traceModeUnit, '(1X,A,I10)')                                      &
          'Particle group: ',simulationData%TraceGroup
        write(traceModeUnit, '(1X,A,I10)')                                      &
          'Particle ID: ',simulationData%TraceID
    end if


    ! Open observation cells files
    if ( simulationData%TrackingOptions%observationSimulation ) then
        ! Open the unit and write the header
        do n = 1, simulationData%TrackingOptions%nObservations
            open( unit=simulationData%TrackingOptions%Observations(n)%outputUnit, &
                  file=simulationData%TrackingOptions%Observations(n)%outputFileName,& 
                  status='replace', form='formatted', access='sequential')
            ! For a normal observation
            if ( simulationData%TrackingOptions%Observations(n)%style .eq. 1 )  then 
              ! Remember to replace by binary once in production
              open( unit=simulationData%TrackingOptions%Observations(n)%auxOutputUnit, &
                    file=simulationData%TrackingOptions%Observations(n)%auxOutputFileName,& 
                    status='replace', form='formatted', access='sequential')
            end if
            ! For a sink observation
            if ( simulationData%TrackingOptions%Observations(n)%style .eq. 2 )  then 
              ! Remember to replace by binary once in production
              open( unit=simulationData%TrackingOptions%Observations(n)%auxOutputUnit, &
                    file=simulationData%TrackingOptions%Observations(n)%auxOutputFileName,& 
                    status='replace', form='formatted', access='sequential')
            end if
            !open( unit=simulationData%TrackingOptions%observationUnits(n),     &
            !      file=simulationData%TrackingOptions%observationFiles(n),     & 
            !      status='replace', form='formatted', access='sequential')
            ! Write the corresponding header
            !call WriteObservationHeader(                              &
            !    simulationData%TrackingOptions%observationUnits(n),   &
            !    simulationData%TrackingOptions%observationCells(n),   &
            !    simulationData%ReferenceTime,                         &
            !    modelGrid%OriginX, modelGrid%OriginY,                 &
            !    modelGrid%RotationAngle)
        end do 
    end if 


    ! Begin time step loop
    pathlineRecordCount = 0
    time = 0.0d0
    nt = 0
    if(allocated(cellData)) deallocate(cellData)
    allocate(cellData)
    
    call ulog('Begin TIME_STEP_LOOP', logUnit)
    ! Call system_clock to get the start of the time step loop
    call system_clock(clockCountStart, clockCountRate, clockCountMax)
    TIME_STEP_LOOP: do ktime= kfirst,klast,kincr
    
    ! Get the stress period and time step from the cummulative time step
    call tdisData%GetPeriodAndStep(ktime, period, step)
    
    ! Load data for the current time step
    call flowModelData%LoadTimeStep(period, step)

    ! Load mass transport data for current time step
    if ( simulationData%TrackingOptions%RandomWalkParticleTracking ) then 
        call transportModelData%LoadTimeStep(period, step)
    end if 
    
    if(flowModelData%SteadyState) then
      write(message,'(A,I5,A,I5,A,1PE12.5,A)') 'Processing Time Step ',step,        &
        ' Period ',period,'.  Time = ',tdisData%TotalTimes(ktime), &
        '  Steady-state flow'
    else
      write(message,'(A,I5,A,I5,A,1PE12.5,A)') 'Processing Time Step ',step,        &
        ' Period ',period,'.  Time = ',tdisData%TotalTimes(ktime), &
        '  Transient flow'
    end if
    message = trim(message)
    write(*,'(A)') message
    write(mplistUnit, *)
    write(mplistUnit,'(1X,A)')                                                  &
      '----------------------------------------------------------------------------------------------'
    write(mplistUnit,'(1X,A)') message
    write(mplistUnit,'(1X,A,I6,A)') '  (Cumulative step = ', ktime,')'
    write(mplistUnit,'(1X,A)')                                                  &
      '----------------------------------------------------------------------------------------------'
    
    ! Check water balance summary for the current time step
    if(simulationData%BudgetOutputOption .gt. 0)                                &
      call WriteWaterBalanceSummary(mplistUnit, trackingEngine, cellData)
    
    ! Check cell-by-cell budgets for this time step
    if(simulationData%BudgetCellsCount .gt. 0) then
        write(mplistUnit, *) 
        write(mplistUnit, '(1X,A,I10,A)') 'Cell data will be printed for',      &
          simulationData%BudgetCellsCount, ' cells.'
        do n = 1, simulationData%BudgetCellsCount
            call trackingEngine%FillCellBuffer(simulationData%BudgetCells(n),   &
              cellData)
            call trackingEngine%WriteCellBuffer(mplistUnit, cellData,           &
              simulationData%TrackingOptions%BackwardTracking)
        end do
    end if
    
    ! Compute the tracking time corresponding to the end or beginning 
    ! of this MODFLOW time step (depending on whether this is a forward 
    ! or backward tracking run.)
    message = 'Compute TSMAX'
    call ulog(message, logUnit)
    if(simulationData%TrackingDirection .eq. 1) then
      ! Forward trackine
      tsMax = tdisData%TotalTimes(ktime) - simulationData%ReferenceTime
      if(simulationData%StoppingTimeOption .eq. 2) then
          if(ktime .eq. tdisData%CumulativeTimeStepCount) tsMax = stoptime
      else if(simulationData%StoppingTimeOption .eq. 3) then
          if(ktime .eq. tdisData%CumulativeTimeStepCount) then
              tsMax = stoptime
          else
              if(tsMax .gt. stoptime) tsMax = stoptime    
          end if 
      end if
    else
      ! Backward tracking
      if(ktime .gt. 1) then
          tsMax = simulationData%ReferenceTime - tdisData%TotalTimes(ktime-1)
      else
          tsMax = simulationData%ReferenceTime
      end if
      if(simulationData%StoppingTimeOption .eq. 2) then
          if(ktime .eq. 1) tsMax = stoptime
      else if(simulationData%StoppingTimeOption .eq. 3) then
          if(ktime .eq. 1) then
              tsMax = stoptime
          else
              if(tsMax .gt. stoptime) tsMax = stoptime    
          end if 
      end if
    end if
    
    ! If simulation type is TIMESERIES, write initial locations of all particles active at tracking time = 0,
    ! or all particles regardless of status if that option is set
    ! RWPT
    if((simulationData%SimulationType .ge. 3) .and. & 
      (simulationData%SimulationType .lt. 7) .and. (ktime .eq. kfirst) ) then
      do groupIndex =1, simulationData%ParticleGroupCount
        do particleIndex = 1, simulationData%ParticleGroups(groupIndex)%TotalParticleCount
          p => simulationData%ParticleGroups(groupIndex)%Particles(particleIndex)
          if(((p%Status .eq. 0) .and. (p%InitialTrackingTime .eq. 0.0d0)) .or.    &
            (simulationData%TimeseriesOutputOption .eq. 1) ) then
            pCoord%CellNumber = p%CellNumber
            pCoord%Layer = p%Layer
            pCoord%LocalX = p%LocalX
            pCoord%LocalY = p%LocalY
            pCoord%LocalZ = p%LocalZ
            pCoord%TrackingTime = 0.0d0
            call modelGrid%ConvertToModelXYZ(pCoord%CellNumber,        &
              pCoord%LocalX, pCoord%LocalY, pCoord%LocalZ,            &
              pCoord%GlobalX, pCoord%GlobalY, pCoord%GlobalZ)
            p%InitialGlobalZ = pCoord%GlobalZ
            p%GlobalZ = p%InitialGlobalZ
            p%GlobalX = pCoord%GlobalX ! GPKDE
            p%GlobalY = pCoord%GlobalY ! GPKDE

            if ( .not. &
              simulationData%TrackingOptions%gpkdeSkipTimeseriesWriter ) then 
              ! If parallel:
              !     - With consolidated output initial positions are 
              !       stored into temporal unit for first thread 
              !       and then consolidated
              !     - Without consolidation then initial positions 
              !       will be stored in file for first thread
              ! If not parallel: remains the same as usual
              call WriteTimeseries(p%SequenceNumber, p%ID, groupIndex, & 
                             ktime, 0, pCoord, geoRef, timeseriesUnit, & 
                          timeseriesRecordCounts, timeseriesTempUnits  )

            end if
          end if
        end do
      end do


      ! GPKDE reconstruction for initial 
      ! particles distribution 
      ! reconstruction grouped by solute
      if ( simulationData%TrackingOptions%GPKDEReconstruction ) then

        do ns=1,transportModelData%nSolutes

          solute => transportModelData%Solutes(ns)

          ! Count how many active, considering
          ! all the particle groups linked to a given solute 
          activeCounter = 0
          do npg=1,solute%nParticleGroups
            groupIndex = solute%pGroups(npg)
            do particleIndex = 1, simulationData%ParticleGroups(groupIndex)%TotalParticleCount
              p => simulationData%ParticleGroups(groupIndex)%Particles(particleIndex)
              ! If active, to the array for GPKDE
              !if( (p%Status .eq. 1) ) then
              if((p%Status .eq. 0) .and. (p%InitialTrackingTime .eq. 0.0d0)) then 
                activeCounter = activeCounter + 1
              end if
            end do
          end do

          ! Allocate active particles coordinates
          if ( allocated( activeParticleCoordinates ) ) deallocate( activeParticleCoordinates )
          allocate( activeParticleCoordinates(activeCounter,3) )
          activeParticleCoordinates = 0d0
          if ( allocated( activeParticleMasses ) ) deallocate( activeParticleMasses )
          allocate( activeParticleMasses(activeCounter) )
          activeParticleMasses = 0d0

          ! Could be parallelized ?
          ! Restart active counter and fill coordinates array
          activeCounter = 0 
          do npg=1,solute%nParticleGroups
            groupIndex = solute%pGroups(npg)
            do particleIndex = 1, simulationData%ParticleGroups(groupIndex)%TotalParticleCount
              p => simulationData%ParticleGroups(groupIndex)%Particles(particleIndex)
              ! If active, to the array for GPKDE
              !if( (p%Status .eq. 1) ) then
              if((p%Status .eq. 0) .and. (p%InitialTrackingTime .eq. 0.0d0)) then 
                activeCounter = activeCounter + 1
                activeParticleCoordinates( activeCounter, 1 ) = p%GlobalX
                activeParticleCoordinates( activeCounter, 2 ) = p%GlobalY
                activeParticleCoordinates( activeCounter, 3 ) = p%GlobalZ
                activeParticleMasses( activeCounter ) = p%Mass
              end if
            end do
          end do

          ! GPKDE
          ! Compute density for the particles linked to a given 
          ! solute. These may have different mass
          call gpkde%ComputeDensity(                                              &
             activeParticleCoordinates,                                           &
             outputFileUnit     = simulationData%TrackingOptions%gpkdeOutputUnit, &
             outputDataId       = 0,                                              & ! timeindex
             particleGroupId    = solute%id,                                      &
             unitVolume         = .true.,                                         &
             weightedHistogram  = .true.,                                         &
             weights            = activeParticleMasses                            &
          )

          ! Reconstruction still needs some review. 
          ! When giving an initial condition for a quasi-2D layer
          ! and the dimension in the "compressed" dimension is 
          ! non-zero, the output of gpkde needs to be normalized 
          ! by this distance. Verify if this happens for other
          ! conditions. 

        end do 

      end if

    end if
  

    ! TRACKING_INTERVAL_LOOP: 
    ! Loop through all the required time points that fall within the
    ! current MODFLOW time step. For runs that do not have any specified 
    ! time points, there will only be one time point that corresponds either 
    ! to the beginning or end of the current MODFLOW time step or to the 
    ! specified stop time for the MODPATH analysis.
    itend = 0
    itcount = 0
    call ulog('Begin TRACKING_INTERVAL_LOOP', logUnit)
    TRACKING_INTERVAL_LOOP: do while (itend .eq. 0)
    itcount = itcount + 1 
    print *, itcount, '-----------------------------------------------------------------------------------'

    itend = 1
    maxTime = tsMax
    isTimeSeriesPoint = .false.
    ! RWPT
    if( (simulationData%SimulationType .gt. 1) .and. (simulationData%SimulationType .lt. 7) ) then     
        ! For timeseries and pathline runs, find out if maxTime should be set to the value of the
        ! next time point or the time at the end of the time step
        if (nt+1 .le. simulationData%TimePointCount) then
            if (simulationData%TimePoints(nt+1) .le. tsMax) then
              nt = nt + 1
              maxTime = simulationData%TimePoints(nt)
              tPoint(1) = maxTime
              itend = 0
              if(maxTime .eq. tsMax) itend = 1
              isTimeSeriesPoint = .true.
            end if
        end if
    end if
   
    ! Track particles
    pendingCount = 0
    activeCount = 0
    if(simulationData%ParticleGroupCount .gt. 0) then
        ! -- Particles groups loop -- !
        do groupIndex = 1, simulationData%ParticleGroupCount
            if ( simulationData%SolutesOption .eq. 1 ) then 
                ! Assign pointers to dispersivities 
                ! in transportModelData
                call transportModelData%SetSoluteDispersion( &
                  simulationData%ParticleGroups(groupIndex)%Solute )
                ! Update dispersion function interface
                ! depending on dispersion model 
                call trackingEngine%UpdateDispersionFunction( &
                  transportModelData%Solutes(&
                  simulationData%ParticleGroups(groupIndex)%Solute )%dispersionModel )
            end if 
            ! -- Particles loop -- !
            !$omp parallel do schedule( dynamic,1 )          &
            !$omp default( none )                            &
            !$omp shared( simulationData, modelGrid )        &
            !$omp shared( flowModelData )                    &
            !$omp shared( geoRef )                           &
            !$omp shared( pathlineUnit, binPathlineUnit )    &
            !$omp shared( timeseriesUnit, traceModeUnit )    &
            !$omp shared( timeseriesTempUnits )              &
            !$omp shared( timeseriesRecordCounts )           &
            !$omp shared( timeseriesTempFiles )              &
            !$omp shared( period, step, ktime, nt )          &
            !$omp shared( time, maxTime, isTimeSeriesPoint ) &
            !$omp shared( tPoint, tPointCount )              &
            !$omp shared( groupIndex )                       &
            !$omp private( p, traceModeOn )                  &
            !$omp private( topActiveCellNumber )             &
            !$omp private( pLoc, plCount, tsCount )          &
            !$omp private( pCoordLast, pCoordFirst )         &
            !$omp private( pCoordTP, pCoord )                &
            !$omp private( trackPathResult, status )         &
            !$omp private( timeseriesRecordWritten )         &
            !$omp private( ompThreadId )                     &
            !$omp private( cellDataBuffer, obs, nobs )       &
            !$omp firstprivate( trackingEngine )             &
            !$omp firstprivate( WriteTimeseries )            &
            !$omp reduction( +:pendingCount )                &
            !$omp reduction( +:activeCount )                 &
            !$omp reduction( +:pathlineRecordCount ) 
            do particleIndex = 1, simulationData%ParticleGroups(groupIndex)%TotalParticleCount
                timeseriesRecordWritten = .false.
                p => simulationData%ParticleGroups(groupIndex)%Particles(particleIndex)
!!                ! Check particle status. 
!!                ! Skip over particles unless they are active or pending release.
!!                if(p%Status .gt. 1) then
!!                    ! Add code here later to deal with advective observations
!!                    ! For now, just cycle to the next particle
!!                    cycle
!!                end if
                
                ! Check to see if trace mode should be turned on for this particle
                traceModeOn = .false.
                if(simulationData%TraceMode .gt. 0) then 
                    if((p%Group .eq. simulationData%TraceGroup) .and.           &
                      (p%ID .eq. simulationData%TraceID)) traceModeOn = .true.
                end if
                
                ! If a particle is pending release (STATUS = 0), check to see if it should
                ! be set to active and released on this pass. If the particle is pending
                ! release and its release time is earlier than the starting time of this
                ! pass, then mark the particle status as permanently unreleased 
                ! (STATUS = 8).
                if(p%Status .eq. 0) then
                    if(p%InitialTrackingTime .lt. time) then
                        p%Status = 8
                    else if(p%InitialTrackingTime .le. maxTime) then
                        p%Status = 1
                        if(p%Drape .eq. 0) then
                            ! Drape option is not in effect.
                            if(trackingEngine%FlowModelData%IBoundTS(p%CellNumber) .eq. 0) then
                                p%Status = 7
                            end if
                        else
                            ! Drape option is in effect. Find the top-most active cell starting with the initial cell number. 
                            ! If no active cell is found, leave the cell number set to its original value and set the Status = 7
                            ! to indicate it is stranded in an inactive cell.
                            topActiveCellNumber = trackingEngine%GetTopMostActiveCell(p%CellNumber)
                            if(topActiveCellNumber .gt. 0) then
                                p%CellNumber = topActiveCellNumber
                            else
                                p%Status = 7
                            end if
                        end if
                        call modelGrid%ConvertToModelZ(p%InitialCellNumber, &
                          p%InitialLocalZ, p%InitialGlobalZ, .true.)
                        p%GlobalZ = p%InitialGlobalZ
                    end if
                end if
        
                ! RWPT
                ! Verify cell not dry anymore
                if (p%Status .eq. 7 ) then 
                    ! Initialize cellBuffer cellNumber
                    call trackingEngine%FillCellBuffer( p%CellNumber, cellDataBuffer )

                    ! Verify dry/partially dried cells
                    call cellDataBuffer%VerifyDryCell()

                    ! If partially dried restore active/track status, otherwise keep Status = 7
                    if ( cellDataBuffer%partiallyDry ) p%Status = 1 ! Track particle

                end if 

                ! Count the number of particles that are currently active or pending
                ! release at the beginning of this pass.         
                if(p%Status .EQ. 0) pendingCount = pendingCount + 1
                if(p%Status .EQ. 1) activeCount = activeCount + 1
                
                ! Track the particle if it is active
                if(p%Status .eq. 1) then
                    ! Set particle location buffer
                    pLoc%CellNumber = p%CellNumber
                    pLoc%Layer = p%Layer
                    pLoc%LocalX = p%LocalX
                    pLoc%LocalY = p%LocalY
                    pLoc%LocalZ = p%LocalZ
                    pLoc%TrackingTime = p%TrackingTime
                   
                    ! Call TrackPath
                    call trackingEngine%TrackPath(trackPathResult, traceModeOn, &
                      traceModeUnit, p%Group, p%ID, p%SequenceNumber, pLoc,     &
                      maxTime, tPoint, tPointCount)
                    
                    ! Update endpoint data. The Face property will only be updated when the endpoint file is written
                    plCount = trackPathResult%ParticlePath%Pathline%GetItemCount()
                    tsCount = trackPathResult%ParticlePath%Timeseries%GetItemCount()
                    pCoordLast => trackPathResult%ParticlePath%Pathline%Items(plCount)
                    pCoordFirst => trackPathResult%ParticlePath%Pathline%Items(1)
                    p%CellNumber =  pCoordLast%CellNumber
                    p%Layer = pCoordLast%Layer
                    p%LocalX = pCoordLast%LocalX
                    p%LocalY = pCoordLast%LocalY
                    p%LocalZ = pCoordLast%LocalZ
                    p%GlobalZ = pCoordLast%GlobalZ
                    p%TrackingTime = pCoordLast%TrackingTime
                    
                    ! Update particle status
                    status = trackPathResult%Status
                    if(  status .eq. trackPathResult%Status_ReachedBoundaryFace()) then
                        p%Status = 2
                    else if(status .eq. trackPathResult%Status_StopAtWeakSink()) then
                        p%Status = 3
                    else if(status .eq. trackPathResult%Status_StopAtWeakSource()) then
                        p%Status = 4
                    else if(status .eq. trackPathResult%Status_NoExitPossible()) then
                        p%Status = 5
                    else if(status .eq. trackPathResult%Status_StopZoneCell()) then
                        p%Status = 6
                    else if(status .eq. trackPathResult%Status_InactiveCell()) then
                        p%Status = 7
                    else if(status .eq. trackPathResult%Status_Undefined()) then
                        p%Status = 9
                    else
                        ! Leave status set to active (status = 1)
                    end if
                    
                    ! Write particle output
                    if((simulationData%SimulationType .eq. 2)  .or.              &
                       (simulationData%SimulationType .eq. 4)  .or.              & 
                       (simulationData%SimulationType .eq. 6)) then
                        ! Write pathline to pathline file
                        if(plCount .gt. 1) then
                            pathlineRecordCount = pathlineRecordCount + 1
                            !$omp critical (pathline)
                            select case (simulationData%PathlineFormatOption)
                                case (1)
                                    call WriteBinaryPathlineRecord(             &
                                      trackPathResult, binPathlineUnit, period, &
                                      step, geoRef)
                                case (2)
                                    call WritePathlineRecord(trackPathResult,   &
                                      pathlineUnit, period, step, geoRef)
                            end select
                            !$omp end critical (pathline)
                        end if
                    end if
                    if(simulationData%SimulationType .ge. 3) then
                      if(tsCount .gt. 0) then
                        ! Write timeseries record to the timeseries file
                        pCoordTP => trackPathResult%ParticlePath%Timeseries%Items(1)
                        p%GlobalX = pCoordTP%GlobalX ! GPKDE
                        p%GlobalY = pCoordTP%GlobalY ! GPKDE
                        if ( .not. &
                          simulationData%TrackingOptions%gpkdeSkipTimeseriesWriter ) then 
                            ! With interface
                            call WriteTimeseries(p%SequenceNumber, p%ID, groupIndex, & 
                                        ktime, nt, pCoordTP, geoRef, timeseriesUnit, & 
                                        timeseriesRecordCounts, timeseriesTempUnits  )
                        end if 
                        timeseriesRecordWritten = .true. ! ?
                        
                        if ( simulationData%anyObservation ) then  
                          ! Write record for resident observations
                          if ( &
                            simulationData%TrackingOptions%isObservation(pCoordTP%CellNumber) ) then 
                            obs => simulationData%TrackingOptions%Observations(&
                                simulationData%TrackingOptions%idObservation(pCoordTP%CellNumber) )
                            if ( obs%style .eq. 1 ) then  
                              do nobs=1,obs%nCells
                                if( obs%cells(nobs) .ne. pCoordTP%CellNumber ) cycle
                                ! If it is part of the cells in the obs, write
                                ! record to auxOutputUnit
                                ! Temp proxy !
                                call WriteTimeseriesRecordCritical(& 
                                    p%SequenceNumber, p%ID, groupIndex, ktime, &
                                      nt, pCoordTP, geoRef, obs%auxOutputUnit, & 
                                  timeseriesRecordCounts, timeseriesTempUnits  )
                              end do 
                            end if
                          end if
                        end if

                      end if
                    end if
                end if
                
                ! If option is set to write timeseries records for all particles
                ! regardless of status, write the record if not done already.
                ! RWPT
                if((simulationData%SimulationType .ge. 3) .and.              &
                   (simulationData%SimulationType .lt. 7) .and.              &
                   (simulationData%TimeseriesOutputOption .eq. 1) .and.      &
                   (isTimeSeriesPoint) .and. (.not.timeseriesRecordWritten)) then
                      pCoord%CellNumber = p%CellNumber
                      pCoord%Layer = p%Layer
                      pCoord%LocalX = p%LocalX
                      pCoord%LocalY = p%LocalY
                      pCoord%LocalZ = p%LocalZ
                      pCoord%TrackingTime = maxTime
                      call modelGrid%ConvertToModelXYZ(pCoord%CellNumber,       &
                        pCoord%LocalX, pCoord%LocalY, pCoord%LocalZ,            &
                        pCoord%GlobalX, pCoord%GlobalY, pCoord%GlobalZ)
                      p%GlobalX = pCoord%GlobalX ! GPKDE
                      p%GlobalY = pCoord%GlobalY ! GPKDE
                      if ( .not. &
                        simulationData%TrackingOptions%gpkdeSkipTimeseriesWriter ) then 
                          ! With interface
                          call WriteTimeseries(p%SequenceNumber, p%ID, groupIndex, & 
                                        ktime, nt, pCoord, geoRef, timeseriesUnit, &
                                        timeseriesRecordCounts, timeseriesTempUnits)
                      end if

                      if ( simulationData%anyObservation ) then  
                        ! Write record for resident observations
                        if ( &
                          simulationData%TrackingOptions%isObservation(pCoordTP%CellNumber) ) then 
                          obs => simulationData%TrackingOptions%Observations(&
                              simulationData%TrackingOptions%idObservation(pCoordTP%CellNumber) )
                          if ( obs%style .eq. 1 ) then  
                            do nobs=1,obs%nCells
                              if( obs%cells(nobs) .ne. pCoordTP%CellNumber ) cycle
                              ! If it is part of the cells in the obs, write
                              ! record to auxOutputUnit
                              call WriteTimeseriesRecordCritical(& 
                                  p%SequenceNumber, p%ID, groupIndex, ktime, &
                                    nt, pCoordTP, geoRef, obs%auxOutputUnit, & 
                                timeseriesRecordCounts, timeseriesTempUnits  )
                            end do 
                          end if
                        end if
                      end if

                end if

            end do
            !$omp end parallel do

        end do


        ! Once it finished transporting all 
        ! particle groups, reconstruction 
        ! grouped by solute
        if ( simulationData%TrackingOptions%GPKDEReconstruction .and. isTimeSeriesPoint ) then

          do ns=1,transportModelData%nSolutes

            solute => transportModelData%Solutes(ns)

            ! Count how many active, considering
            ! all the particle groups linked to a given solute 
            activeCounter = 0
            do npg=1,solute%nParticleGroups
              groupIndex = solute%pGroups(npg)
              do particleIndex = 1, simulationData%ParticleGroups(groupIndex)%TotalParticleCount
                p => simulationData%ParticleGroups(groupIndex)%Particles(particleIndex)
                ! If active, to the array for GPKDE
                if( (p%Status .eq. 1) ) then
                  activeCounter = activeCounter + 1
                end if
              end do
            end do

            ! Allocate active particles coordinates
            if ( allocated( activeParticleCoordinates ) ) deallocate( activeParticleCoordinates )
            allocate( activeParticleCoordinates(activeCounter,3) )
            activeParticleCoordinates = 0d0
            if ( allocated( activeParticleMasses ) ) deallocate( activeParticleMasses )
            allocate( activeParticleMasses(activeCounter) )
            activeParticleMasses = 0d0

            ! Could be parallelized ?
            ! Restart active counter and fill coordinates array
            activeCounter = 0 
            do npg=1,solute%nParticleGroups
              groupIndex = solute%pGroups(npg)
              do particleIndex = 1, simulationData%ParticleGroups(groupIndex)%TotalParticleCount
                p => simulationData%ParticleGroups(groupIndex)%Particles(particleIndex)
                ! If active, to the array for GPKDE
                if( (p%Status .eq. 1) ) then
                  activeCounter = activeCounter + 1
                  activeParticleCoordinates( activeCounter, 1 ) = p%GlobalX
                  activeParticleCoordinates( activeCounter, 2 ) = p%GlobalY
                  activeParticleCoordinates( activeCounter, 3 ) = p%GlobalZ
                  activeParticleMasses( activeCounter ) = p%Mass
                end if
              end do
            end do

            ! GPKDE
            ! Compute density for the particles linked to a given 
            ! solute. These may have different mass
            call gpkde%ComputeDensity(                                              &
               activeParticleCoordinates,                                           &
               outputFileUnit     = simulationData%TrackingOptions%gpkdeOutputUnit, &
               outputDataId       = nt,                                             & ! timeindex
               particleGroupId    = solute%id,                                      &
               unitVolume         = .true.,                                         &
               weightedHistogram  = .true.,                                         &
               weights            = activeParticleMasses                            &
            )
            
            ! And needs volume correction for 
            ! cells, considering porosities and so on
            ! It would work only for structured grids
            ! sharing the same discretization than GPKDE

            ! if( modelGrid%isUnstructured ) then 
            !   ! It could do some kind of histogram reconstruction
            ! else
            !  ! Run over the active 
            !  read(inUnit, *) layer, row, column
            !  call ugetnode(layerCount, rowCount, columnCount, layer, row, column,cellNumber)
            !  obs%cells(no) = cellNumber
            ! end if

          end do 

        end if

    end if
 

    ! Observation cells: if any sink cell observation, 
    ! writing flow rates might be done here.
    ! This works for the approach were the timeseries
    ! points determine the obs records
    if( ( simulationData%TrackingOptions%anySinkObservation ) & 
                                    .and. isTimeSeriesPoint ) then
        ! Write sink flow rates only for obs of this kind
        do nobs=1,simulationData%TrackingOptions%nObservations
          obs =>simulationData%TrackingOptions%Observations(nobs)
          if ( obs%style .ne. 2 ) cycle
          if (allocated(qSinkBuffer))deallocate(qSinkBuffer)
          allocate(qSinkBuffer(obs%nCells))
          do n=1,obs%nCells
            qSinkBuffer(n) = flowModelData%SinkFlows(obs%cells(n))
          end do
          write(qSinkFormat,*) '(1I10,', obs%nCells, 'es18.9e3)'
          write(obs%auxOutputUnit,qSinkFormat) nt, qSinkBuffer(:)
        end do

    end if


    ! If timeseries simulation and parallel and consolidated output
    ! Consolidation is done at this stage to preserve
    ! sorting of time indexes
    if( ((simulationData%SimulationType .eq. 3)  .or. (simulationData%SimulationType .eq. 4) ) .and. & 
        parallel .and. (tsOutputType .eq. 2) ) then
        if ( .not. all( timeseriesRecordCounts .eq. 0 ) ) then
            call ConsolidateParallelTimeseriesRecords( timeseriesTempUnits, timeseriesUnit, timeseriesRecordCounts, lastRecord )
            ! Restart temporal binary units and record counters
            do m = 1, ompNumThreads
                if ( timeseriesRecordCounts(m) .gt. 0 ) then
                    close( timeseriesTempUnits( m ) )
                    open(unit=timeseriesTempUnits( m ), status='scratch', form='unformatted', &
                                                        access='stream', action='readwrite'   )
                end if
            end do
            timeseriesRecordCounts = 0 
        end if
    end if

    ! Update tracking time
    time = maxTime
    
    ! Check to see if there are any particles remaining to track in the next
    ! pass. If not, exit the loop.
    if(simulationData%ParticleGroupCount .gt. 0) then
        if(activeCount .eq. 0 .and. pendingCount .eq. 0) then
            call ulog('No active particles remain. Exit TRACKING_INTERVAL_LOOP.', logUnit)
            exit TIME_STEP_LOOP
        end if
    end if

    end do TRACKING_INTERVAL_LOOP   
    call ulog('Exit TRACKING_INTERVAL_LOOP', logUnit)
       
    ! Exit TIME_STEP_LOOP if the tracking time has reached the specified stop time.
    IF(time .ge. stoptime) exit TIME_STEP_LOOP

    end do TIME_STEP_LOOP
    call ulog('Exit TIME_STEP_LOOP', logUnit)
    
    ! Stop timer
    call system_clock(clockCountStop, clockCountRate, clockCountMax)
    
    ! Write endpoint file
    if(simulationData%ParticleGroupCount .gt. 0) then
        call ulog('Write endpoint file.', logUnit)
        call WriteEndpoints(simulationData, modelGrid, geoRef, endpointUnit)
    end if
    
    ! Finalize and process binary pathline file if pathline format option = 1
    if((simulationData%SimulationType .eq. 2) .or. & 
       (simulationData%SimulationType .eq. 4) .or. &
       (simulationData%SimulationType .eq. 6)) then
        if(simulationData%PathlineFormatOption .eq. 1) then
            call ulog('Consolidating pathline segments.', logUnit)
            call ConsolidatePathlines(binPathlineUnit, pathlineUnit,            &
              pathlineRecordCount, simulationData%TotalParticleCount) 
        end if
    end if
 

    ! RWPT
    ! Process observation cells for reconstruction
    if ( simulationData%TrackingOptions%observationSimulation ) then


      ! If so, reset gpkde
      call ulog('Reset and reinitialize GPKDE for observation cells ', logUnit)
      ! If it was allocated from spatial reconstruction
      if( allocated( gpkde ) ) then 
         call gpkde%Reset()
      else
         allocate( gpkde )
      end if  


      ! Loop over observations
      do nobs=1, simulationData%TrackingOptions%nObservations
        obs => simulationData%TrackingOptions%Observations(nobs)

        ! If this is a normal obs cell
        if ( obs%style .eq. 1 ) then 

          ! The length of the timeseries is needed

          ! This would work only for timeseries simulations
          ! with regular timestep definition
          dTObsSeries = simulationData%TimePoints(2)-simulationData%TimePoints(1)  ! binSize

          ! Consider irregular binsize GPKDE Reconstruction

          ! If no timeseries run, it would be possible to 
          ! create an histogram based on stoptime for example 
          ! and given number of bins. Still, this 
          ! should probably be done before, while reading
          ! simulation data, in order to force a timeseries
          ! run, hand craft a timeseries sim writing records
          ! at the specified times

          ! Initialize gpkde for timeseries reconstruction
          call gpkde%Initialize(& 
              (/maxval(simulationData%TimePoints(:)),0d0,0d0/),                  &
              (/dtObsSeries,0d0,0d0/),                                           &
              domainOrigin=(/simulationData%ReferenceTime,0d0,0d0/),             &
              nOptimizationLoops=simulationData%TrackingOptions%gpkdeNOptLoops,  &
              databaseOptimization=.false.                                       &
          )

          ! Then read the observation file
          ! and pass it to gpkde for reconstruction
      
          ! It needs some obs record count or something
          rewind( obs%auxOutputUnit )
          nlines = 0
          do
            read(obs%auxOutputUnit,*,iostat=io)
            if (io/=0) exit
            nlines = nlines + 1
          end do

          ! Allocate active particles coordinates (temporary)
          if ( allocated( activeParticleCoordinates ) ) deallocate( activeParticleCoordinates )
          allocate( activeParticleCoordinates(nlines,2) )
          activeParticleCoordinates = 0d0
          if ( allocated( activeParticleMasses ) ) deallocate( activeParticleMasses )
          allocate( activeParticleMasses(nlines) )
          activeParticleMasses = 0d0


          ! It seems that the most reasonable 
          ! approach would be to read solute id
          ! from records. 


          ! Load file records into array
          rewind( obs%auxOutputUnit )
          do n = 1, nlines
              ! Read from obs file
              ! Based on TS record, it could be reduced 
              read( obs%auxOutputUnit, '(2I8,es18.9e3,i10,i5,2i10,6es18.9e3,i10)')          &
                timePointIndex, timeStep, initialTime, sequenceNumber, groupIndex,  &
                particleID, pCoord%CellNumber, pCoord%LocalX, pCoord%LocalY, pCoord%LocalZ, &
                modelX, modelY, pCoord%GlobalZ, pCoord%Layer
              ! Needs some kind of understanding of the particle group, and that it
              ! means another solute ( column? )
              activeParticleCoordinates(n,1) = initialTime 
              activeParticleCoordinates(n,2) = groupIndex

              ! A similar access could be used for getting soluteId from 
              ! the particle directly, avoiding the identification stage 
              ! coming further down
              activeParticleMasses(n) = &
                simulationData%ParticleGroups(groupIndex)%Particles(particleID)%Mass
          end do 

          ! idColFormat for observations output
          idColFormat = 2

          ! For storing the BTCS
          if (allocated(BTCPerSolute)) deallocate(BTCPerSolute)
          allocate(BTCPerSolute(simulationData%TimePointCount, transportModelData%nSolutes))
          BTCPerSolute = 0d0
          if (idColFormat.eq.2) then 
            if (allocated(BTCHistPerSolute)) deallocate(BTCHistPerSolute)
            allocate(BTCHistPerSolute(simulationData%TimePointCount, transportModelData%nSolutes))
            BTCHistPerSolute = 0d0
          end if

          ! Loop over solutes
          do ns=1, transportModelData%nSolutes
            solute => transportModelData%Solutes(ns)

            ! Count how many for this solute
            solCount = 0
            do npg=1,solute%nParticleGroups
              solCount = solCount + count(activeParticleCoordinates(:,2).eq.solute%pGroups(npg) ) 
            end do
            if ( allocated(gpkdeDataCarrier) ) deallocate(gpkdeDataCarrier) 
            allocate( gpkdeDataCarrier(solCount,3) )
            gpkdeDataCarrier = 0d0
            if ( allocated(gpkdeWeightsCarrier) ) deallocate(gpkdeWeightsCarrier) 
            allocate( gpkdeWeightsCarrier(solCount) )
            gpkdeWeightsCarrier = 0d0

            ! Not necessarily the most efficient,
            ! think about cases with lots of pgroups per
            ! solute. Is either this or write the solute id 
            ! to the observation record. This would require
            ! modification of the TrackPath interfaces at trackingEngine
            irow = 0
            do n=1,nlines
              anyFromThisSolute = .false.
              do npg=1,solute%nParticleGroups
                if (activeParticleCoordinates(n,2).eq.solute%pGroups(npg)) then 
                    anyFromThisSolute = .true.
                    exit
                end if
              end do
              if ( .not. anyFromThisSolute ) cycle
              irow = irow + 1
              gpkdeDataCarrier(irow,1) = activeParticleCoordinates(n,1)
              gpkdeWeightsCarrier(irow) = activeParticleMasses(n)
            end do

            ! Timeseries reconstruction    
            call gpkde%ComputeDensity(   &
              gpkdeDataCarrier,          &
              unitVolume = .true.,       &
              histogramScalingFactor=1d0,&
              weightedHistogram = .true.,&
              weights = gpkdeWeightsCarrier )

            BTCPerSolute(:,ns) = gpkde%densityEstimateGrid(:,1,1)

            if (idColFormat.eq.2) then 
              BTCHistPerSolute(:,ns) = gpkde%rawDensityEstimateGrid(:,1,1)
            end if 

          end do


          ! Accumulate volumes
          ! Compute porous volume for each cell in the 
          ! observation and add them for writing the 
          ! observation record. 
          ! Assume perfectly saturated cell for all times and 
          ! volume is computed with dZ = Top - Bottom
          obsAccumPorousVolume = 0d0
          do n=1, obs%nCells
            dX       = modelGrid%DelX(obs%cells(n))
            dY       = modelGrid%DelY(obs%cells(n))
            dZ       = abs(modelGrid%Top(obs%cells(n))-modelGrid%Bottom(obs%cells(n)))
            porosity = flowModelData%Porosity(obs%cells(n))
            obsAccumPorousVolume = & 
               obsAccumPorousVolume + dX*dY*dZ*porosity 
          end do 

          ! Once the accumulated porous volume is known, compute resident
          ! concentration

          ! Apply the logic to determine where to write the obs records
          ! Case 1: write to the same file as before: close it, open again and dump
          close( obs%outputUnit )
          open( unit=obs%outputUnit, &
                file=obs%outputFileName,& 
                status='replace', form='formatted', access='sequential')
          ! Case 2: close the previous file and open a new one using the same
          ! unit number for the obs, but with a different filename: keep the arrival records
         

          ! And write
          ! Remember to decide what to do for different species, pgroups (?)
          if ( obsAccumPorousVolume .ne. 0d0 ) then 
            ! Probably should be as a property from obs 
            select case(idColFormat)
            case (1)
              ! idTime, time, C-GPKDE(:)
              ! Needs the format 
              write (colFormat,*) '(1I8,',&
                  1 + transportModelData%nSolutes, 'es18.9e3)'
              ! And write
              do nit = 1, simulationData%TimePointCount
                ! idTime, time, C-GPKDE(t,:)
                write(obs%outputUnit, colFormat) nit, simulationData%TimePoints(nit), &
                                            BTCPerSolute(nit,:)/obsAccumPorousVolume
              end do 
            case (2)
              ! idTime, time, C-GPKDE(:), C-HIST(:)
              ! Needs the format 
              write (colFormat,*) '(1I8,',&
                  1 + 2*transportModelData%nSolutes, 'es18.9e3)'
              ! And write
              do nit = 1, simulationData%TimePointCount
                ! idTime, time, C-GPKDE(t,:), C-HIST(t,:)
                write(obs%outputUnit, colFormat) nit, simulationData%TimePoints(nit), &
                      BTCPerSolute(nit,:)/obsAccumPorousVolume, &
                        BTCHistPerSolute(nit,:)/obsAccumPorousVolume
              end do
            case default 
                continue
            end select 

          end if 


          ! And reset gpkde 
          call gpkde%Reset()


        end if ! If obs%style.eq.1


        ! If this is is a sink obs cell
        if ( obs%style .eq. 2 ) then 

          ! The length of the timeseries is needed

          ! This would work only for timeseries simulations
          dTObsSeries = simulationData%TimePoints(2)-simulationData%TimePoints(1)  ! binSize

          ! Consider irregular binsize GPKDE Reconstruction

          ! If no timeseries run, it would be possible to 
          ! create an histogram based on stoptime for example 
          ! and given number of bins

          ! Initialize gpkde for timeseries reconstruction
          call gpkde%Initialize(& 
              (/maxval(simulationData%TimePoints(:)),0d0,0d0/),                  &
              (/dtObsSeries,0d0,0d0/),                                           &
              domainOrigin=(/simulationData%ReferenceTime,0d0,0d0/),             &
              nOptimizationLoops=simulationData%TrackingOptions%gpkdeNOptLoops,  &
              databaseOptimization=.false.                                       &
          )

          ! Then read the observation file to extract arrival times
          ! and pass it to gkde for reconstruction
      
          ! It needs some obs record count or something
          rewind( obs%outputUnit )
          nlines = 0
          do
            read(obs%outputUnit,*,iostat=io)
            if (io/=0) exit
            nlines = nlines + 1
          end do


          ! Allocate active particles coordinates (temporary)
          if ( allocated( activeParticleCoordinates ) ) deallocate( activeParticleCoordinates )
          allocate( activeParticleCoordinates(nlines,2) )
          activeParticleCoordinates = 0d0
          if ( allocated( activeParticleMasses ) ) deallocate( activeParticleMasses )
          allocate( activeParticleMasses(nlines) )
          activeParticleMasses = 0d0

          ! Load file records into array
          rewind( obs%outputUnit )
          do n = 1, nlines
              ! Read from obs file
              read( obs%outputUnit, '(3I8,5es18.9e3)' ) &
                groupIndex, particleID, cellNumber,     &
                initialTime, initialGlobalX, initialGlobalY, initialGlobalZ, QSinkCell
              ! Needs some kind of understanding of the particle group, and that it
              ! means another solute ( column? )
              activeParticleCoordinates(n,1) = initialTime 
              activeParticleCoordinates(n,2) = groupIndex ! remember these indexes

              ! A similar access could be used for getting soluteId from 
              ! the particle directly, avoiding the identification stage 
              ! coming further down
              activeParticleMasses(n) = &
                simulationData%ParticleGroups(groupIndex)%Particles(particleID)%Mass
          end do 

          ! idColFormat for observations output
          idColFormat = 2

          ! For storign the BTCS
          if (allocated(BTCPerSolute)) deallocate(BTCPerSolute)
          allocate(BTCPerSolute(simulationData%TimePointCount, transportModelData%nSolutes))
          BTCPerSolute = 0d0
          if (idColFormat.eq.2) then 
            if (allocated(BTCHistPerSolute)) deallocate(BTCHistPerSolute)
            allocate(BTCHistPerSolute(simulationData%TimePointCount, transportModelData%nSolutes))
          end if
          BTCHistPerSolute = 0d0

          ! Loop over solutes
          do ns=1, transportModelData%nSolutes
            solute => transportModelData%Solutes(ns)

            ! Count how many for this solute
            solCount = 0
            do npg=1,solute%nParticleGroups
              solCount = solCount + count(activeParticleCoordinates(:,2).eq.solute%pGroups(npg)) 
            end do
            if ( allocated(gpkdeDataCarrier) ) deallocate(gpkdeDataCarrier) 
            allocate( gpkdeDataCarrier(solCount,3) )
            gpkdeDataCarrier = 0d0
            if ( allocated(gpkdeWeightsCarrier) ) deallocate(gpkdeWeightsCarrier) 
            allocate( gpkdeWeightsCarrier(solCount) )
            gpkdeWeightsCarrier = 0d0

            irow = 0
            do n=1,nlines
              anyFromThisSolute = .false.
              do npg=1,solute%nParticleGroups
                if (activeParticleCoordinates(n,2).eq.solute%pGroups(npg)) then 
                    anyFromThisSolute = .true.
                    exit
                end if
              end do
              if ( .not. anyFromThisSolute ) cycle
              irow = irow + 1
              gpkdeDataCarrier(irow,1) = activeParticleCoordinates(n,1)
              gpkdeWeightsCarrier(irow) = activeParticleMasses(n)
            end do

            ! Timeseries reconstruction    
            call gpkde%ComputeDensity(   &
              gpkdeDataCarrier,          &
              unitVolume = .false.,      &
              histogramScalingFactor=1d0,&
              weightedHistogram = .true.,&
              weights = gpkdeWeightsCarrier )

            BTCPerSolute(:,ns) = gpkde%densityEstimateGrid(:,1,1)

            if (idColFormat.eq.2) then 
              BTCHistPerSolute(:,ns) = gpkde%rawDensityEstimateGrid(:,1,1)
            end if 

          end do 
    

          ! Flow rates were written to obs file
          ! And are now sorted by arrival time
          if ( allocated(obsSinkFlowInTime) ) deallocate(obsSinkFlowInTime)
          ! With as many columns as cells
          ! composing the observation
          allocate(obsSinkFlowInTime(simulationData%TimePointCount, obs%nCells))
          obsSinkFlowInTime = 0d0

          ! Fill flow-rate timeseries for each cell
          rewind( obs%auxOutputUnit ) 
          do n =1, simulationData%TimePointCount
            read(obs%auxOutputUnit,*) timeIndex, obsSinkFlowInTime( n, : )
          end do

          ! Accumulate flow rates, absolute values 
          obsAccumSinkFlowInTime = sum( abs(obsSinkFlowInTime), dim=2 )
          ! Once flow-rates are known, can compute flux-concentration

          ! Something to verify if sink flows were zero the whole time

          ! Apply the logic to determine where to write the obs records
          ! Case 1: write to the same file as before: close it, open again and dump
          close( obs%outputUnit )
          open( unit=obs%outputUnit, &
                file=obs%outputFileName,& 
                status='replace', form='formatted', access='sequential')
          ! Case 2: close the previous file and open a new one using the same
          ! unit number for the obs, but with a different filename: keep the arrival records
      
          ! Probably should be as a property from  obs 
          select case(idColFormat)
          case (1)
            ! idTime, time, QSink, CFlux-GPKDE(:)
            ! Needs the format 
            write (colFormat,*) '(1I8,',&
                2 + transportModelData%nSolutes, 'es18.9e3)'
            ! And write
            do nit = 1, simulationData%TimePointCount
               if ( obsAccumSinkFlowInTime(nit) .gt. 0d0 ) then 
                  ! idTime, time, QSink, CFlux-GPKDE(t,:)
                  write(obs%outputUnit, colFormat) nit, simulationData%TimePoints(nit), &
                        obsAccumSinkFlowInTime(nit), BTCPerSolute(nit,:)/obsAccumSinkFlowInTime(nit)
               else
                 ! No concentrations
                 ! idTime, time, QSink, CFlux-GPKDE(t,:), CFlux-HIST(t,:)
                 write(obs%outputUnit, colFormat) nit, simulationData%TimePoints(nit), &
                       0d0, spread(0d0,1,transportModelData%nSolutes) 
               end if 
            end do 
          case (2)
            ! idTime, time, QSink, CFlux-GPKDE(:), CFlux-HIST(:)
            ! Needs the format 
            write (colFormat,*) '(1I8,',&
                2 + 2*transportModelData%nSolutes, 'es18.9e3)'
            ! And write
            do nit = 1, simulationData%TimePointCount
               if ( obsAccumSinkFlowInTime(nit) .gt. 0d0 ) then 
                  ! idTime, time, QSink, CFlux-GPKDE(t,:), CFlux-HIST(t,:)
                  write(obs%outputUnit, colFormat) nit, simulationData%TimePoints(nit), &
                        obsAccumSinkFlowInTime(nit), BTCPerSolute(nit,:)/obsAccumSinkFlowInTime(nit), &
                                                     BTCHistPerSolute(nit,:)/obsAccumSinkFlowInTime(nit) 
               else
                 ! No concentrations
                 ! idTime, time, QSink, CFlux-GPKDE(t,:), CFlux-HIST(t,:)
                 write(obs%outputUnit, colFormat) nit, simulationData%TimePoints(nit), &
                       0d0, spread(0d0,1,2*transportModelData%nSolutes) 
               end if 
            end do
          case default 
              continue
          end select 


          ! And reset gpkde 
          call gpkde%Reset()


        end if ! If obs%style.eq.2


      end do ! obsLoop


    end if ! process obs cells 


    ! Write particle summary information
    call WriteParticleSummaryInfo(simulationData, mplistUnit)


100 continue    

    ! RWPT
    ! Close observation units if any
    if ( simulationData%TrackingOptions%observationSimulation ) then
        do n = 1, simulationData%TrackingOptions%nObservations
            close( simulationData%TrackingOptions%Observations(n)%outputUnit )
            if ( simulationData%TrackingOptions%Observations(n)%style .eq. 1 ) then 
              close( simulationData%TrackingOptions%Observations(n)%auxOutputUnit )
            end if 
            if ( simulationData%TrackingOptions%Observations(n)%style .eq. 2 ) then 
              close( simulationData%TrackingOptions%Observations(n)%auxOutputUnit )
            end if 
        end do 
        ! And deallocate arrays of observation information
        call simulationData%TrackingOptions%Reset()
    end if 


    ! Deallocate major components
    call ulog('Begin memory deallocation.', logUnit)
    if(allocated(headReader)) deallocate(headReader)
    if(allocated(budgetReader)) deallocate(budgetReader)
    if(allocated(tdisData)) deallocate(tdisData)
    if(allocated(trackingEngine)) deallocate(trackingEngine)
    if(allocated(flowModelData)) deallocate(flowModelData)
    if(allocated(transportModelData)) deallocate(transportModelData)
    if(allocated(basicData)) deallocate(basicData)
    if(allocated(simulationData)) deallocate(simulationData)
    ! GPKDE
    if(allocated(gpkde)) deallocate(gpkde)
    call ulog('Memory deallocation complete.', logUnit)
    
    write(*, '(a)') terminationMessage
    write(mplistUnit, '(1x/,a)', err=200) terminationMessage
    elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
    write(mplistUnit, '(1X,A,E15.5,A)') 'Elapsed time = ', elapsedTime, ' seconds'
    write(*, '(1X,A,E15.5,A)') 'Elapsed time = ', elapsedTime, ' seconds'

    ! Close files
200 continue    
    close(mplistUnit)
    if (logType /= 0) close(logUnit)

    ! Uncomment the following pause statement when running in debug mode within Visual Studio.
    ! The pause statement keeps the command window from immediately closing when the MODPATH run completes.
    ! For a compiled executable, do not use the pause statement, but run the executable from inside a batch
    ! file to keep the window from closing immediately.
    
    !pause
    
    contains



    subroutine ParseCommandLine(mpsimFile, mplogFile, logType, parallel, tsOutputType)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
    implicit none
    character*(*),intent(inout) :: mpsimFile
    character*(*),intent(inout) :: mplogFile
    integer,intent(inout) :: logType
    logical,intent(inout) :: parallel
    integer,intent(inout) :: tsOutputType
    integer :: defaultTsOutputType
    character*200 comlin
    integer :: narg, length, status, ndot, nlast, na, nc
    integer :: nprocs
    character*200 nprocschar
    character*200 tsoutchar
    logical :: exists
!---------------------------------------------------------------------------------------------------------------
    
    ! Get the number of command-line arguments
    narg = command_argument_count()
    
    ! Initialize mpsimFile, mplogFile, and logType
    mpsimFile    = ""
    mplogFile    = ""
    logType      = 1
    parallel     = .false.
    nprocs       = 0
    tsOutputType = 0
    defaultTsOutputType = 1

    ! Loop through the command-line arguments (if any)
    na = 1
    do while (na <= narg)
        call get_command_argument(na, comlin, length, status)
        if ((na == 1) .and. (comlin(1:1) /= "-")) then
            na = na + 1
            mpsimFile = comlin(1:length)
            ! Check for existence of the file, adding .mpsim extension if necessary
            call CheckSimulationFile(mpsimFile)
        else
            na = na + 1
            select case (comlin(1:length))
            case ("-shortlog")
                ! -shortlog option
                if (logType.ne.1) then
                    call ustop('Conflicting or redundant log options on the command line. Stop.')
                else
                    logType = 2
                end if
            case ("-nolog")
                ! -nolog option
                if (logType.ne.1) then
                    call ustop('Conflicting or redundant log options on the command line. Stop.')
                else
                    logType = 0
                end if
            case ("-logname")
                ! -logname option
                if (mplogFile == "") then
                    call get_command_argument(na, comlin, length, status)
                    na = na + 1
                    if ((status /= 0) .or. (comlin(1:1) == "-")) then
                        call ustop('Invalid or missing log file name on the command line. Stop.')
                    else
                        mplogFile = comlin(1:length)
                    end if
                else
                    call ustop('Conflicting log file names on the command line. Stop.')
                end if
            case ("-parallel")
                ! -parallel option
                parallel = .true.
            case ("-np")
                ! -np option
                call get_command_argument(na, comlin, length, status)
                na = na + 1
                if ((status /= 0) .or. (comlin(1:1) == "-")) then
                    call ustop('Invalid or missing number of processes from command line. Stop.')
                else
                    nprocschar = comlin(1:length)
                    read(nprocschar,*) nprocs
                    if( nprocs .gt. 1 ) parallel = .true.
                end if
            case ("-tsoutput")
                ! -tsoutput option
                call get_command_argument(na, comlin, length, status)
                na = na + 1
                if ((status /= 0) .or. (comlin(1:1) == "-")) then
                    call ustop('Invalid or missing output for timeseries on commmand line. Stop.')
                else
                    tsoutchar = comlin(1:length)
                    read(tsoutchar,*) tsOutputType
                end if
            case default
                if (comlin(1:1).eq."-") then
                    call ustop('Unrecognized option on the command line. Stop.')
                else
                    call ustop('An error occurred processing the command line. Stop.')
                end if
            end select
        end if
    end do
    if ((mplogFile /= "") .and. (logType.eq.0))                                   &
        call ustop('Options -logname and -nolog both on the command line. Stop.')
    
    ! If log file name not specified, set to default
    if (mplogFile == "") mplogFile = "mpath7.log"

    ! Set parallel processes
    if ( parallel .and. (nprocs .eq. 0) ) then 
        ! If parallel and no np specified, set number of processors 
        call omp_set_num_threads( omp_get_num_procs() )
        ! Set default outputType to parallel consolidated
        if ( tsOutputType .eq. 0 ) tsOutputType = defaultTsOutputType
    else if ( parallel .and. (nprocs .gt. 1 ) ) then 
        ! If parallel and np specified, set np processes
        call omp_set_num_threads( nprocs )
        if ( tsOutputType .eq. 0 ) tsOutputType = defaultTsOutputType
    else if ( nprocs .eq. 1 ) then
        ! If nprocs equal one, then serial
        parallel = .false.
        call omp_set_num_threads( nprocs )
    else if ( omp_get_max_threads() .gt. 1 ) then 
        ! If max threads defined through OMP_NUM_THREADS, use it
        parallel = .true.
        nprocs = omp_get_max_threads()
        call omp_set_num_threads( nprocs )
        if ( tsOutputType .eq. 0 ) tsOutputType = defaultTsOutputType
    end if 


    return
   

    end subroutine ParseCommandLine

    subroutine PromptSimulationFile(mpsimFile)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
    use utl7module,only : urword
    implicit none
    character*(*),intent(inout) :: mpsimFile
    integer :: icol, istart, istop, n
    real(kind=4) :: r
    logical :: exists
!---------------------------------------------------------------------------------------------------------------

    ! Prompt user to enter simulation file name
    icol = 1
    write(*, *) 'Enter the MODPATH simulation file: '
    read(*, '(a)') mpsimFile
    call urword(mpsimFile,icol,istart,istop,0,n,r,0,0)
    mpsimFile = mpsimFile(istart:istop)
    
    ! Check for existence of simulation file, adding .mpsim extension if necessary
    call CheckSimulationFile(mpsimFile)

    return

    end subroutine PromptSimulationFile

    subroutine CheckSimulationFile(mpsimFile)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
    implicit none
    character*(*),intent(inout) :: mpsimFile
    integer :: nc
    logical :: exists
!---------------------------------------------------------------------------------------------------------------
    
    ! Check for existence of the file
    inquire (file=mpsimFile, exist=exists)
    if(.not. exists) then
        ! Add .mpsim extension, check again, and stop if not found
        nc = index(mpsimFile,' ')
        mpsimFile(nc:nc+5)='.mpsim'
        inquire (file=mpsimFile, exist=exists)
        if(.not. exists) call ustop('The specified simulation file could not be found. Stop.')
    end if
    mpsimFile = trim(mpsimFile)

    return
    
    end subroutine CheckSimulationFile
    
    subroutine ReadNameFile(filename, outUnit, gridFileType)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
    use UTL8MODULE,only : urword, ustop
    implicit none
    character*(*),intent(in) :: filename
    integer,intent(in) :: outUnit
    integer,intent(inout) :: gridFileType
    character(len=200) :: line
    character(len=150) :: fname
    character(len=16) :: filtyp
    character(len=30) :: gridFileTypeString
    character(len=80) :: message
    character(len=132) :: errMessage
    integer,dimension(6) :: nfiltyp
    integer :: inUnit, n, icol, ityp1, ityp2, inam1, inam2, nc, iflen, numflag, istart, istop
    doubleprecision :: r
    logical :: complete
!---------------------------------------------------------------------------------------------------------------
    
    errMessage = ' '
    
    do n = 1, 6
        nfiltyp(n) = 0
    end do
    gridFile = ' '
    tdisFile = ' '
    mpbasFile = ' '
    headFile = ' '
    budgetFile = ' '
    gridMetaFile = ' '
    
    inUnit = 99
    open(unit=inUnit, file=filename, status='old', form='formatted', access='sequential')
    
    write(outUnit, '(1x/a)') 'MODPATH name file data'
    write(outUnit, '(a)')    '----------------------'
    
        gridFileType = 0
    do
        read(inUnit, '(a)', end=1000) line
        ! Check for comment lines or blank lines and skip over them if present
        if(len_trim(line) .eq. 0) cycle
        if(line(1:1) .eq. '#') cycle
        if(line(1:1) .eq. '!') cycle
        if(line(1:2) .eq. '//') cycle
        
        ! Check for unit numbers
        numflag = 0
        read(line, *, err=200) filtyp, n
        numflag = 1
        goto 200
200     continue
        
        icol=1
        call urword(line, icol, istart, istop, 1, n, r, outUnit, inUnit)
        filtyp = line(istart:istop)
        if(numflag .eq. 1) then
            call urword(line, icol, istart, istop, 2, n, r, outUnit, inUnit)
        end if
        call urword(line, icol, istart, istop, 0, n, r, outUnit, inUnit)
        iflen = istop - istart + 1
        fname = line(istart:istop)
        
        if(filtyp .eq. 'DIS') then
            gridFile = fname(1:iflen)
            open(unit=disUnit,file=gridFile,status='old', form='formatted', access='sequential')
            write(outUnit,'(A15,A)') 'MODFLOW-2005/MODFLOW-USG Structured Grid File (DIS): ', gridFile(1:iflen)
            nfiltyp(1) = 1
            nfiltyp(2) = 1
            gridFileType = 1
        else if(filtyp .eq. 'DISU') then
            gridFile = fname(1:iflen)
            open(unit=disUnit,file=gridFile,status='old', form='formatted', access='sequential', err=500, iomsg=errMessage)
            write(outUnit,'(A15,A)') 'DISU File: ', gridFile(1:iflen)
            nfiltyp(1) = 2
            nfiltyp(2) = 1
            gridFileType = 2
        else if(filtyp .eq. 'GRBDIS') then
            gridFile = fname(1:iflen)
            open(unit=disUnit,file=gridFile,form='unformatted',access='stream',status='old',action='read', err=500, &
                iomsg=errMessage)
            write(outUnit,'(A15,A)') 'GRBDIS File: ', gridFile(1:iflen)
            nfiltyp(1) = 3
            gridFileType = 3
        else if(filtyp .eq. 'GRBDISV') then
            gridFile = fname(1:iflen)
            open(unit=disUnit,file=gridFile,form='unformatted',access='stream',status='old',action='read', err=500, &
                iomsg=errMessage)
            write(outUnit,'(A15,A)') 'GRBDISV File: ', gridFile(1:iflen)
            nfiltyp(1) = 4
            gridFileType = 4
        else if(filtyp .eq. 'GRBDISU') then
            call ustop('Binary grid file type DISU not yet supported. Stop.')
        else if(filtyp .eq. 'TDIS') then
            if(nfiltyp(1) .ge. 3) then
                tdisFile = fname(1:iflen)
                open(unit=tdisUnit,file=tdisFile,status='old', form='formatted', access='sequential', err=500, iomsg=errMessage)
                write(outUnit,'(A15,A)') 'TDIS File: ', tdisFile(1:iflen)
                nfiltyp(2) = 1
            end if
        else if(filtyp .eq. 'MPBAS') then
            mpbasFile = fname(1:iflen)
            open(unit=mpbasUnit,file=mpbasFile,status='old', form='formatted', access='sequential', err=500, iomsg=errMessage)
            write(outUnit,'(A15,A)') 'MPBAS File: ', mpbasFile(1:iflen)
            nfiltyp(3) = 1
        else if(filtyp .eq. 'HEAD') then
            headFile = fname(1:iflen)
            write(outUnit,'(A15,A)') 'HEAD File: ', headFile(1:iflen)
            nfiltyp(4) = 1
        else if(filtyp .eq. 'BUDGET') then
            budgetFile = fname(1:iflen)
            write(outUnit,'(A15,A)') 'BUDGET File: ', budgetFile(1:iflen)
            nfiltyp(5) = 1
        else if(filtyp .eq. 'GRIDMETA') then
            gridMetaFile = fname(1:iflen)
            open(unit=gridMetaUnit,file=gridMetaFile,status='old', form='formatted', access='sequential', err=500, iomsg=errMessage)
            write(outUnit,'(A15,A)') 'GRIDMETA File: ', gridMetaFile(1:iflen)
            nfiltyp(6) = 1
        end if
          
        cycle
        
    end do
    
1000 continue
     
     if(gridFileType .eq. 0) then
        message = 'No valid grid file type was specified in the name file. Stop.'
        call ustop(message)
     else
         complete = .true.
         do n = 1, 5
             if(nfiltyp(n) .eq. 0) then
                 complete = .false.
                 message = 'The MODPATH name file is not complete. Stop.'
             else
                 if((nfiltyp(n) .eq. 2) .and. (nfiltyp(6) .eq. 0)) then
                     complete = .false.
                     message = 'A GRIDMETA file must be specified for grid type DISU. Stop.'
                 else if((nfiltyp(n) .eq. 5) .and. (nfiltyp(6) .eq. 0)) then
                     complete = .false.
                     message = 'A GRIDMETA file must be specified for grid type GRBDISU. Stop.'
                 end if
             end if
         end do
         if(.not. complete) then
             write(outUnit, '(1x,a)') message
             call ustop(message)
         end if
     end if
    
100 continue
    return
500 continue
    call ustop(errMessage)
    
    end subroutine
    
    subroutine WriteBudgetFileInfo(outUnit, budgetReader)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
    use BudgetReaderModule,only : BudgetReaderType
    use BudgetRecordHeaderModule,only : BudgetRecordHeaderType
    integer,intent(in) :: outUnit
    type(BudgetReaderType),intent(in) :: budgetReader
    
    write(outUnit,*)
    write(outUnit, '(1x,a)') 'Budget File Data'
    write(outUnit, '(1x,a)') '----------------'
    if(budgetReader%GetBudgetType() .eq. 1) then
        write(outUnit, '(1x,a20,2x,a)') 'Budget file type:',  'Structured grid'
    else if(budgetReader%GetBudgetType() .eq. 2) then
        write(outUnit, '(1x,a20,2x,a)') 'Budget file type:',  'Unstructured grid' 
    else
        write(outUnit, '(1x,a20,2x,a)') 'Budget file type:',  'Undetermined type'
        return
    end if
    if(budgetReader%GetBudgetFileFormat() .eq. 1) then
        write(outUnit, '(1x,a20,2x,a)') 'Budget file format:',  'Standard'
    else if(budgetReader%GetBudgetFileFormat() .eq. 2) then
        write(outUnit, '(1x,a20,2x,a)') 'Budget file format:',  'Compact'           
    end if
    if(budgetReader%GetPrecisionType() .eq. 1) then
            write(outUnit, '(1x,a20,2x,a)') 'Budget precision:',  'Single'
    else if(budgetReader%GetPrecisionType() .eq. 2) then
            write(outUnit, '(1x,a20,2x,a)') 'Budget precision:',  'Double'            
    end if
    
    return
    
    end subroutine WriteBudgetFileInfo
    
    subroutine WriteBudgetRecordHeaders(outUnit, budgetReader)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
    use BudgetReaderModule,only : BudgetReaderType
    use BudgetRecordHeaderModule,only : BudgetRecordHeaderType
    integer,intent(in) :: outUnit
    type(BudgetReaderType),intent(in) :: budgetReader
    integer :: n, m, auxCount, recordHeaderCount
    type(BudgetRecordHeaderType) :: budgetRecordHeader
    
    write(outUnit,*)
    write(outUnit, '(1x,a)') 'Budget Record Headers:'
    if(budgetReader%GetBudgetFileFormat() .eq. 1) then
        write(outUnit, '(1x,a)') '    Record    Period      Step      Text label'
        recordHeaderCount = budgetReader%GetRecordHeaderCount()
        do n = 1, recordHeaderCount
            budgetRecordHeader = budgetReader%GetRecordHeader(n)
            write(outUnit, '(1x, 3i10,2x,a16)')                                 &
              n, budgetRecordHeader%StressPeriod, budgetRecordHeader%TimeStep,  &
              budgetRecordHeader%TextLabel
        end do
    else if(budgetReader%GetBudgetFileFormat() .eq. 2) then
        write(outUnit, '(1x,a)')                                                &
          '    Record    Period      Step      Text label      Method         Step length       Period length          Total time'
        recordHeaderCount = budgetReader%GetRecordHeaderCount()
        do n = 1, recordHeaderCount
            budgetRecordHeader = budgetReader%GetRecordHeader(n)
            write(outUnit, '(1x, 3i10,2x,a16,i10,3E20.12)')                     &
              n, budgetRecordHeader%StressPeriod, budgetRecordHeader%TimeStep,  &
              budgetRecordHeader%TextLabel, budgetRecordHeader%Method,          &
              budgetRecordHeader%TimeStepLength,                                &
              budgetRecordHeader%StressPeriodLength,                            &
              budgetRecordHeader%TotalTime
            if(budgetRecordHeader%Method .eq. 5) then
                auxCount = budgetRecordHeader%GetAuxiliaryNamesCount()
                write(outUnit, '(58x,a,i10,5x,a,i5)')                           &
                  'List item count = ', budgetRecordHeader%ListItemCount,       &
                  'Auxiliary item count = ', auxCount
                do m = 1, auxCount
                    write(outUnit, '(58x,a)')                                   &
                      budgetRecordHeader%AuxiliaryNames(m)
                end do
            else if(budgetRecordHeader%Method .eq. 2) then
                write(outUnit, '(58x,a,i10)')                                   &
                  'List item count = ', budgetRecordHeader%ListItemCount                    
            end if
        end do
    end if
    
    return
        
    end subroutine WriteBudgetRecordHeaders
    
    subroutine WriteWaterBalanceSummary(outUnit, trackingEngine, cellData)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
    use ParticleTrackingEngineModule,only : ParticleTrackingEngineType
    use ModpathCellDataModule,only : ModpathCellDataType
    integer,intent(in) :: outUnit
    type(ParticleTrackingEngineType),intent(in) :: trackingEngine
    type(ModpathCellDataType),intent(inout) :: cellData
    integer :: n, budgetIntervalBreakCount, maxErrorCell
    doubleprecision :: maxError
    doubleprecision,dimension(6) :: budgetIntervalBreaks
    integer,dimension(7) :: budgetIntervalBins
    
        budgetIntervalBreakCount = 6
        budgetIntervalBreaks(1) = 1.0d-02
        budgetIntervalBreaks(2) = 1.0d-01
        budgetIntervalBreaks(3) = 1.0d0
        budgetIntervalBreaks(4) = 1.0d+01
        budgetIntervalBreaks(5) = 5.0d+01
        budgetIntervalBreaks(6) = 1.0d+02
        maxErrorCell = 0
        maxError = 0.0d0
        do n = 1, budgetIntervalBreakCount  + 1
            budgetIntervalBins(n) = 0
        end do 
        
        call trackingEngine%GetVolumetricBalanceSummary(                        &
          budgetIntervalBreakCount, budgetIntervalBreaks, budgetIntervalBins,   &
          maxError, maxErrorCell)
        
        write(outUnit, *)
        write(outUnit, '(1X,A)') 'Volumetric water balance summary:'
        write(outUnit, *)
        
        write(outUnit, '(1X,I10,A,F8.2,A)')                                     &
          budgetIntervalBins(1), ' cells had errors less than or equal to',     &
          budgetIntervalBreaks(1), ' percent'
        do n = 2, budgetIntervalBreakCount
            write(outUnit, '(1X,I10,A,F8.2,A,F8.2,A)')                          &
              budgetIntervalBins(n), ' cells had errors between ',              &
              budgetIntervalBreaks(n-1), ' and ',                               &
              budgetIntervalBreaks(n), ' percent'                
        end do
        write(outUnit, '(1X,I10,A,F8.2,A)')                                     &
          budgetIntervalBins(budgetIntervalBreakCount + 1),                     &
          ' cells had errors greater than ',                                    &
          budgetIntervalBreaks(budgetIntervalBreakCount), ' percent'
        
        write(outUnit, *)
        write(outUnit, '(1X,A,E12.5,A,I10)')                                    &
          'A maximum error of ', maxError, ' percent occurred in cell ',        &
          maxErrorCell
        write(outUnit, *)
        call trackingEngine%FillCellBuffer(maxErrorCell, cellData)
        call trackingEngine%WriteCellBuffer(outUnit, cellData,                  &
          simulationData%TrackingOptions%BackwardTracking)
    
        return
        
    end subroutine WriteWaterBalanceSummary
    
    subroutine WriteParticleSummaryInfo(simulationData, outUnit)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
    use ParticleModule,only : ParticleType
    implicit none
    type(ModpathSimulationDataType),target,intent(in) :: simulationData
    type(ParticleType),pointer :: p
    integer,intent(in) :: outUnit
    integer :: groupIndex, particleIndex, n
    integer,dimension(0:9) :: statusBins
    
    do n = 0, 9
        statusBins(n) = 0
    end do
    
    do groupIndex = 1, simulationData%ParticleGroupCount
        do particleIndex = 1, simulationData%ParticleGroups(groupIndex)%TotalParticleCount
            p => simulationData%ParticleGroups(groupIndex)%Particles(particleIndex)
            if(p%Status.ge.0 .and. p%Status.le.9) then
                statusBins(p%Status) = statusBins(p%Status) + 1
            else
                statusBins(9) = statusBins(9) + 1
            end if
        end do
    end do
    
    ! Write to listing file
    write(outUnit, '(1x/a)') 'Particle Summary:'
    write(outUnit, '(i10,1x,a)') statusBins(0), 'particles are pending release.'
    write(outUnit, '(i10,1x,a)') statusBins(1), 'particles remain active.'
    write(outUnit, '(i10,1x,a)') statusBins(2), 'particles terminated at boundary faces.'
    write(outUnit, '(i10,1x,a)') statusBins(4), 'particles terminated at weak source cells.'
    write(outUnit, '(i10,1x,a)') statusBins(5),                                 &
      'particles terminated at strong source/sink cells or other cells with no potential exit face.'
    write(outUnit, '(i10,1x,a)') statusBins(6), 'particles terminated in cells with a specified zone number.'    
    write(outUnit, '(i10,1x,a)') statusBins(7), 'particles were stranded in inactive or dry cells.'
    write(outUnit, '(i10,1x,a)') statusBins(8), 'particles were unreleased.'
    write(outUnit, '(i10,1x,a)') statusBins(9), 'particles have an unknown status.'
    write(outUnit, '(a)') ' '

    ! Write to screen
    write(*, '(1x/a)') 'Particle Summary:'
    write(*, '(i10,1x,a)') statusBins(0), 'particles are pending release.'
    write(*, '(i10,1x,a)') statusBins(1), 'particles remain active.'
    write(*, '(i10,1x,a)') statusBins(2), 'particles terminated at boundary faces.'
    write(*, '(i10,1x,a)') statusBins(3), 'particles terminated at weak sink cells.'
    write(*, '(i10,1x,a)') statusBins(4), 'particles terminated at weak source cells.'
    write(*, '(i10,1x,a)') statusBins(5), 'particles terminated at strong source/sink cells.'
    write(*, '(i10,1x,a)') statusBins(6), 'particles terminated in cells with a specified zone number.'    
    write(*, '(i10,1x,a)') statusBins(7), 'particles were stranded in inactive or dry cells.'
    write(*, '(i10,1x,a)') statusBins(8), 'particles were unreleased.'
    write(*, '(i10,1x,a)') statusBins(9), 'particles have an unknown status.'
    write(*, '(a)') ' '
    
    end subroutine WriteParticleSummaryInfo


    ! OBSERVATION CELLS
    subroutine WriteObservationHeader(outUnit, cellNumber, referenceTime,   &
                                             originX, originY, rotationAngle)
        !---------------------------------------------------------------------
        ! Doc me
        !---------------------------------------------------------------------
        ! Specifications
        !---------------------------------------------------------------------
        implicit none
        integer,intent(in) :: outUnit, cellNumber
        doubleprecision,intent(in) :: referenceTime, originX, originY, rotationAngle
        integer :: version, subversion
        !----------------------------------- 
        version = 7
        subversion = 2
        write(outUnit, '(a,2i10)') 'MODPATH_CELL_OBSERVATION_FILE', version, subversion
        write(outUnit, '(i8,1x,4e18.10)') cellNumber, referenceTime, originX, originY, &
        rotationAngle
        write(outUnit, '(a)') 'END HEADER'

    end subroutine WriteObservationHeader



    end program MPath7
