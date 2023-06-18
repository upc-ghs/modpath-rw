module ModpathSimulationDataModule
  use ParticleTrackingOptionsModule,only : ParticleTrackingOptionsType
  use ParticleGroupModule,only : ParticleGroupType
  use ModflowRectangularGridModule,only : ModflowRectangularGridType
  use StartingLocationReaderModule,only : ReadAndPrepareLocations, &
                                               MassParticlesArray, &
                               CreateMassParticlesAsInternalArray, &
                    CreateMassParticlesAsQuasiRandomInternalArray, &
                         CreateMassParticlesAsRandomInternalArray, &
                                       CreateMassParticlesOnFaces, &
                                pr_CreateParticlesAsInternalArray
  use TimeDiscretizationDataModule,only : TimeDiscretizationDataType
  use FlowModelDataModule,only : FlowModelDataType
  implicit none
  
! Set default access status to private
  private
  
! Public derived data type definitions
!--------------------------------------
! type: 
!--------------------------------------
  type,public :: ModpathSimulationDataType
    character(len=200) :: NameFile
    character(len=200) :: ListingFile
    integer :: TraceMode, TraceGroup, TraceID
    integer :: SimulationType
    integer :: TrackingDirection
    integer :: WeakSinkOption
    integer :: WeakSourceOption
    integer :: ReferenceTimeOption
    integer :: StoppingTimeOption
    integer :: BudgetOutputOption
    integer :: TimeseriesOutputOption
    integer :: PathlineFormatOption
    integer :: ZoneDataOption
    integer :: RetardationFactorOption
    integer :: AdvectiveObservationsOption
    integer :: TimePointOption
    integer :: ParticleGroupCount
    integer :: TotalParticleCount
    integer :: TimePointCount
    integer :: StopZone
    integer :: BudgetCellsCount
    doubleprecision :: StopTime
    doubleprecision :: ReferenceTime
    doubleprecision :: TimePointInterval
    character(len=200) :: EndpointFile
    character(len=200) :: PathlineFile
    character(len=200) :: TimeseriesFile
    character(len=200) :: TraceFile
    character(len=200) :: AdvectiveObservationsFile
    character(len=200) :: DispersionFile                   ! RWPT
    integer            :: ParticlesMassOption              ! RWPT
    integer            :: SolutesOption                    ! RWPT
    integer            :: EndpointOutputOption             ! RWPT
    logical            :: shouldUpdateDispersion = .false. ! RWPT
    logical            :: isMF6 = .false.                  ! RWPT
    integer,dimension(:),allocatable :: BudgetCells
    integer,dimension(:),allocatable :: Zones
    doubleprecision,dimension(:),allocatable :: Retardation
    doubleprecision,dimension(:),allocatable :: TimePoints
    type(ParticleGroupType),dimension(:),allocatable :: ParticleGroups
    type(ParticleTrackingOptionsType),allocatable :: TrackingOptions
    logical :: isUniformPorosity =.false.                  ! RWPT
    logical :: isUniformRetardation = .false.              ! RWPT
    doubleprecision :: uniformPorosity = 1d0               ! RWPT
    doubleprecision :: uniformRetardation = 1d0            ! RWPT
    class(TimeDiscretizationDataType), pointer :: tdisData ! RWPT
    integer :: currentSeqNumber                            ! RWPT
  contains
    procedure :: ReadFileHeaders=>pr_ReadFileHeaders
    procedure :: ReadData=>pr_ReadData
    procedure :: ReadGPKDEData=>pr_ReadGPKDEData   ! RWPT
    procedure :: ReadOBSData=>pr_ReadOBSData       ! RWPT
    procedure :: ReadRWOPTSData=>pr_ReadRWOPTSData ! RWPT
    procedure :: ReadICData=>pr_ReadICData         ! RWPT
    procedure :: ReadSRCData=>pr_ReadSRCData       ! RWPT
    procedure :: SetUniformPorosity=>pr_SetUniformPorosity ! RWPT
  end type


contains


  subroutine pr_ReadFileHeaders(this, inUnit)
    use UTL8MODULE,only : u8rdcom
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    class(ModpathSimulationDataType) :: this
    integer,intent(in) :: inUnit
    integer :: outUnit, errorCode
    character(len=200) line
    !--------------------------------------------------------------
  
    outUnit = 0
    call u8rdcom(inUnit, outUnit, line, errorCode)
    
    ! Assign the name file
    this%NameFile = line
    
    ! Read MODPATH listing file filename
    read(inUnit, '(a)') this%ListingFile
  
  end subroutine pr_ReadFileHeaders


  ! Inform simulationData about uniform porosity
  subroutine pr_SetUniformPorosity(this, basicData)
    use ModpathBasicDataModule,only : ModpathBasicDataType
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    class(ModpathSimulationDataType) :: this
    type(ModpathBasicDataType),intent(in) :: basicData
    !--------------------------------------------------------------
 
    this%isUniformPorosity = basicData%isUniformPorosity
    if ( this%isUniformPorosity ) then 
      ! Needs something to handle the case were ibound(1) != 0 ?
      this%uniformPorosity = basicData%Porosity(1)
    end if

    ! At trackingoptions, goes to tracksubcell
    this%TrackingOptions%isUniformPorosity = basicData%isUniformPorosity

  end subroutine pr_SetUniformPorosity


  ! Read simulation data 
  subroutine pr_ReadData(this, inUnit, outUnit, ibound, timeDiscretization, grid)
    use UTL8MODULE,only : urword, ustop, u1dint, u1drel, u1ddbl, u8rdcom, &
                          u3ddblmpusg, u3dintmp, u3dintmpusg, u3ddblmp, ugetnode
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    class(ModpathSimulationDataType), target :: this
    class(ModflowRectangularGridType),intent(in) :: grid
    integer,intent(in) :: inUnit, outUnit
    integer,dimension(:),allocatable :: cellsPerLayer
    integer,dimension(grid%CellCount),intent(in) :: ibound
    class(TimeDiscretizationDataType),intent(in),target :: timeDiscretization
    integer :: icol, istart, istop, n, kper, kstp, seqNumber, particleCount, nn, slocUnit, errorCode
    integer :: releaseOption, releaseTimeCount
    doubleprecision :: initialReleaseTime, releaseInterval
    doubleprecision,dimension(:),allocatable :: releaseTimes
    doubleprecision :: frac, r
    character(len=24) aname(2)
    character(len=200) line
    DATA aname(1) /'              ZONE ARRAY'/
    DATA aname(2) /'                 RFACTOR'/
    !---------------------------------------------

    ! Deallocate arrays
    if(allocated(this%Zones)) deallocate(this%Zones)
    if(allocated(this%Retardation)) deallocate(this%Retardation)
    if(allocated(this%TimePoints)) deallocate(this%TimePoints)
    if(allocated(this%ParticleGroups)) deallocate(this%ParticleGroups)
    if(allocated(this%TrackingOptions)) deallocate(this%TrackingOptions)
    allocate(this%Zones(grid%CellCount))
    allocate(this%Retardation(grid%CellCount))
    allocate(cellsPerLayer(grid%LayerCount))
    do n = 1, grid%LayerCount
        cellsPerLayer(n) = grid%GetLayerCellCount(n)
    end do
    ! Allocate TrackingOptions
    allocate(this%TrackingOptions)
    
    ! Write header to the listing file
    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW simulation file data'
    write(outUnit, '(1x,a)') '-------------------------------'
    
    ! Rewind simulation file, then re-read comment lines and the first two non-comment
    ! lines containing the name file and listing file names that were read previously.
    rewind(inUnit)
    call u8rdcom(inUnit, outUnit, line, errorCode)
    read(inUnit, '(a)') line
    
    ! Read simulation options line, then parse line using subroutine urword
    read(inUnit, '(a)') line
    
    ! Simulation type
    icol = 1
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    this%SimulationType = n
    
    ! Tracking direction
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    this%TrackingDirection = n
    
    ! Weak sink option
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    this%WeakSinkOption = n
    
    ! Weak source option
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    this%WeakSourceOption = n
    
    ! Budget output option
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    this%BudgetOutputOption = n
    
    ! Trace mode
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    this%TraceMode = n

    ! Timeseries output option
    ! 0: Original behavior, timeseries records for active particles
    ! 1: Timeseries records for all particles
    ! 2: No timeseries records for any particle ! RWPT
    call urword(line, icol, istart, istop, 2, n, r, -1, 0)
    ! If error while reading the last option (could be triggered by # comments ) 
    if ( line(len(line):len(line)).eq.'E' ) then
      ! Continue as zero
      this%TimeseriesOutputOption = 0
    else
      ! Read from input
      if (istart.eq.len(line)) then
        this%TimeseriesOutputOption = 0
      else
        this%TimeseriesOutputOption = n
      end if
    end if

    ! Endpoint output option
    call urword(line, icol, istart, istop, 2, n, r, -1, 0)
    ! If error while reading the last option (could be triggered by # comments ) 
    if ( line(len(line):len(line)).eq.'E' ) then
      ! Continue as zero
      this%EndpointOutputOption = 0
    else
      ! Read from input
      if (istart.eq.len(line)) then
        this%EndpointOutputOption = 0
      else
        this%EndpointOutputOption = n
      end if
    end if

    ! Particles mass option
    call urword(line, icol, istart, istop, 2, n, r, -1, 0)
    ! If error while reading the last option (could be triggered by # comments ) 
    if ( line(len(line):len(line)).eq.'E' ) then
      ! Continue as zero
      this%ParticlesMassOption = 0
    else
      ! Read from input
      if (istart.eq.len(line)) then
        this%ParticlesMassOption = 0
      else
        this%ParticlesMassOption = n
      end if
    end if

    ! Solutes option
    call urword(line, icol, istart, istop, 2, n, r, -1, 0)
    ! If error while reading the last option (could be triggered by # comments ) 
    if ( line(len(line):len(line)).eq.'E' ) then
      ! Continue as zero
      this%SolutesOption = 0
    else
      ! Read from input
      if (istart.eq.len(line)) then
        this%SolutesOption = 0
      else
        this%SolutesOption = n
      end if
    end if


    ! Pathline format option (hardwire value 1 = consolidate)
    this%PathlineFormatOption = 1
    
    ! Advective observations option (hardwire value 1 = do not use advective observations)
    this%AdvectiveObservationsOption = 1
    
    ! Read coordinate output file names based on simulation type
    select case (this%SimulationType)
      case (1)
        write(outUnit,'(A,I2,A)') 'Endpoint Analysis (Simulation type =',this%SimulationType,')'
        read(inUnit,'(a)') this%EndpointFile
        icol=1
        call urword(this%EndpointFile, icol, istart, istop, 0, n, r, 0, 0)
        this%Endpointfile=this%EndpointFile(istart:istop)
      case (2)
        write(outUnit,'(A,I2,A)') 'Pathline Analysis (Simulation type =', this%SimulationType, ')'
        read(inUnit, '(a)') this%EndpointFile
        icol = 1
        call urword(this%EndpointFile,icol,istart,istop,0,n,r,0,0)
        this%EndpointFile = this%EndpointFile(istart:istop)
        read(inUnit, '(a)') this%PathlineFile
        icol = 1
        call urword(this%PathlineFile, icol, istart, istop, 0, n, r, 0, 0)
        this%PathlineFile = this%PathlineFile(istart:istop)
      case (3)
        write(outUnit,'(A,I2,A)') 'Timeseries Analysis (Simulation type =',this%SimulationType,')'
        read(inUnit, '(a)') this%EndpointFile
        icol=1
        call urword(this%EndpointFile, icol, istart, istop, 0, n, r, 0, 0)
        this%Endpointfile=this%EndpointFile(istart:istop)
        read(inUnit, '(a)') this%TimeseriesFile
        icol = 1
        call urword(this%TimeseriesFile, icol, istart, istop, 0, n, r, 0, 0)
        this%TimeseriesFile = this%TimeseriesFile(istart:istop)
        if(this%AdvectiveObservationsOption.eq.2) then
          read(inUnit, '(a)') this%AdvectiveObservationsFile
          icol = 1
          call urword(this%AdvectiveObservationsFile, icol, istart, istop, 0, n, r,0,0)
          this%AdvectiveObservationsFile = this%AdvectiveObservationsFile(istart:istop)
        end if
      case (4)
        write(outUnit,'(A,I2,A)') 'Combined Pathline and Timeseries Analysis (Simulation type =', this%SimulationType, ')'
        read(inUnit, '(a)') this%EndpointFile
        icol = 1
        call urword(this%EndpointFile,icol,istart,istop,0,n,r,0,0)
        this%EndpointFile = this%EndpointFile(istart:istop)
        read(inUnit, '(a)') this%PathlineFile
        icol = 1
        call urword(this%PathlineFile, icol, istart, istop, 0, n, r, 0, 0)
        this%PathlineFile = this%PathlineFile(istart:istop)
        read(inUnit, '(a)') this%TimeseriesFile
        icol = 1
        call urword(this%TimeseriesFile, icol, istart, istop, 0, n, r, 0, 0)
        this%TimeseriesFile = this%TimeseriesFile(istart:istop)
        if(this%AdvectiveObservationsOption.eq.2) then
          read(inUnit, '(a)') this%AdvectiveObservationsFile
          icol = 1
          call urword(this%AdvectiveObservationsFile, icol, istart, istop, 0, n, r,0,0)
          this%AdvectiveObservationsFile = this%AdvectiveObservationsFile(istart:istop)
        end if
      ! RWPT
      case(5)
        write(outUnit,'(A,I2,A)') 'RWPT with Timeseries Analysis (Simulation type =', this%SimulationType, ')'
        read(inUnit, '(a)') this%EndpointFile
        icol = 1
        call urword(this%EndpointFile,icol,istart,istop,0,n,r,0,0)
        this%EndpointFile = this%EndpointFile(istart:istop)
        read(inUnit, '(a)') this%TimeseriesFile
        icol = 1
        call urword(this%TimeseriesFile, icol, istart, istop, 0, n, r, 0, 0)
        this%TimeseriesFile = this%TimeseriesFile(istart:istop)
        if(this%AdvectiveObservationsOption.eq.2) then
          read(inUnit, '(a)') this%AdvectiveObservationsFile
          icol = 1
          call urword(this%AdvectiveObservationsFile, icol, istart, istop, 0, n, r,0,0)
          this%AdvectiveObservationsFile = this%AdvectiveObservationsFile(istart:istop)
        end if
        this%TrackingOptions%RandomWalkParticleTracking = .true.
      case(6)
        write(outUnit,'(A,I2,A)') 'RWPT with Pathline and Timeseries Analysis (Simulation type =', this%SimulationType, ')'
        read(inUnit, '(a)') this%EndpointFile
        icol = 1
        call urword(this%EndpointFile,icol,istart,istop,0,n,r,0,0)
        this%EndpointFile = this%EndpointFile(istart:istop)
        read(inUnit, '(a)') this%PathlineFile
        icol = 1
        call urword(this%PathlineFile, icol, istart, istop, 0, n, r, 0, 0)
        this%PathlineFile = this%PathlineFile(istart:istop)
        read(inUnit, '(a)') this%TimeseriesFile
        icol = 1
        call urword(this%TimeseriesFile, icol, istart, istop, 0, n, r, 0, 0)
        this%TimeseriesFile = this%TimeseriesFile(istart:istop)
        if(this%AdvectiveObservationsOption.eq.2) then
          read(inUnit, '(a)') this%AdvectiveObservationsFile
          icol = 1
          call urword(this%AdvectiveObservationsFile, icol, istart, istop, 0, n, r,0,0)
          this%AdvectiveObservationsFile = this%AdvectiveObservationsFile(istart:istop)
        end if
        this%TrackingOptions%RandomWalkParticleTracking = .true.
      case(7)
        write(outUnit,'(A,I2,A)') 'RWPT Endpoint Analysis (Simulation type =', this%SimulationType, ')'
        read(inUnit, '(a)') this%EndpointFile
        icol = 1
        call urword(this%EndpointFile,icol,istart,istop,0,n,r,0,0)
        this%EndpointFile = this%EndpointFile(istart:istop)
        this%TrackingOptions%RandomWalkParticleTracking = .true.
      case default
        call ustop('Invalid simulation type. Stop.')
    end select
    
    ! Read trace mode filename if trace mode is on
    if(this%TraceMode .gt. 0) then
      read(inUnit,'(a)') this%TraceFile
      icol=1
      call urword(this%EndpointFile, icol, istart, istop, 0, n, r, 0, 0)
      this%TraceFile=this%TraceFile(istart:istop)
      read(inUnit,*) this%TraceGroup, this%TraceID
    end if
    
    ! Read budget cells
    read(inUnit, *) this%BudgetCellsCount
    if(allocated(this%BudgetCells)) then
        deallocate(this%BudgetCells)
    end if
    allocate(this%BudgetCells(this%BudgetCellsCount))
    if(this%BudgetCellsCount .gt. 0) then
      read(inUnit, *) (this%BudgetCells(n), n = 1, this%BudgetCellsCount)
    end if
 
    ! RWPT 
    ! Only allow forward tracking for RWPT simulations
    !if ((( this%SimulationType .eq. 5 ) .or.  &
    !     ( this%SimulationType .eq. 6 ) .or.  & 
    !     ( this%SimulationType .eq. 7 )    )  &
    !  .and. ( this%TrackingDirection .eq. 2 ) ) then 
    !  call ustop('Random Walk Particle Tracking only accepts Forward tracking. Stop.')
    !end if

    ! Tracking direction
    select case(this%TrackingDirection)
      case(1)
        write(outUnit,'(A,I2,A)') 'Forward tracking (Tracking direction = ', this%TrackingDirection,')'
      case(2)
        write(outUnit,'(A,I2,A)') 'Backward tracking (Tracking direction =', this%TrackingDirection,')'
      case default
        call ustop('Invalid tracking direction code. Stop.')
    end select
    
    ! Weak sink option
    select case(this%WeakSinkOption)
      case (1)
        write(outUnit, '(A)') 'Let particles pass through weak sink cells (Weak sink option = 1)'
      case (2)
        write(outUnit,'(A)') 'Stop particles when they enter weak sink cells. (Weak sink option = 2)'
      case default
        call ustop('Invalid weak sink option.')
    end select

    ! Weak source option   
    select case(this%WeakSourceOption)
      case(1)
      write(outUnit,'(A)') 'Let particles pass through weak source cells for backtracking simulations (Weak source option = 1)'
      case(2)
      write(outUnit,'(A)') 'Stop particles when they enter weak source cells for backtracking simulations (Weak source option = 2)'
      case default
      call ustop('Invalid weak source option.')
    end select

    ! Timeseries output option
    select case(this%TimeseriesOutputOption)
      case (0)
        write(outUnit, '(A)') 'Timeseries output for active particles only (Timeseries output option = 0)'
      case (1)
        write(outUnit,'(A)') 'Timeseries output for all particles (Timeseries output option = 1)'
      case (2)
        write(outUnit,'(A)') 'No timeseries output, skip TimeseriesWriter (Timeseries output option = 2)'
      case default
        call ustop('Invalid timeseries output option.')
    end select

    ! EndpointOutputOption
    select case(this%EndpointOutputOption)
      case (0)
        write(outUnit,'(A)') 'Endpoint file will be written with MODPATH-v7 format (EndpointOutputOption= 0)'
      case (1)
        write(outUnit,'(A)') 'Endpoint file will be written with MODPATH-RW format (EndpointOutputOption= 1)'
      case (2)
        write(outUnit,'(A)') 'Endpoint file will not be written (EndpointOutputOption= 2)'
      case default
        call ustop('Invalid endpoint output option.')
    end select

    ! Particles mass option
    select case(this%ParticlesMassOption)
      case (0)
        write(outUnit, '(A)') 'Particle groups with default mass and soluteid (ParticlesMassOption = 0)'
      case (1)
        write(outUnit,'(A)') 'Particle groups will read mass (ParticlesMassOption = 1)'
      case (2)
        write(outUnit,'(A)') 'Particle groups will read both mass and soluteid (ParticlesMassOption = 2)'
      case default
        call ustop('Invalid particles mass option.')
    end select
 
    ! Solutes dispersion option
    select case(this%SolutesOption)
      case (0)
        write(outUnit, '(A)') 'Solutes with the same dispersion (SolutesOption= 0)'
      case (1)
        write(outUnit,'(A)') 'Solutes with specific dispersion (SolutesOption = 1)'
      case default
        call ustop('Invalid solutes option.')
    end select

    ! Reference time option
    read(inUnit, '(a)') line
    icol = 1
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    this%ReferenceTimeOption = n
  
    select case(this%ReferenceTimeOption)
      case(1)
        read(inUnit,*) this%ReferenceTime
        write(outUnit,'(A,E15.7)') 'Reference time = ', this%ReferenceTime
      case(2)
        read(inUnit, *) kper, kstp, frac
        this%ReferenceTime = timeDiscretization%GetTimeFromPeriodAndStep(kper, kstp, frac)
        write(outUnit,'(A,I6,A,I6)') 'Reference time will be calculated for: Period ', KPER,' Step ', KSTP
        write(outUnit,'(A,F10.7)') 'The relative time position within the time step is =',FRAC
        write(outUnit,'(A,E15.7)') 'Computed reference time = ', this%ReferenceTime
      case default
        call ustop('Invalid reference time option.')
    end select
    ! Regardless of time specifications
    ! save the pointer to timeDiscretization
    this%tdisData => timeDiscretization


    ! Read stopping option
    this%StopTime = 1.0E+30
    read(inUnit, '(a)') line
    icol = 1
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    this%StoppingTimeOption = n
    select case(this%StoppingTimeOption)
      case(1)
        write(outUnit,'(A,I2,A)')                                           &
          'Stop tracking at the beginning or end of the MODFLOW simulation (Stopping time option = ',  &
          this%StoppingTimeOption,')'
      case(2)
        write(outUnit,'(A,I2,A)')                                           &
          'Extend initial or final steady-state time step and continue tracking (Stopping time option = ', &
          this%StoppingTimeOption,')'
      case(3)
        write(outUnit,'(A,I2,A)')                                           &
          'Specify a limit for tracking time (Stoping time option = ',this%StoppingTimeOption,')'
        read(inUnit, *) this%StopTime
        write(outUnit,'(A,E15.7)') 'Stop time = ', this%StopTime
      case default
        call ustop('Invalid stop time code. Stop.')
    end select
  
    ! RWPT
    ! Time point data
    if((this%SimulationType .eq. 3) .or. (this%SimulationType .eq. 4) .or.  &
      (this%SimulationType .eq. 5) .or. (this%SimulationType .eq. 6)) then
      read(inUnit, *) this%TimePointOption
      if(this%TimePointOption .eq. 1) then
        read(inUnit, *) this%TimePointCount, this%TimePointInterval
        allocate(this%TimePoints(this%TimePointCount))
        if(this%TimePointCount .gt. 0) then
          this%TimePoints(1) = this%TimePointInterval
          do n = 2, this%TimePointCount
            this%TimePoints(n) = this%TimePoints(n-1) + this%TimePointInterval
          end do   
        end if
      else if(this%TimePointOption .eq. 2) then
        read(inUnit, *) this%TimePointCount
        allocate(this%TimePoints(this%TimePointCount))
        if(this%TimePointCount .gt. 0) then
          read(inUnit, *) (this%TimePoints(n), n = 1, this%TimePointCount)
        end if
      else
        ! write an error message and stop
      end if
    else
      this%TimePointOption = 0
      this%TimePointCount = 0
      this%TimePointInterval = 0.0d0
      allocate(this%TimePoints(0))      
    end if
  
    ! Zone array
    read(inUnit, '(a)') line
    icol = 1
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    this%ZoneDataOption = n
    if(this%ZoneDataOption .gt. 1) then
      write(outUnit, '(/a)') 'A zone array will be read.'
      read(inUnit,*) this%StopZone
      if(this%StopZone .lt. 1) then
        write(outUnit,'(A,I5)')                                               &
          'Particles will be allowed to pass through all zones. StopZone = ', this%StopZone
      else
        write(outUnit,'(A,I5)')                                               &
          'Particles will be terminated when they enter cells with a zone numbers equal to ', this%StopZone
      end if
      if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
        call u3dintmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
          grid%ColumnCount, grid%CellCount, this%Zones, ANAME(1))            
      else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
        call u3dintmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount, this%Zones,&
          ANAME(1), cellsPerLayer)
      else
        write(outUnit,*) 'Invalid grid type specified when reading zone array data.'
        write(outUnit,*) 'Stopping.'
        call ustop(' ')
      end if
    else
      write(outUnit,'(A)') 'The zone value for all cells = 1'
      this%StopZone = 0
      do n = 1, grid%CellCount
          this%Zones(n) = 1
      end do
    end if
      
    ! Retardation array
    read(inUnit, '(a)') line
    icol = 1
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    this%RetardationFactorOption = n  
    if(this%RetardationFactorOption .gt. 1) then
      write(outUnit,'(A)') 'The retardation factor array will be read.'
      if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
        call u3ddblmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,     &
          grid%ColumnCount, grid%CellCount, this%Retardation, aname(2)) 
      else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
        call u3ddblmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount, &
          this%Retardation, aname(2), cellsPerLayer)
      else
        write(outUnit,*) 'Invalid grid type specified when reading retardation array data.'
        write(outUnit,*) 'Stopping.'
        call ustop(' ')            
      end if
      ! RWPT
      ! Check if all cells have the same retardation factor
      if (all(this%Retardation.eq.this%Retardation(1))) then
          this%isUniformRetardation = .true.
          this%uniformRetardation = this%Retardation(1) 
      end if
    else
      write(outUnit,'(A)') 'The retardation factor for all cells = 1'
      do n = 1, grid%CellCount
        this%Retardation(n) = 1.0d0
      end do
      this%isUniformRetardation = .true.
    end if
      
    ! Particle data
    read(inUnit, *) this%ParticleGroupCount
    write(outUnit,'(A,I5)') 'Number of particle groups in simulation file = ', this%ParticleGroupCount
  
    seqNumber = 0
    this%currentSeqNumber = seqNumber ! Save seqNumber, see below
    this%TotalParticleCount = 0
    particleCount = 0
    if(this%ParticleGroupCount .gt. 0) then
      allocate(this%ParticleGroups(this%ParticleGroupCount))
      do n = 1, this%ParticleGroupCount
        this%ParticleGroups(n)%Group = n
        read(inUnit, '(a)') this%ParticleGroups(n)%Name
        read(inUnit, *) releaseOption
        
        select case (releaseOption)
          case (1)
            read(inUnit, *) initialReleaseTime
            call this%ParticleGroups(n)%SetReleaseOption1(initialReleaseTime)
          case (2)
            read(inUnit, *) releaseTimeCount, initialReleaseTime, releaseInterval
            call this%ParticleGroups(n)%SetReleaseOption2(initialReleaseTime, &
              releaseTimeCount, releaseInterval)
          case (3)
            read(inUnit, *) releaseTimeCount
            if(allocated(releaseTimes)) deallocate(releaseTimes)
            allocate(releaseTimes(releaseTimeCount))
            read(inUnit, *) (releaseTimes(nn), nn = 1, releaseTimeCount)
            call this%ParticleGroups(n)%SetReleaseOption3(releaseTimeCount,   &
              releaseTimes)
          case default
            ! write error message and stop
        end select
      
        read(inUnit, '(a)') line
        icol = 1
        call urword(line,icol,istart,istop,1,n,r,0,0)
        if(line(istart:istop) .eq. 'EXTERNAL') then
          call urword(line,icol,istart,istop,0,n,r,0,0)
          this%ParticleGroups(n)%LocationFile = line(istart:istop)
          slocUnit = 0
        else if(line(istart:istop) .eq. 'INTERNAL') then
          this%ParticleGroups(n)%LocationFile = ''
          slocUnit = inUnit
        else
          call ustop('Invalid starting locations file name. stop.')
        end if
        call ReadAndPrepareLocations(slocUnit, outUnit, this%ParticleGroups(n),   &
          ibound, grid%CellCount, grid, seqNumber)
        write(outUnit, '(a,i4,a,i10,a)') 'Particle group ', n, ' contains ',      &
          this%ParticleGroups(n)%TotalParticleCount, ' particles.'
        particleCount = particleCount + this%ParticleGroups(n)%TotalParticleCount

        ! RWPT
        if ( this%ParticlesMassOption .ge. 1 ) then 
          write(outUnit, '(a)') 'Will read ParticlesMass'
          ! Read group mass, is a proxy for concentrations
          ! when mass is uniform for a pgroup
          read(inUnit, '(a)') line
          icol = 1
          call urword(line,icol,istart,istop,3,n,r,0,0)
          if ( r.lt.0d0 ) then 
            write(outUnit,'(A)')'Given particles mass is negative, it should be positive. Stop.'
            call ustop('Given particles mass is negative, it should be positive. Stop.')
          end if
          this%ParticleGroups(n)%Mass = r
          this%ParticleGroups(n)%Particles(:)%Mass = this%ParticleGroups(n)%Mass
          ! Read the solute id for this group 
          if ( this%ParticlesMassOption .eq. 2 ) then 
            write(outUnit, '(a)') 'Will read a SpeciesID'
            read(inUnit, '(a)') line
            icol = 1
            call urword(line,icol,istart,istop,2,this%ParticleGroups(n)%Solute,r,0,0)
            if ( this%ParticleGroups(n)%Solute.lt.1 ) then 
              write(outUnit,'(A)')'Given SpeciesID is less than 1. Should be at least 1. Stop.'
              call ustop('Given SpeciesID is less than 1. Should be at least 1. Stop.')
            end if
          end if
        end if
      end do

      this%TotalParticleCount = particleCount
      write(outUnit, '(a,i10)') 'Total number of particles = ', this%TotalParticleCount
      write(outUnit, *)
    end if
        
    ! Needs to save sequence number. This is unique 
    ! across all particle groups and while reading 
    ! SRC and IC packages, sequence number should 
    ! start from the already defined value.
    this%currentSeqNumber = seqNumber

    ! TrackingOptions data
    !allocate(this%TrackingOptions) ! Moved up 
    ! Initialize defaults
    this%TrackingOptions%DebugMode = .false.
    this%TrackingOptions%BackwardTracking = .false.
    this%TrackingOptions%CreateTrackingLog = .false.
    this%TrackingOptions%StopAtWeakSinks = .false.
    this%TrackingOptions%StopAtWeakSources = .false.
    this%TrackingOptions%ExtendSteadyState = .true.
    this%TrackingOptions%SpecifyStoppingTime = .false.
    this%TrackingOptions%SpecifyStoppingZone = .false.
    this%TrackingOptions%StopTime = this%StopTime
    this%TrackingOptions%StopZone = this%StopZone
    ! Set specific option values
    if(this%TrackingDirection .eq. 2) this%TrackingOptions%BackwardTracking = .true.
    if(this%WeakSinkOption .eq. 2) this%TrackingOptions%StopAtWeakSinks = .true.
    if(this%WeakSourceOption .eq. 2) this%TrackingOptions%StopAtWeakSources = .true.
    if(this%StoppingTimeOption .ne. 2) this%TrackingOptions%ExtendSteadyState = .false.
    if(this%StoppingTimeOption .eq. 3) this%TrackingOptions%SpecifyStoppingTime = .true.
    if(this%ZoneDataOption .eq. 1) this%TrackingOptions%SpecifyStoppingZone = .true. 
    if(this%TimeseriesOutputOption .eq. 2) this%TrackingOptions%skipTimeseriesWriter = .true.

    ! Set flag to indicate whether dispersion 
    ! should be updated for different particles or not.
    ! It should be done only if SolutesOption indicates 
    ! multidispersion and simulation is RWPT. 
    if ( & 
      ( this%SolutesOption .eq. 1 ) .and. & 
      ( this%TrackingOptions%RandomWalkParticleTracking ) ) then 
      this%shouldUpdateDispersion = .true.
    end if 

    ! Set flag to indicate whether flow model is MF6
    if ( &
      ( grid%GridType .eq. 3 ) .or. &
      ( grid%GridType .eq. 4 ) ) then 
      this%isMF6 = .true.
    end if 

    ! flush
    flush(outUnit)

  end subroutine pr_ReadData


  ! Read specific GPKDE data
  subroutine pr_ReadGPKDEData( this, gpkdeFile, gpkdeUnit, outUnit )
    use UTL8MODULE,only : urword, ustop
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    ! input 
    class(ModpathSimulationDataType), target :: this
    character(len=200), intent(in)           :: gpkdeFile
    integer, intent(in)                      :: gpkdeUnit
    integer, intent(in)                      :: outUnit
    ! local
    integer :: isThisFileOpen
    integer :: icol,istart,istop,n,m
    doubleprecision    :: r
    character(len=200) :: line
    integer :: io
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW GPKDE file data'
    write(outUnit, '(1x,a)') '--------------------------'

    ! Verify if GPKDE unit is open 
    isThisFileOpen = -1
    inquire( file=gpkdeFile, number=isThisFileOpen )
    if ( isThisFileOpen .lt. 0 ) then 
      ! No gpkde 
      write(outUnit,'(A)') 'GPKDE reconstruction is disabled'
      ! flush
      flush(outUnit)
      return
    end if

    ! Yes gpkde 
    ! Requires a timeseries simulation
    if ( &
      (this%SimulationType .eq. 3) .or. (this%SimulationType .eq. 4) .or. &
      (this%SimulationType .eq. 5) .or. (this%SimulationType .eq. 6) ) then

      write(outUnit,'(A)') 'GPKDE reconstruction is enabled'
      this%TrackingOptions%GPKDEReconstruction = .true.
    
      ! Read gpkde output file
      read(gpkdeUnit, '(a)') line 
      icol = 1
      call urword(line,icol,istart,istop,0,n,r,0,0)
      this%TrackingOptions%gpkdeOutputFile = line(istart:istop)
      write(outUnit,'(A,A)') 'GPKDE output will be written to file: ', adjustl(trim(this%TrackingOptions%gpkdeOutputFile))
 
      ! Look for output column format 
      ! 0: bin ids, density data
      ! 1: bin ids, cell coordinates, density data
      ! 2: cell coordinates, density data
      this%TrackingOptions%gpkdeOutColFormat = 0 
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      if ( n.le.0 ) then
        write(outUnit, '(a)') 'GPKDE output column format default to bin ids and density data.' 
        this%TrackingOptions%gpkdeOutColFormat = 0 
      else
        select case(n)
        case(1)
          write(outUnit, '(a,I10)') 'GPKDE output column format write bin ids, cell coordinates and density data.' 
          this%TrackingOptions%gpkdeOutColFormat = n 
        case(2)
          write(outUnit, '(a,I10)') 'GPKDE output column format write bin ids, cell coordinates and density data.' 
          this%TrackingOptions%gpkdeOutColFormat = n 
        case default
          write(outUnit, '(a,I10)') 'GPKDE output column format not valid, default to bin ids and density data.' 
        end select 
      end if 

      ! Look for output file format 
      ! 0: text plain 
      ! 1: binary 
      this%TrackingOptions%gpkdeOutFileFormat = 0 
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      if ( n.le.0 ) then
        write(outUnit, '(a,I10)') 'GPKDE output file is written as text-plain.' 
        this%TrackingOptions%gpkdeOutFileFormat = 0 
      else
        select case(n)
        case(1)
          write(outUnit, '(a,I10)') 'GPKDE output file is binary' 
          this%TrackingOptions%gpkdeOutFileFormat = n
        case default
          write(outUnit, '(a,I10)') 'GPKDE output file defaults to text-plain.' 
        end select 
      end if 

      ! Read domainOrigin
      read(gpkdeUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      this%TrackingOptions%gpkdeDomainOrigin(1) = r
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      this%TrackingOptions%gpkdeDomainOrigin(2) = r
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      this%TrackingOptions%gpkdeDomainOrigin(3) = r

      ! Read domainSize
      read(gpkdeUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      this%TrackingOptions%gpkdeDomainSize(1) = r
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      this%TrackingOptions%gpkdeDomainSize(2) = r
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      this%TrackingOptions%gpkdeDomainSize(3) = r

      ! Health control
      if ( any(this%TrackingOptions%gpkdeDomainSize .lt. 0.0) ) then 
        write(outUnit,'(A)') 'One of the GPKDE domain sizes is negative. They should be positive.'
        call ustop('One of the GPKDE domain sizes is negative. They should be positive. Stop.')
      end if 
      if ( all(this%TrackingOptions%gpkdeDomainSize .eq. 0.0) ) then
        ! No gpkde 
        write(outUnit,'(A)') 'All the domain dimensions for GPKDE are zero, will disable spatial reconstruction.'
        write(outUnit,'(A)') 'GPKDE reconstruction is disabled'
        this%TrackingOptions%GPKDEReconstruction = .false.
        ! flush
        flush(outUnit)
        return
      end if 

      ! Look for grid allocation format 
      ! 0: allocate with domain grid size
      ! 1: allocate according to particle positions
      this%TrackingOptions%gpkdeGridAllocFormat = 0
      this%TrackingOptions%gpkdeAdaptGridToCoords = .false.
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      if ( n.le.0 ) then
        write(outUnit, '(a)') 'GPKDE grid allocated according to domain size.' 
      else
        select case(n)
        case(1)
          write(outUnit, '(a)') 'GPKDE grid allocated according to the particle distribution.' 
          this%TrackingOptions%gpkdeGridAllocFormat = n
          this%TrackingOptions%gpkdeAdaptGridToCoords = .true.

          ! Look for border fraction only if grid allocation format = 1
          this%TrackingOptions%gpkdeGridBorderFraction= 0.05
          call urword(line, icol, istart, istop, 3, n, r, 0, 0)
          if ( r.ne.0.0 ) then
            ! Bound border fraction to be between 0 and 1 
            if ( (r.lt.0.0).or.(r.gt.1.0) ) then 
              write(outUnit,'(a)') 'Border fraction should be between 0 and 1. Will remain as default.'
            else
              ! Ok 
              this%TrackingOptions%gpkdeGridBorderFraction = r
            end if
            write(outUnit,'(a,es18.9e3)') 'Domain border fraction set to :', this%TrackingOptions%gpkdeGridBorderFraction
          end if 
        case default
          write(outUnit, '(a)') 'GPKDE grid allocated according to domain size.' 
        end select 
      end if 

      ! Read binSize
      read(gpkdeUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      this%TrackingOptions%gpkdeBinSize(1) = r
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      this%TrackingOptions%gpkdeBinSize(2) = r
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      this%TrackingOptions%gpkdeBinSize(3) = r
     
      ! Health control
      if ( any(this%TrackingOptions%gpkdeBinSize.lt.0d0) ) then 
        write(outUnit,'(A)') 'One of the GPKDE bin sizes is negative. They should be positive.'
        call ustop('One of the GPKDE bin sizes is negative. They should be positive. Stop.')
      end if 
      if ( all(this%TrackingOptions%gpkdeBinSize.eq.0d0) ) then
        ! No gpkde 
        write(outUnit,'(A)') 'All the bin sizes are zero, will disable spatial reconstruction.'
        write(outUnit,'(A)') 'GPKDE reconstruction is disabled'
        this%TrackingOptions%GPKDEReconstruction = .false.
        ! flush
        flush(outUnit)
        return
      end if 
      ! Set binVolume, cannot be zero
      this%TrackingOptions%gpkdeBinVolume = product(&
          this%TrackingOptions%gpkdeBinSize, mask=this%TrackingOptions%gpkdeBinSize.ne.0d0)

      ! Read if sliced reconstruction
      read(gpkdeUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      select case(n)
      case(1)
        this%TrackingOptions%gpkdeSlicedReconstruction = .true.
      case default
        this%TrackingOptions%gpkdeSlicedReconstruction = .false.
      end select

      ! Read the sliced dimension
      if ( this%TrackingOptions%gpkdeSlicedReconstruction ) then
        call urword(line, icol, istart, istop, 2, n, r, 0, 0)
        if ( n.gt.0 ) then 
          select case(n)
          case(1,2,3)
            this%TrackingOptions%gpkdeSlicedDimension = n
          case default
            ! disable
            this%TrackingOptions%gpkdeSlicedDimension = 0
            this%TrackingOptions%gpkdeSlicedReconstruction = .false.
          end select
        else if (n.eq.0) then
          ! default to dimension z 
          this%TrackingOptions%gpkdeSlicedDimension = 3
        else
          ! disable
          this%TrackingOptions%gpkdeSlicedDimension = 0
          this%TrackingOptions%gpkdeSlicedReconstruction = .false.
        end if 
      end if 
      if ( this%TrackingOptions%gpkdeSlicedReconstruction ) then 
        write(outUnit,'(A,I2,A)') 'Reconstruction is configured to be sliced in dimension ', &
                this%TrackingOptions%gpkdeSlicedDimension, ' if possible.' 
      else
        write(outUnit,'(A)') 'Reconstruction is not specifically sliced.'
      end if 

      ! Read nOptimizationLoops
      read(gpkdeUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      if (n.lt.0) then
        write(outUnit,'(a)') 'Given number of optimization loops is less than 0. Defaults to 0.'
        this%TrackingOptions%gpkdeNOptLoops = 0
      else
        this%TrackingOptions%gpkdeNOptLoops = n
      end if

      ! Skip error convergence ?
      ! 0: Break if convergence criteria is met 
      ! 1: Skip and run nOptLoops optimization loops 
      this%TrackingOptions%gpkdeSkipError = .false.
      read(gpkdeUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      select case(n)
      case(0)
        this%TrackingOptions%gpkdeSkipError = .false.
        call urword(line, icol, istart, istop, 3, n, r, 0, 0)
        if ( r.lt.0.0 ) then 
          this%TrackingOptions%gpkdeRelErrorConvergence = 0.02
          write(outUnit,'(a,es18.9e3)') 'Relative error convergence cannot be negative, default to : ', &
          this%TrackingOptions%gpkdeRelErrorConvergence
        else 
          this%TrackingOptions%gpkdeRelErrorConvergence = r
          write(outUnit,'(a,es18.9e3)') 'Relative error convergence set to: ',& 
          this%TrackingOptions%gpkdeRelErrorConvergence
        end if 
      case(1)
        this%TrackingOptions%gpkdeSkipError = .true.
        write(outUnit,'(a)') 'GPKDE will run until the maximum number of optimization loops.'
      case default
        write(outUnit,'(a)') 'Invalid skip error convergence parameter, will remain as false.'
      end select

      ! Read the kernel database parameter
      ! 0: without kernel database, brute force
      ! 1: with kernel database and read parameters
      read(gpkdeUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      if (n.eq.0) then 
        this%TrackingOptions%gpkdeKernelDatabase = .false.
      else
        this%TrackingOptions%gpkdeKernelDatabase = .true.
      end if
      
      ! Read the bound kernel size format
      ! 0: bounded by domain restrictions
      ! 1: bounded by minhl and maxhl 
      ! 2: unbounded
      this%TrackingOptions%gpkdeBoundKernelSize = 0
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      select case(n)
      case(1)
        this%TrackingOptions%gpkdeBoundKernelSize = n
        write(outUnit,'(A)') 'Kernel size will be bounded by MinHL and MaxHL.'
      case(2)
        this%TrackingOptions%gpkdeBoundKernelSize = n
        write(outUnit,'(A)') 'Kernel size will be unbounded.'
      case default
        write(outUnit,'(A)') 'Kernel size will be bounded by domain restrictions.'
      end select 

      ! Read the kernels format 
      ! 0: anisotropic kernels
      ! 1: isotropic kernels 
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      select case(n)
      case(1)
        this%TrackingOptions%gpkdeIsotropicKernels = .true.
      case default 
        this%TrackingOptions%gpkdeIsotropicKernels = .false.
      end select 

      ! Interpret KDB params
      if ( this%TrackingOptions%gpkdeKernelDatabase ) then 
        write(outUnit,'(A)') 'GPKDE reconstruction with kernel database'
        ! Read kernel database params
        ! - min   h/lambda
        ! - delta h/lambda
        ! - max   h/lambda
        read(gpkdeUnit, '(a)') line
        icol = 1
        call urword(line, icol, istart, istop, 3, n, r, 0, 0)
        this%TrackingOptions%gpkdeKDBParams(1) = r
        call urword(line, icol, istart, istop, 3, n, r, 0, 0)
        this%TrackingOptions%gpkdeKDBParams(2) = r
        call urword(line, icol, istart, istop, 3, n, r, 0, 0)
        this%TrackingOptions%gpkdeKDBParams(3) = r
      else
        write(outUnit,'(A)') 'GPKDE reconstruction with brute force, no kernel database'
        this%TrackingOptions%gpkdeKernelDatabase = .false.
        if ( this%TrackingOptions%gpkdeBoundKernelSize .eq. 1 ) then 
          ! Read kernel database params
          ! - min   h/lambda
          ! - max   h/lambda
          read(gpkdeUnit, '(a)') line
          icol = 1
          call urword(line, icol, istart, istop, 3, n, r, 0, 0)
          this%TrackingOptions%gpkdeKDBParams(1) = r
          ! call urword(line, icol, istart, istop, 3, n, r, 0, 0) ! don't read deltaHL
          this%TrackingOptions%gpkdeKDBParams(2) = 0.0
          call urword(line, icol, istart, istop, 3, n, r, 0, 0)
          this%TrackingOptions%gpkdeKDBParams(3) = r
        else
          this%TrackingOptions%gpkdeKDBParams(:) = 0.0 ! all as zero
        end if 
      end if 
      if ( this%TrackingOptions%gpkdeIsotropicKernels ) then 
        write(outUnit,'(A)') 'GPKDE reconstruction with isotropic kernels' 
      end if 

      ! Read the initial smoothing format 
      ! 0: automatic selection from the global expression of Silverman (1986)
      ! 1: as a factor multiplying the bin size 
      read(gpkdeUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      this%TrackingOptions%gpkdeInitialSmoothingFormat = 0
      this%TrackingOptions%gpkdeBinSizeFactor = 5.0
      select case(n)
      case(1)
        this%TrackingOptions%gpkdeInitialSmoothingFormat =  n
        write(outUnit,'(A)') 'Initial kernel size as a factor multiplying bin size.'
        call urword(line, icol, istart, istop, 3, n, r, 0, 0)
        if ( r.le.0.0 ) then 
          write(outUnit,'(A)') 'Given initial bin size factor is less or equal to zero, will default to automatic global selection.' 
          this%TrackingOptions%gpkdeInitialSmoothingFormat = 0
        else
          this%TrackingOptions%gpkdeBinSizeFactor = r
          write(outUnit,'(A,es18.9e3)') 'Bin size factor set to : ', this%TrackingOptions%gpkdeBinSizeFactor
        end if 
      case default 
        write(outUnit,'(A)') 'Initial kernel size selected from automatic global selection of Silverman (1986).' 
      end select

      ! Read kind of reconstruction output
      ! 0: as total mass density. Smoothed phi*R*c_r
      ! 1: as resident concentration
      read(gpkdeUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      if (n.eq.0) then 
        write(outUnit,'(A)') 'GPKDE output is expressed as smoothed total mass density.'
        this%TrackingOptions%gpkdeAsConcentration = .false.
        this%TrackingOptions%gpkdeScalingFactor =&
          1.0/(this%TrackingOptions%gpkdeBinVolume)
      else
        ! If requested as resident concentration, 
        ! verifies whether porosities and retardation 
        ! are spatially uniform.
        write(outUnit,'(A)') 'GPKDE output is requested to be expressed as resident concentration.'
        if ( this%isUniformPorosity .and. this%isUniformRetardation ) then 
          write(outUnit,'(A)') 'Porosity and retardation are spatially uniform, GPKDE output is given as concentration.'
          this%TrackingOptions%gpkdeAsConcentration = .true.
          this%TrackingOptions%gpkdeScalingFactor =&
            1.0/(this%uniformPorosity*this%uniformRetardation*this%TrackingOptions%gpkdeBinVolume)
        else
          write(outUnit,'(A)') 'Porosity and retardation are NOT spatially uniform, GPKDE output is total mass density.'
          this%TrackingOptions%gpkdeAsConcentration = .false.
          this%TrackingOptions%gpkdeScalingFactor =&
            1.0/(this%TrackingOptions%gpkdeBinVolume)
        end if
      end if

      ! effectiveWeightFormat
      ! 0: compute effective number of points at domain-level (Kish 1965,1992)
      ! 1: compute average particles weight 
      ! 2: bandwidth selection based on particle positions and final mass density reconstruction
      ! 3: bandwidth selection based on local effective particles and final mass density reconstruction
      read(gpkdeUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      write(outUnit,'(a)') 'Determining the method for reconstruction of weighted particles. '
      this%TrackingOptions%gpkdeEffectiveWeightFormat = 0
      if ( n.gt.0 ) then 
        select case(n)
        case(1)
          write(outUnit,'(a)') 'Effective weight obtained as the average over particles.'
          this%TrackingOptions%gpkdeEffectiveWeightFormat = n
        case(2)
          write(outUnit,'(a)') 'Histogram calculates both counts and weights, bandwidth selected with counts.'
          this%TrackingOptions%gpkdeEffectiveWeightFormat = n
        case(3)
          write(outUnit,'(a)') 'Histogram calculates both counts and weights, bandwidth selected with cell effective counts.'
          this%TrackingOptions%gpkdeEffectiveWeightFormat = n
        case default
          write(outUnit,'(a)') 'Default to domain-level effective weight from effective number of points (Kish, 1965,1992).'
          this%TrackingOptions%gpkdeEffectiveWeightFormat = 0
        end select
      else
        write(outUnit,'(a)') 'Default to domain-level effective weight from effective number of points (Kish, 1965,1992).'
      end if

      ! Time point option 
      ! 0: follows the timeseries output stages
      ! 1: populate output timepoints based on timepointcount and a timeinterval
      ! 2: populate output timepoints by reading an array of times
      read(gpkdeUnit, '(a)', iostat=io) line
      if ( io.lt.0 ) then
        ! If no further options given, defaults
        write(outUnit,'(a)') 'GPKDE reconstruction will be performed consistently with timeseries simulation points. '
        this%TrackingOptions%gpkdeTimePointOption = 0
        this%TrackingOptions%gpkdeSkipInitialCondition = .false.
      else
        icol = 1
        call urword(line, icol, istart, istop, 2, n, r, 0, 0)
        this%TrackingOptions%gpkdeTimePointOption = n

        ! Skip initial condition
        ! 0: don't skip, as usual
        ! 1: skip initial condition
        call urword(line, icol, istart, istop, 2, n, r, 0, 0)
        select case(n)
        case(1)
          write(outUnit,'(a)') 'Reconstruction will not be performed for tracking time zero. '
          this%TrackingOptions%gpkdeSkipInitialCondition = .true.
        case default 
          this%TrackingOptions%gpkdeSkipInitialCondition = .false.
        end select

        ! Read additional data according to timepointoption
        select case(this%TrackingOptions%gpkdeTimePointOption)
        case(1)
          write(outUnit,'(a)') 'Reconstruction will be performed on specific time points.'
          write(outUnit,'(a)') 'Will interpret timepointcount and a time interval.'
          ! time interval read in the variable r
          read(gpkdeUnit, *) this%TrackingOptions%gpkdeTimePointCount, r 
          if(this%TrackingOptions%gpkdeTimePointCount .gt. 0) then
            allocate(this%TrackingOptions%gpkdeTimePoints(this%TrackingOptions%gpkdeTimePointCount))
            this%TrackingOptions%gpkdeTimePoints(1) = r 
            do m = 2, this%TrackingOptions%gpkdeTimePointCount
              this%TrackingOptions%gpkdeTimePoints(m) = this%TrackingOptions%gpkdeTimePoints(m-1) + r
            end do
          else
            write(outUnit,'(a)') 'Given timepointcount is less than zero, will default to reconstruction following timeseries.'
            this%TrackingOptions%gpkdeTimePointOption = 0
          end if
        case(2)
          write(outUnit,'(a)') 'Reconstruction will be performed on specific time points.'
          write(outUnit,'(a)') 'Will read an array of times.'
          read(gpkdeUnit, *) this%TrackingOptions%gpkdeTimePointCount
          if(this%TrackingOptions%gpkdeTimePointCount .gt. 0) then
            allocate(this%TrackingOptions%gpkdeTimePoints(this%TrackingOptions%gpkdeTimePointCount))
            read(gpkdeUnit, *) (this%TrackingOptions%gpkdeTimePoints(m), m = 1, this%TrackingOptions%gpkdeTimePointCount)
          else
            write(outUnit,'(a)') 'Given timepointcount is less than zero, will default to reconstruction following timeseries.'
            this%TrackingOptions%gpkdeTimePointOption = 0
          end if
        case default
          write(outUnit,'(a)') 'GPKDE reconstruction will be performed consistently with timeseries simulation points. '
          this%TrackingOptions%gpkdeTimePointOption = 0
        end select
      end if
    else
      ! If simulation is not timeseries
      write(outUnit,'(A)') 'GPKDE reconstruction requires a timeseries. Will remain disabled.'
    end if

    ! flush
    flush(outUnit)

    ! Close gpkde data file
    close( gpkdeUnit )

  end subroutine pr_ReadGPKDEData


  ! Read specific OBS data
  subroutine pr_ReadOBSData( this, obsFile, obsUnit, outUnit, grid )
    use UTL8MODULE,only : urword,ustop,ugetnode,u3dintmp, u3dintmpusg
    use ObservationModule, only: ObservationType
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    ! input 
    class(ModpathSimulationDataType), target :: this
    character(len=200), intent(in)           :: obsFile
    integer, intent(in)                      :: obsUnit
    integer, intent(in)                      :: outUnit
    class(ModflowRectangularGridType),intent(in) :: grid
    ! local
    integer :: nObservations  = 0
    type( ObservationType ), pointer :: obs => null()
    integer :: readStyle, nobs, no, cellNumber, layer, row, column, ocount
    integer,dimension(:),allocatable :: obsCells
    integer,dimension(:),allocatable :: cellsPerLayer
    integer :: isThisFileOpen
    integer :: icol,istart,istop,n
    integer :: ioInUnit = 0
    doubleprecision    :: r
    character(len=200) :: line
    character(len=24)  :: aname(1)
    DATA aname(1) /'                OBSCELLS'/
    integer :: baseObsUnit    = 1000
    integer :: baseObsAuxUnit = 4000 
    integer :: baseObsRecUnit = 7000 
    character(len=16) :: tempChar='tmp'
    character(len=16) :: recordChar='rec'
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW OBS file data'
    write(outUnit, '(1x,a)') '--------------------------'

    ! Verify if OBS unit is open 
    isThisFileOpen = -1
    inquire( file=obsFile, number=isThisFileOpen )
    if ( isThisFileOpen .lt. 0 ) then 
      ! No obs
      write(outUnit,'(A)') 'No observations were specified'
      ! flush
      flush(outUnit)
      return
    end if

    ! OBS unit is open

    ! Read the number of observation cells
    read(obsUnit, *, iostat=ioInUnit) line
    icol = 1
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    if ( n .le. 0 ) then 
      ! No obs
      write(outUnit,'(1X,A,I6,A)') 'Given number of observations: ', n, '. Default to no observations'
    else
      ! ok, initialize
      nObservations = n
      this%TrackingOptions%anyObservation = .true.
      write(outUnit,'(A,I5)')'Given number of observations: ', nObservations

      ! Allocate observation arrays
      call this%TrackingOptions%InitializeObservations( nObservations )

      ! Allocate id arrays in tracking options
      if(allocated(this%TrackingOptions%isObservation)) & 
          deallocate(this%TrackingOptions%isObservation)
      allocate(this%TrackingOptions%isObservation(grid%CellCount))
      if(allocated(this%TrackingOptions%idObservation)) & 
          deallocate(this%TrackingOptions%idObservation)
      allocate(this%TrackingOptions%idObservation(grid%CellCount))
      this%TrackingOptions%isObservation(:) = .false.
      this%TrackingOptions%idObservation(:) = -999

      ! Read observation cells and assign 
      ! proper variables
      do nobs = 1, nObservations

        ! A pointer
        obs => this%TrackingOptions%Observations(nobs) 

        ! Read obs name/stringid
        read(obsUnit, '(a)') line
        icol = 1
        call urword(line,icol,istart,istop,0,n,r,0,0)
        obs%stringid = line(istart:istop)

        ! Assign it increasingly, internal id
        obs%id = nobs 

        ! Report which OBS is being read 
        write(outUnit,'(A,A)') 'Reading OBS specification: ', trim(adjustl(obs%stringid))

        ! Interpret observation kind
        read(obsUnit, '(a)') line
        icol = 1
        call urword(line,icol,istart,istop,1,n,r,0,0)
        obs%stylestringid = line(istart:istop)

        ! Assing the observation style
        select case(obs%stylestringid)
        case('RES','RESIDENT')
          obs%style = 1
          ! Report
          write(outUnit,'(A)') 'Observation is of RESIDENT concentration type.'
          ! Enable flag
          this%TrackingOptions%anyResObservation = .true.
        case('SINK','FLUX')
          obs%style = 2
          ! Report 
          write(outUnit,'(A)') 'Observation is of SINK type for flux concentrations.'
          ! Enable flag
          this%TrackingOptions%anySinkObservation = .true.
        end select

        ! Read observation filename 
        read(obsUnit, '(a)') obs%outputFileName
        icol = 1
        call urword(obs%outputFileName,icol,istart,istop,0,n,r,0,0)
        obs%outputFileName = obs%outputFileName(istart:istop)

        ! Read output option
        read(obsUnit, '(a)', iostat=ioInUnit) line
        icol = 1
        call urword(line, icol, istart, istop, 2, n, r, 0, 0)
        obs%outputOption = n 

        ! Do postprocess ?
        obs%doPostprocess = .true.
        if ( obs%outputOption .eq. 0 ) obs%doPostprocess = .false.
        if ( obs%doPostprocess ) then 
          ! Read postprocess option
          call urword(line, icol, istart, istop, 2, n, r, 0, 0)
          obs%postprocessOption = n 
        end if 


        ! Continue reading if do postprocess
        if ( obs%doPostprocess ) then

          ! Histogram option for sink observations
          if ( obs%style.eq.2 ) then 
            call urword(line, icol, istart, istop, 2, n, r, 0, 0)
            select case(n)
            case (0,1)
              obs%histogramOptions = n
              write(outUnit,'(a)') 'Will interpret histogram options for the observation.'
            case default
              write(outUnit,'(a)') 'Invalid histogram options value, will remain disabled.'
              obs%histogramOptions = 0
            end select 
          end if 


          ! Look for reconstruction options
          call urword(line, icol, istart, istop, 2, n, r, 0, 0)
          select case(n)
          case (0,1)
            obs%reconstructionOptions = n
            write(outUnit,'(a)') 'Will interpret reconstruction options for the observation.'
          case default
            write(outUnit,'(a)') 'Invalid reconstruction options value, will remain disabled.'
            obs%reconstructionOptions = 0
          end select


          ! Process histogram options
          select case(obs%histogramOptions)
          case(1)
            write(outUnit,'(a)') 'Interpreting histogram options.'

            ! Histogram option for sink observations
            ! 0: Scott's rule for bin size
            ! 1: Freedman-Diaconis rule for bin size
            ! 2: Given bin size
            ! 3: Follow time step from timeseries run
            read(obsUnit, '(a)') line
            icol = 1
            call urword(line, icol, istart, istop, 2, n, r, 0, 0)
            select case(n)
            case(0)
              write(outUnit,'(a)') "Will estimate the histogram bin with Scott's rule."
              obs%histogramBinFormat = n
            case(1)
              write(outUnit,'(a)') "Will estimate the histogram bin with Freedman-Diaconis rule."
              obs%histogramBinFormat = n
            case(2)
              write(outUnit,'(a)') "Will interpret bin size for histogram bin."
              obs%histogramBinFormat = n
            case(3)
              write(outUnit,'(a)') "Will employ the timeseries time step for histogram bin."
              obs%histogramBinFormat = n
            case default
              write(outUnit,'(a)') 'Invalid bin option, will remain as default.'
              obs%histogramBinFormat = 1 ! Freedman diaconis rule for bin size selection
            end select

            select case(obs%histogramBinFormat)
            case(0,1)
             ! If format is 0,1 loof for bin size fraction
             call urword(line, icol, istart, istop, 3, n, r, 0, 0)
             if ( (r.gt.0d0).and.(r.le.1d0) ) then
              obs%binSizeFraction = r
              write(outUnit,'(a,es18.9e3)') 'Bin size fraction for automatic bin selection ', & 
                      obs%binSizeFraction
             else
              write(outUnit,'(a,es18.9e3)') 'Bin size fraction for automatic bin selection will remain as default ', & 
                      obs%binSizeFraction
             end if
            case(2)
             ! If format is 2, look for bin size for histogram
             call urword(line, icol, istart, istop, 3, n, r, 0, 0)
             if ( r.gt.0d0 ) then
               obs%histogramBin = r
             else
              write(outUnit,'(a)') 'Invalid value for histogram bin size, will fallback to default bin option.'
              obs%histogramBinFormat = 1
             end if
            case(3)
              ! do nothing
            end select


            ! Look for time step out, in all cases of bin option
            call urword(line, icol, istart, istop, 3, n, r, 0, 0)
            if ( r.gt.0d0 ) then
              write(outUnit,'(a,es18.9e3)') 'Will try to interpolate the observation series to the time step ', r 
              obs%timeStepOut = r 
            else
              write(outUnit,'(a)') 'Time step out is zero, will not perform interpolation for observation.' 
              obs%timeStepOut = 0d0
            end if

          case default
            write(outUnit,'(a)') 'Will not interpret histogram options for this observation.'
          end select

          
          ! interpret reconstruction options 
          !call urword(line, icol, istart, istop, 2, n, r, 0, 0)
          select case(obs%reconstructionOptions)
          case(1)
            write(outUnit,'(a)') 'Interpret reconstruction parameters for the observation.'
            ! Read opt loops
            read(obsUnit, '(a)', iostat=ioInUnit) line
            icol = 1
            call urword(line, icol, istart, istop, 2, n, r, 0, 0)
            if (n.ge.0) then
              obs%nOptLoops = n  
            else
              obs%nOptLoops = 10 
              write(outUnit,'(a,I3)') 'Invalid number of optimization loops, will default to ', obs%nOptLoops
            end if 

            ! Error convergence
            call urword(line, icol, istart, istop, 3, n, r, 0, 0)
            if (r.gt.0d0) then
              obs%errorConvergence = r
            else 
              obs%errorConvergence = 0.01d0
              write(outUnit,'(a,es18.9e3)') 'Invalid error convergence, will default to ', obs%errorConvergence
            end if
 
            ! initialsmoothingformat
            read(obsUnit, '(a)', iostat=ioInUnit) line
            icol = 1
            call urword(line, icol, istart, istop, 2, n, r, 0, 0)
            select case(n)
            case(0)
              obs%initialSmoothingFormat = n
            case(1)
              obs%initialSmoothingFormat = n
            case default
              obs%initialSmoothingFormat = 0
              write(outUnit,'(a,I3)') 'Invalid smoothing format will default to ', obs%initialSmoothingFormat
            end select

            ! Read bin size factor
            if ( obs%initialSmoothingFormat.eq.1) then 
              call urword(line, icol, istart, istop, 3, n, r, 0, 0)
              if ( r.gt.0d0 ) then 
                obs%binSizeFactor = r
              else
                obs%binSizeFactor = 3d0
                write(outUnit,'(a,es18.9e3)') 'Bin size factor will default to ', obs%binSizeFactor
              end if 
            end if

            ! effectiveweightformat
            read(obsUnit, '(a)', iostat=ioInUnit) line
            icol = 1
            call urword(line, icol, istart, istop, 2, n, r, 0, 0)
            select case(n)
            case(0,1,2,3)
              obs%effectiveWeightFormat = n
            case default
              write(outUnit,'(a,I3)') 'Invalid effective weight format will default to ', obs%effectiveWeightFormat
            end select

          case default
            write(outUnit,'(A)') 'Will not interpret reconstruction parameters for observation.' 
          end select

        end if 


        ! Assign output units, according to 
        ! obs output option
        select case ( obs%outputOption )
        case(0)
          ! Only records, no postprocess
          ! Copy outputFileName to recOutputFileName
          obs%recOutputUnit = baseObsRecUnit + nobs
          obs%recOutputFileName = obs%outputFileName
        case(1,2)
          ! Records and postprocess, also initialize aux params
          ! cases 1 and 2, only aux unit is initialized as scratch
          ! case 2 both records and aux units are initailized as scratch
          obs%outputUnit = baseObsUnit + nobs
          obs%recOutputUnit = baseObsRecUnit + nobs
          write( unit=obs%recOutputFileName, fmt='(a)')&
              trim(adjustl(recordChar))//'_'//trim(adjustl(obs%outputFileName))
          obs%auxOutputUnit = baseObsAuxUnit + nobs
          write( unit=obs%auxOutputFileName, fmt='(a)')&
              trim(adjustl(tempChar))//'_'//trim(adjustl(obs%outputFileName))
        case default
          call ustop('Invalid output option for observation. Stop.')
        end select


        ! Read observation cell option
        ! Determine how to read cells
        read(obsUnit, '(a)', iostat=ioInUnit) line
        icol = 1
        call urword(line, icol, istart, istop, 2, n, r, 0, 0)
        obs%cellOption = n 


        ! Load observation cells
        select case( obs%cellOption )
          ! In case 0, a list of cell ids is specified, that 
          ! compose the observation.  
          case (0)
            ! Read number of observation cells 
            call urword(line, icol, istart, istop, 2, n, r, 0, 0)
            if ( n.lt.1 ) then 
              call ustop('Given number of cells for observation is .lt. 1. Should be at least 1. Stop.')
            end if 
            obs%nCells = n 
            
            ! Are these ids as (lay,row,col) or (cellid) ?
            call urword(line, icol, istart, istop, 2, n, r, 0, 0)
            readStyle = n

            ! Depending on the number of cells 
            ! allocate array for cell ids
            if ( allocated( obs%cells ) ) deallocate( obs%cells )
            allocate( obs%cells(obs%nCells) )

            ! Observation of RESIDENT concentration cannot directly
            ! count inside particles loop. Needs some reduction 
            ! spec for consistent counting do to OpenMP
            !if ( allocated( obs%nRecordsCell ) ) deallocate( obs%nRecordsCell )
            !allocate( obs%nRecordsCell(obs%nCells) )
            !obs%nRecordsCell(:) = 0

            ! Load the observation cells

            if( readStyle .eq. 0) then
              ! Validate grid type for reading as layer, row, column
              ! GridType .eq. 1 = RectangularGridDis       (MF2005)
              ! GridType .eq. 2 = RectangularGridDisuMfusg (MFUSG)
              ! GridType .eq. 3 = RectangularGridDisMf6    (MF6DIS)
              ! GridType .eq. 4 = RectangularGridDisvMf6   (MF6DISV)
              ! GridType .eq. 5 = RectangularGridDisuMf6   (MF6DISU)
              select case(grid%GridType)
              case(0,2,4,5) ! mfusg, disv, disu
              write(outUnit,'(A)') 'Error: reading cells as (lay,row,col) only allowed for structured grids (DIS).'
              call ustop('Error: reading cells as (lay,row,col) only allowed for structured grids (DIS).')
              end select
              ! Read as layer, row, column
              do no = 1, obs%nCells
                read(obsUnit, *) layer, row, column
                call ugetnode(& 
                  grid%LayerCount, & 
                  grid%RowCount,   & 
                  grid%ColumnCount,&
                  layer, row, column,cellNumber)
                obs%cells(no) = cellNumber
              end do 
            else if ( readStyle .eq. 1 ) then 
              ! Read as cell number
              do no = 1, obs%nCells
                read(obsUnit,*)  cellNumber
                obs%cells(no) = cellNumber
              end do 
            else
              call ustop('Invalid cell reading style for this observation. Stop.')
            end if

          case (1)
            ! In case 1, observation cells are given by specifying a 3D array
            ! with 0 (not observation) and 1 (observation) 

            ! Required for u3d
            if(allocated(obsCells)) deallocate(obsCells)
            allocate(obsCells(grid%CellCount))
            obsCells(:) = 0
         
            ! Read cells
            if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
              call u3dintmp(obsUnit, outUnit, grid%LayerCount, grid%RowCount,      &
                grid%ColumnCount, grid%CellCount, obsCells, aname(1)) 
            else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
              if ( .not. allocated( cellsPerLayer ) ) then
                allocate(cellsPerLayer(grid%LayerCount))
                do n = 1, grid%LayerCount
                  cellsPerLayer(n) = grid%GetLayerCellCount(n)
                end do
              end if 
              call u3dintmpusg(obsUnit, outUnit, grid%CellCount, grid%LayerCount,  &
                obsCells, aname(1), cellsPerLayer)
            else
              write(outUnit,*) 'Invalid grid type specified when reading OBSCELLS array data.'
              write(outUnit,*) 'Stopping.'
              call ustop(' ')          
            end if

            ! Count how many obs cells specified 
            obs%nCells = count(obsCells/=0)

            if ( obs%nCells .eq. 0 ) then 
              write(outUnit,*) 'No observation cells in the array of cells for observation:', obs%id
              write(outUnit,*) 'Stopping.'
              call ustop('No observation cells in the array of cells. Stop.')
            end if

            ! Depending on the number of cells 
            ! allocate array for cell ids
            if ( allocated( obs%cells ) ) deallocate( obs%cells )
            allocate( obs%cells(obs%nCells) )

            !if ( allocated( obs%nRecordsCell ) ) deallocate( obs%nRecordsCell )
            !allocate( obs%nRecordsCell(obs%nCells) )
            !obs%nRecordsCell(:) = 0

            ! Fill obs%cells with the corresponding cell numbers
            ocount = 0
            do n =1,grid%CellCount
              if(obsCells(n).eq.0) cycle
              ocount = ocount + 1
              obs%cells(ocount) = n 
            end do

          case default
            ! Invalid option
            call ustop('Invalid observation cells reading option. Stop.')
        end select

        ! Assign into id arrays
        do no =1, obs%nCells
          this%TrackingOptions%isObservation(obs%cells(no)) = .true.
          ! The id on the list of cells !
          this%TrackingOptions%idObservation(obs%cells(no)) = nobs
        end do

      end do 

      ! Close the OBS unit
      close( obsUnit ) 

    end if  

    ! flush
    flush(outUnit)


  end subroutine pr_ReadOBSData


  ! Read specific RWOPTS data
  subroutine pr_ReadRWOPTSData( this, rwoptsFile, rwoptsUnit, outUnit )
    use UTL8MODULE,only : urword,ustop,u3dintmpusg, u3dintmp
    use CompilerVersion,only : get_compiler
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    ! input 
    class(ModpathSimulationDataType), target     :: this
    character(len=200), intent(in)               :: rwoptsFile
    integer, intent(in)                          :: rwoptsUnit
    integer, intent(in)                          :: outUnit
    ! local
    type(ParticleTrackingOptionsType), pointer :: trackingOptions
    integer :: isThisFileOpen
    integer :: icol,istart,istop,n,nd,currentDim,dcount
    doubleprecision    :: r
    character(len=200) :: line
    character(len=24),dimension(1) :: aname
    data aname(1) /'       ICBOUND'/
    character(len=10) :: compiler
    integer :: io
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW RWOPTS file data'
    write(outUnit, '(1x,a)') '---------------------------'

    ! Verify if unit is open 
    isThisFileOpen = -1
    inquire( file=rwoptsFile, number=isThisFileOpen )
    if ( isThisFileOpen .lt. 0 ) then 
      ! No rwopts file
      write(outUnit,'(A)') 'RWOPTS were not specified in name file and are required for RW simulation.'
      call ustop('RWOPTS were not specified in name file and are required for RW simulation. Stop.')
    end if

    ! Pointer to this%TrackingOptions
    trackingOptions => this%TrackingOptions

    ! Time Step kind 
    read(rwoptsUnit, '(a)') line
    icol = 1
    call urword(line,icol,istart,istop,1,n,r,0,0)
    line = line(istart:istop)
    ! Advection 
    if ( line .eq. 'ADV' ) then
      trackingOptions%timeStepKind = 1
      read( rwoptsUnit, * ) line
      icol = 1
      call urword(line,icol,istart,istop,3,n,r,0,0)
      trackingOptions%timeStepParameters(1) = r
      write(outUnit,'(A)') 'RW time step will be selected with the ADV criteria.'
    ! Dispersion 
    else if ( line .eq. 'DISP' ) then
      trackingOptions%timeStepKind = 2
      read( rwoptsUnit, * ) line
      icol = 1
      call urword(line,icol,istart,istop,3,n,r,0,0)
      trackingOptions%timeStepParameters(2) = r
      write(outUnit,'(A)') 'RW time step will be selected with the DISP criteria.'
    ! Minimum between advection and dispersion 
    else if ( line .eq. 'MIN_ADV_DISP' ) then
      trackingOptions%timeStepKind = 3
      read( rwoptsUnit, * ) line
      icol = 1
      call urword(line,icol,istart,istop,3,n,r,0,0)
      trackingOptions%timeStepParameters(1) = r
      read( rwoptsUnit, * ) line
      icol = 1
      call urword(line,icol,istart,istop,3,n,r,0,0)
      trackingOptions%timeStepParameters(2) = r
      write(outUnit,'(A)') 'RW time step will be selected with the MIN_ADV_DISP criteria.'
    ! Fixed 
    else if ( line .eq. 'FIXED' ) then
      trackingOptions%timeStepKind = 4
      read( rwoptsUnit, * ) line
      icol = 1
      call urword(line,icol,istart,istop,3,n,r,0,0)
      trackingOptions%timeStepParameters(1) = r
      write(outUnit,'(A)') 'RW time step is given with FIXED criteria.'
    else
      call ustop('Invalid option for time step selection. Stop.')
    end if


    ! Advection Integration Kind
    read(rwoptsUnit, '(a)') line
    icol = 1
    call urword(line,icol,istart,istop,1,n,r,0,0)
    line = line(istart:istop)
    select case(line)
      case('EXPONENTIAL')
        trackingOptions%advectionKind = 1
        write(outUnit,'(A)') 'RW advection integration is EXPONENTIAL.'
      case('EULERIAN')
        trackingOptions%advectionKind = 2
        write(outUnit,'(A)') 'RW advection integration is EULERIAN.'
      case default
        trackingOptions%advectionKind = 2
        write(outUnit,'(A)') 'Given RW advection integration is not valid. Defaults to EULERIAN.'
    end select

    ! Read RW dimensionsmask. Determines to which dimensions apply RW displacements
    read(rwoptsUnit, '(a)') line
    
    ! X
    icol = 1
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    trackingOptions%dimensionMask(1) = n

    ! Y
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    trackingOptions%dimensionMask(2) = n

    ! Z
    call urword(line, icol, istart, istop, 2, n, r, 0, 0)
    trackingOptions%dimensionMask(3) = n

    ! Health check
    if ( any( trackingOptions%dimensionMask.gt.1 ) ) then 
      ! Invalid dimensions
      write(outUnit,'(A)') 'Invalid value for dimensions mask. Should 0 or 1.'
      call ustop('Invalid value for dimensions mask. Should 0 or 1. Stop.')
    end if 
    if ( any( trackingOptions%dimensionMask.lt.0 ) ) then 
      ! Invalid dimensions
      write(outUnit,'(A)') 'Invalid value for dimensions mask. Should 0 or 1.'
      call ustop('Invalid value for dimensions mask. Should 0 or 1. Stop.')
    end if 


    ! Set nDim
    trackingOptions%nDim = sum(trackingOptions%dimensionMask)
    if ( trackingOptions%nDim .le. 0 ) then
      ! No dimensions
      write(outUnit,'(A)') 'No dimensions were given for RW displacements at RWOPTS, nDim .eq. 0.'
      call ustop('No dimensions were given for RW displacements at RWOPTS, nDim .eq. 0. Stop.')
    end if 


    ! Save dim mask into dimensions 
    if ( allocated( trackingOptions%dimensions ) ) deallocate( trackingOptions%dimensions ) 
    allocate( trackingOptions%dimensions( trackingOptions%nDim  ) )
    dcount= 0
    do nd = 1, 3
      if ( trackingOptions%dimensionMask(nd) .eq. 0 ) cycle
      dcount = dcount + 1 
      trackingOptions%dimensions(dcount) = nd
    end do 


    ! Detect idDim and report dimensions
    ! where displacements will be applied
    select case(trackingOptions%nDim)
      ! 1D
      case(1)
        write(outUnit,'(A)') 'Random displacements for 1 dimension.'
        ! Relate x,y,z dimensions to 1 dimensions
        do nd = 1,3
          if ( trackingOptions%dimensionMask( nd ) .eq. 0 ) cycle
          select case(nd) 
            case (1)
              trackingOptions%idDim1 = nd
              write(outUnit,'(A)') 'Random displacements for X dimension.'
            case (2)
              trackingOptions%idDim1 = nd
              write(outUnit,'(A)') 'Random displacements for Y dimension.'
            case (3)
              trackingOptions%idDim1 = nd
              write(outUnit,'(A)') 'Random displacements for Z dimension.'
          end select   
          ! Use the first found
          exit
        end do
      ! 2D
      case(2)
        write(outUnit,'(A)') 'Random displacements for 2 dimensions.'
        ! Relate x,y,z dimensions to 1,2 dimensions
        do nd = 1,3
          if ( trackingOptions%dimensionMask( nd ) .eq. 0 ) cycle
          currentDim = sum( trackingOptions%dimensionMask(1:nd) )
          select case(nd) 
            case (1)
              trackingOptions%idDim1 = n
              write(outUnit,'(A)') 'Random displacements for X dimension.'
            case (2)
              if ( currentDim .eq. 1 ) then 
                trackingOptions%idDim1 = nd
              else if ( currentDim .eq. 2 ) then
                trackingOptions%idDim2 = nd
              end if
              write(outUnit,'(A)') 'Random displacements for Y dimension.'
            case (3)
              trackingOptions%idDim2 = nd
              write(outUnit,'(A)') 'Random displacements for Z dimension.'
          end select   
        end do
      ! 3D
      case(3)
        write(outUnit,'(A)') 'Random displacements for 3 dimensions.'
    end select


    ! Read the random number generator function
    read(rwoptsUnit, '(a)', iostat=io) line
    if ( io.lt.0 ) then 
      ! assume not given, choose by compiler
      call get_compiler(compiler) 
      select case( compiler ) 
      case('GFORTRAN')
        trackingOptions%randomGenFunction = 1 
      case('IFORT')
        trackingOptions%randomGenFunction = 2 
      case default
        trackingOptions%randomGenFunction = 2 
      end select
    else
      icol = 1
      call urword(line,icol,istart,istop,2,n,r,0,0)
      select case(n)
      case(0)
        ! assume not given, choose by compiler
        call get_compiler(compiler) 
        select case( compiler ) 
        case('GFORTRAN')
          trackingOptions%randomGenFunction = 1 
        case('IFORT')
          trackingOptions%randomGenFunction = 2 
        case default
          trackingOptions%randomGenFunction = 2 
        end select
      case(1)
        ! random_number
        trackingOptions%randomGenFunction = 1 
      case(2)
        ! rgn_par_zig
        trackingOptions%randomGenFunction = 2
      case default
        write(outUnit,'(A)') 'Selected random generator function not implemented. Stop.'
        call ustop('Selected random generator function not implemented. Stop.')
      end select 
    end if 

    ! flush
    flush(outUnit)

    ! Close rwopts data file
    close( rwoptsUnit )


  end subroutine pr_ReadRWOPTSData


  ! Read specific IC data
  subroutine pr_ReadICData( this, icFile, icUnit, outUnit, grid, porosity )
    use UTL8MODULE,only : urword,ustop,u3ddblmpusg, u3ddblmp
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    ! input 
    class(ModpathSimulationDataType), target     :: this
    character(len=200), intent(in)               :: icFile
    integer, intent(in)                          :: icUnit
    integer, intent(in)                          :: outUnit
    class(ModflowRectangularGridType),intent(in) :: grid
    doubleprecision, dimension(:), intent(in)    :: porosity
    ! local
    integer :: isThisFileOpen
    integer :: icol,istart,istop,n
    integer :: nc, nic, nd, m
    doubleprecision    :: r
    character(len=200) :: line
    integer :: nInitialConditions, nValidInitialConditions 
    integer :: initialConditionFormat
    integer :: particleArrayFormat
    integer :: newParticleGroupCount
    integer :: pgCount
    integer :: soluteId
    integer, dimension(:), pointer :: dimensionMask
    integer, pointer               :: nDim
    type(ParticleGroupType),dimension(:),allocatable :: particleGroups
    type(ParticleGroupType),dimension(:),allocatable :: newParticleGroups
    doubleprecision, dimension(:), allocatable :: densityDistribution
    integer, dimension(:), allocatable :: cellsPerLayer
    doubleprecision :: initialReleaseTime
    doubleprecision :: particleMass, effParticleMass
    doubleprecision :: cellTotalMass, totalAccumulatedMass
    doubleprecision :: totalAccumulatedParticlesMass
    doubleprecision :: cellDissolvedMass, totalDissolvedMass
    doubleprecision :: cellVolume,sX,sY,sZ,nPX,nPY,nPZ
    doubleprecision :: nParticlesCell
    doubleprecision :: avgParticleMass, varParticleMass
    doubleprecision, parameter :: nParticlesCellMin = 0.5d0 
    integer :: totalParticleCount, seqNumber, idmax, particleCount
    integer :: iNPX,iNPY,iNPZ,NPCELL
    integer :: validCellCounter, cellCounter, cellNumber
    integer, allocatable, dimension(:,:) :: subDivisions
    integer, allocatable, dimension(:)   :: validCellNumbers
    doubleprecision, allocatable, dimension(:) :: activeCellTotalMass
    logical :: saveValidCellNumbers = .false.
    character(len=24),dimension(1) :: aname
    data aname(1) /'            IC'/
    ! Interface
    procedure(MassParticlesArray), pointer :: CreateMassParticlesArray => null()
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW IC file data'
    write(outUnit, '(1x,a)') '-----------------------'

    ! Verify if unit is open
    isThisFileOpen = -1
    inquire( file=icFile, number=isThisFileOpen )
    if ( isThisFileOpen .lt. 0 ) then 
      ! No ic 
      write(outUnit,'(A)') 'IC package was not specified in name file.'
      return
    end if

    ! Preparations for interpreting IC's

    ! RW dimensionality vars
    dimensionMask => this%TrackingOptions%dimensionMask
    nDim => this%TrackingOptions%nDim

    ! Process IC's
    read(icUnit, *) nInitialConditions
    write(outUnit,'(A,I5)') 'Given number of initial conditions = ', nInitialConditions
    nValidInitialConditions = 0
    particleCount = 0 

    if(nInitialConditions .le. 0) then
      ! No ic 
      write(outUnit,'(A)') 'Number of given initial conditions is .le. 0. Leaving the function.'
      return
    end if

    ! CellsPerLayer, required for u3d reader
    allocate(cellsPerLayer(grid%LayerCount))
    do n = 1, grid%LayerCount
      cellsPerLayer(n) = grid%GetLayerCellCount(n)
    end do

    ! Carrier for candidate particle groups 
    allocate(particleGroups(nInitialConditions))

    ! Loop over initial conditions
    do nic = 1, nInitialConditions
      
      ! Report which IC will be processed
      write(outUnit,'(A,I5)') 'Processing initial condition: ', nic

      ! Increase pgroup counter
      particleGroups(nic)%Group = this%ParticleGroupCount + nic

      ! Set release time for initial condition.
      ! It is an initial condition, then
      ! assumes release at the beginning
      initialReleaseTime = 0d0
      call particleGroups(nic)%SetReleaseOption1(initialReleaseTime)

      ! Read id 
      read(icUnit, '(a)') particleGroups(nic)%Name

      ! Initial condition format
      ! All read a concentration distribution with u3d
      ! 0: consistency of mass at global level, uniform mass particles
      ! 1: consistency of mass at a cell level, mass particles with per-cell mass
      read(icUnit, '(a)') line
      icol = 1
      call urword(line,icol,istart,istop,2,n,r,0,0)
      select case(n)
        case(0,1)
          initialConditionFormat = n
        case default
          write(outUnit,*) 'Invalid initial condition format. Stop.'
          call ustop('Invalid initial condition format. Stop.')
      end select
      
      ! Read the particles location format
      ! 0: structured, equispaced
      ! 1: quasi-random
      ! 2: random
      !read(icUnit, '(a)') line
      !icol = 1
      call urword(line,icol,istart,istop,2,n,r,0,0)
      select case(n)
        ! uniform/structured
        case(0)
          CreateMassParticlesArray => CreateMassParticlesAsInternalArray
          particleArrayFormat = n
        ! quasi random
        case(1)
          CreateMassParticlesArray => CreateMassParticlesAsQuasiRandomInternalArray
          particleArrayFormat = n
        ! random
        case(2)
          CreateMassParticlesArray => CreateMassParticlesAsRandomInternalArray
          particleArrayFormat = n
        case default
          write(outUnit,'(A)') 'Given format for creating mass particles internal array is not valid. Stop.' 
          call ustop('Given format for creating mass particles internal array is not valid. Stop.' )
      end select


      ! Initialize sequence number to the current value
      seqNumber = this%currentSeqNumber
    

      ! And process initial condition 
      select case ( initialConditionFormat )
      ! Read initial condition as resident concentration (ML^-3)
      case (0) 
        ! Given a value for the mass of particles, use flowModelData to compute cellvolume
        ! and a shape factor from which the number of particles per cell is estimated.
        ! Establishes mass consistency at the domain level and all particles remain
        ! with the same mass.

        if(allocated(densityDistribution)) deallocate(densityDistribution)
        allocate(densityDistribution(grid%CellCount))

        ! Read particles mass
        read(icUnit, '(a)') line
        icol = 1
        call urword(line,icol,istart,istop,3,n,r,0,0)
        if ( r.lt.0d0 ) then 
          write(outUnit,'(A)') 'Given particle mass is negative. Should be positive. Stop.'
          call ustop('Given particle mass is negative. Should be positive. Stop.')
        end if
        particleMass = r 

        if ( ( this%ParticlesMassOption .eq. 2 ) ) then 
          ! Read solute id
          read(icUnit, '(a)') line
          icol = 1
          call urword(line,icol,istart,istop,2,n,r,0,0)
          if ( n.lt.1 ) then 
            write(outUnit,'(A)')'Given SpeciesID is less than 1. Minimum is 1. Stop.'
            call ustop('Given SpeciesID is less than 1. Minimum is 1. Stop.')
          end if
          soluteId = n
        end if

        ! Read concentration
        if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
          call u3ddblmp(icUnit, outUnit, grid%LayerCount, grid%RowCount,      &
            grid%ColumnCount, grid%CellCount, densityDistribution, aname(1))
        else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
          call u3ddblmpusg(icUnit, outUnit, grid%CellCount, grid%LayerCount,  &
            densityDistribution, aname(1), cellsPerLayer)
        else
          write(outUnit,*) 'Invalid grid type specified when reading IC array ', & 
              particleGroups(nic)%Name, ' name. Stop.'
          call ustop('Invalid grid type specified when reading IC array') 
        end if

        ! Validity of initial condition
        ! Used to allocate subdivisions 
        validCellCounter = count( densityDistribution /= 0d0 )
        if( validCellCounter .eq.0 ) then
          write(outUnit,'(a,a,a)') 'Warning: initial condition ',&
            trim(adjustl(particleGroups(nic)%Name)),' has a distribution only with zeros'
          write(outUnit,'(a)') 'It will not create a particle group. Continue to the next.'
          ! Process the next one
          cycle
        end if 

        ! Allocate subdivisions
        if (allocated(subDivisions)) deallocate(subDivisions)
        allocate( subDivisions(validCellCounter,3) )
        subDivisions(:,:) = 0

        ! If validCellCounter is not equal to the 
        ! size of the initial concentration array, 
        ! and somehow much smaller,  then it could 
        ! be of benefit to save the valid cell ids  
        ! and deallocate densityDistribution 
        saveValidCellNumbers = .false. 
        if ( validCellCounter .lt. 0.5*grid%CellCount ) then
          if ( allocated( validCellNumbers ) ) deallocate( validCellNumbers ) 
          allocate( validCellNumbers(validCellCounter) ) 
          saveValidCellNumbers = .true. 
        end if 

        ! Loop over densityDistribution and compute:
        !  - totalDissolvedMass
        !  - totalAccumulatedMass: considers retardation factor
        !  - totalParticleCount
        totalDissolvedMass   = 0d0
        totalAccumulatedMass = 0d0
        totalParticleCount   = 0
        cellCounter          = 0
        select case(particleArrayFormat)
        ! for equispaced and quasi-random
        case(0,1)
          do nc = 1, grid%CellCount
            ! If no concentration, next
            if ( densityDistribution(nc) .eq. 0d0 ) cycle

            ! Increase cellCounter
            cellCounter = cellCounter + 1

            ! Compute cell volume
            cellVolume = 1d0
            do nd=1,3
              if( dimensionMask(nd).eq.0 ) cycle
              select case(nd)
              case(1)
                cellVolume = cellVolume*grid%DelX(nc)
              case(2)
                cellVolume = cellVolume*grid%DelY(nc)
              case(3)
                ! simple dZ
                cellVolume = cellVolume*(grid%Top(nc)-grid%Bottom(nc))
              end select
            end do
            cellDissolvedMass = 0d0
            cellTotalMass = 0d0
            ! Absolute value is required for the weird case that 
            ! densityDistribution contains negative values
            cellDissolvedMass = abs(densityDistribution(nc))*porosity(nc)*cellVolume
            totalDissolvedMass = totalDissolvedMass + cellDissolvedMass 
            cellTotalMass = cellDissolvedMass*this%Retardation(nc)
            totalAccumulatedMass = totalAccumulatedMass + cellTotalMass
            
            ! nParticlesCell: estimate the number of particles 
            ! for the cell using the specified mass 
            nParticlesCell = cellTotalMass/particleMass

            ! If less than a threshold, the cell has no particles
            if ( nParticlesCell .lt. nParticlesCellMin ) cycle

            ! Compute shapeFactors only if dimension is active
            ! If not, will remain as zero
            sX = 0
            sY = 0
            sZ = 0
            do nd=1,3
              if( dimensionMask(nd).eq.0 ) cycle
              select case(nd)
              case(1)
                sX = grid%DelX(nc)/(cellVolume**(1d0/dble(nDim)))
              case(2)
                sY = grid%DelY(nc)/(cellVolume**(1d0/dble(nDim)))
              case(3)
                ! simple dZ
                sZ = (grid%Top(nc)-grid%Bottom(nc))/(cellVolume**(1d0/dble(nDim)))
              end select
            end do

            ! Estimate subdivisions
            nPX    = sX*( (nParticlesCell)**(1d0/dble(nDim)) ) 
            nPY    = sY*( (nParticlesCell)**(1d0/dble(nDim)) )
            nPZ    = sZ*( (nParticlesCell)**(1d0/dble(nDim)) )
            iNPX   = max(int( nPX + 0.5d0 ),1) 
            iNPY   = max(int( nPY + 0.5d0 ),1) 
            iNPZ   = max(int( nPZ + 0.5d0 ),1) 
            NPCELL = iNPX*iNPY*iNPZ
            totalParticleCount = totalParticleCount + NPCELL

            ! Save in subdivisions
            subDivisions(cellCounter,1) = iNPX
            subDivisions(cellCounter,2) = iNPY
            subDivisions(cellCounter,3) = iNPZ

            ! Save valid cell numbers
            if ( saveValidCellNumbers ) then
              validCellNumbers(cellCounter) = nc
            end if 
          end do ! end loop over densityDistribution

        ! for random
        case(2)
          do nc = 1, grid%CellCount
            ! If no concentration, next
            if ( densityDistribution(nc) .eq. 0d0 ) cycle

            ! Increase cellCounter
            cellCounter = cellCounter + 1

            ! Compute cell volume
            cellVolume = 1d0
            do nd=1,3
              if( dimensionMask(nd).eq.0 ) cycle
              select case(nd)
              case(1)
                cellVolume = cellVolume*grid%DelX(nc)
              case(2)
                cellVolume = cellVolume*grid%DelY(nc)
              case(3)
                ! simple dZ
                cellVolume = cellVolume*(grid%Top(nc)-grid%Bottom(nc))
              end select
            end do
            cellDissolvedMass = 0d0
            cellTotalMass = 0d0
            ! Absolute value is required for the weird case that 
            ! densityDistribution contains negative values
            cellDissolvedMass = abs(densityDistribution(nc))*porosity(nc)*cellVolume
            totalDissolvedMass = totalDissolvedMass + cellDissolvedMass 
            cellTotalMass = cellDissolvedMass*this%Retardation(nc)
            totalAccumulatedMass = totalAccumulatedMass + cellTotalMass
            
            ! nParticlesCell: estimate the number of particles 
            ! for the cell using the specified mass 
            nParticlesCell = cellTotalMass/particleMass

            ! If less than a threshold, the cell has no particles
            if ( nParticlesCell .lt. nParticlesCellMin ) cycle

            ! For random distribution, does not need to 
            ! compute subdivisions, can simply take the integer
            NPCELL = int(nParticlesCell) + 1
            totalParticleCount = totalParticleCount + NPCELL

            ! Save in subdivisions
            ! NOTE: the first position contains
            ! all the particles
            subDivisions(cellCounter,1) = NPCELL
            !subDivisions(cellCounter,2) = 0 
            !subDivisions(cellCounter,3) = 0

            ! Save valid cell numbers
            if ( saveValidCellNumbers ) then
              validCellNumbers(cellCounter) = nc
            end if 
          end do ! end loop over densityDistribution

        end select

        ! Assign totalParticleCount and continue to the next IC if no particles
        particleGroups(nic)%TotalParticleCount = totalParticleCount
        if ( totalParticleCount .eq. 0 ) then 
          write(outUnit,*) ' Warning: initial condition ',&
              particleGroups(nic)%Name,' has zero particles, it will skip this group.'
          ! Process the next one
          cycle
        end if 

        ! effective particles mass
        effParticleMass = totalAccumulatedMass/dble(totalParticleCount)
        write(outUnit,'(A,es18.9e3)') 'Original particle mass for initial condition = ', particleMass
        write(outUnit,'(A,es18.9e3)') 'Effective particle mass for initial condition = ', effParticleMass

        ! Both the total number of particles and the effective mass are known, 
        ! The mass is uniform for all particles
        totalAccumulatedParticlesMass = effParticleMass*dble(totalParticleCount)
        if ( saveValidCellNumbers ) then
          ! If valid cell numbers were saved...

          ! Deallocate densityDistribution
          deallocate( densityDistribution ) 

          ! Allocate particles for this IC 
          if(allocated(particleGroups(nic)%Particles)) deallocate(particleGroups(nic)%Particles)
          allocate(particleGroups(nic)%Particles(totalParticleCount))

          ! Assign to the particle group 
          particleGroups(nic)%Mass = effParticleMass 

          ! Create particles
          m = 0
          cellCounter = 0
          do nc=1,validCellCounter

            cellCounter = cellCounter + 1

            ! Skip this cell if all subDivisions remained as zero
            if ( all( subDivisions( cellCounter, : ) .eq. 0 ) ) cycle

            ! Assign as effective mass, the one
            ! that guarantees consistency of total mass
            !particleMass = effParticleMass

            !! For the weird requirement where density
            !! might be negative...
            !if ( densityDistribution(nc) .gt. 0d0 ) then 
            !  particleMass = effParticleMass
            !else ! If zero already cycled 
            !  particleMass = -1*effParticleMass
            !end if

            cellNumber = validCellNumbers(cellCounter)

            ! 0: is for drape. TEMPORARY
            ! Drape = 0: particle placed in the cell. If dry, status to unreleased
            ! Drape = 1: particle placed in the uppermost active cell
            call CreateMassParticlesArray(& 
              particleGroups(nic), cellNumber, m, &
              subDivisions(cellCounter,1),&
              subDivisions(cellCounter,2),&
              subDivisions(cellCounter,3),& 
              0, effParticleMass, particleGroups(nic)%GetReleaseTime(1) )
          end do

        else
          ! Allocate particles without deallocating density
          ! distribution, might be memory demanding for large models

          ! Allocate particles for this IC 
          if(allocated(particleGroups(nic)%Particles)) deallocate(particleGroups(nic)%Particles)
          allocate(particleGroups(nic)%Particles(totalParticleCount))

          ! Assign to the particle group 
          particleGroups(nic)%Mass = effParticleMass 

          ! Create particles
          m = 0
          cellCounter = 0
          do nc=1,grid%CellCount

            ! If no concentration, next
            if ( densityDistribution(nc) .eq. 0d0 ) cycle

            cellCounter = cellCounter + 1

            ! Skip this cell if all subDivisions remained as zero
            if ( all( subDivisions( cellCounter, : ) .eq. 0 ) ) cycle

            ! Assign as effective mass, the one
            ! that guarantees consistency of total mass
            !particleMass = effParticleMass

            !! For the weird requirement where density
            !! might be negative...
            !if ( densityDistribution(nc) .gt. 0d0 ) then 
            !  particleMass = effParticleMass
            !else ! If zero already cycled 
            !  particleMass = -1*effParticleMass
            !end if

            ! 0: is for drape. TEMPORARY
            ! Drape = 0: particle placed in the cell. If dry, status to unreleased
            ! Drape = 1: particle placed in the uppermost active cell
            call CreateMassParticlesArray(& 
              particleGroups(nic), nc, m, &
              subDivisions(cellCounter,1),&
              subDivisions(cellCounter,2),&
              subDivisions(cellCounter,3),& 
              0, effParticleMass, particleGroups(nic)%GetReleaseTime(1) )
          end do

        end if 

        ! Assign for each particle of this group
        !  - id
        !  - group
        !  - sequence number
        !  - layer 
        idmax = 0
        do m = 1, totalParticleCount
          seqNumber = seqNumber + 1
          if(particleGroups(nic)%Particles(m)%Id .gt. idmax) idmax = particleGroups(nic)%Particles(m)%Id
          particleGroups(nic)%Particles(m)%Group = particleGroups(nic)%Group
          particleGroups(nic)%Particles(m)%SequenceNumber = seqNumber
          particleGroups(nic)%Particles(m)%InitialLayer =                                   &
            grid%GetLayer(particleGroups(nic)%Particles(m)%InitialCellNumber)
          particleGroups(nic)%Particles(m)%Layer =                                          &
            grid%GetLayer(particleGroups(nic)%Particles(m)%CellNumber)
        end do

        ! Done with this IC kind

      ! Read initial condition as resident concentration (ML^-3)
      case(1)
        ! Similar to the previous case, but now enforces consistency at a cell level. 
        ! Given a magnitude for particles' mass, determines a number of particles 
        ! inside each cell. From here, establish a per-cell value for the mass, 
        ! obtained from the total mass inside the cell. Particles with non-uniform mass. :O 

        if(allocated(densityDistribution)) deallocate(densityDistribution)
        allocate(densityDistribution(grid%CellCount))

        ! Read particles mass
        read(icUnit, '(a)') line
        icol = 1
        call urword(line,icol,istart,istop,3,n,r,0,0)
        if ( r.lt.0d0 ) then 
          write(outUnit,'(A)') 'Given particle mass is negative. Should be positive. Stop.'
          call ustop('Given particle mass is negative. Should be positive. Stop.')
        end if
        particleMass = r 

        if ( ( this%ParticlesMassOption .eq. 2 ) ) then 
          ! Read solute id
          read(icUnit, '(a)') line
          icol = 1
          call urword(line,icol,istart,istop,2,n,r,0,0)
          if ( n.lt.1 ) then 
            write(outUnit,'(A)')'Given SpeciesID is less than 1. Minimum is 1. Stop.'
            call ustop('Given SpeciesID is less than 1. Minimum is 1. Stop.')
          end if
          soluteId = n
        end if

        ! Read concentration
        if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
          call u3ddblmp(icUnit, outUnit, grid%LayerCount, grid%RowCount,      &
            grid%ColumnCount, grid%CellCount, densityDistribution, aname(1))
        else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
          call u3ddblmpusg(icUnit, outUnit, grid%CellCount, grid%LayerCount,  &
            densityDistribution, aname(1), cellsPerLayer)
        else
          write(outUnit,*) 'Invalid grid type specified when reading IC array ', & 
              particleGroups(nic)%Name, ' name. Stop.'
          call ustop('Invalid grid type specified when reading IC array') 
        end if

        ! Validity of initial condition
        ! Used to allocate subdivisions 
        validCellCounter = count( densityDistribution /= 0d0 )
        if( validCellCounter .eq.0 ) then
          write(outUnit,'(a,a,a)') 'Warning: initial condition ',&
            trim(adjustl(particleGroups(nic)%Name)),' has a distribution only with zeros'
          write(outUnit,'(a)') 'It will not create a particle group. Continue to the next.'
          ! Process the next one
          cycle
        end if 

        ! Allocate subdivisions
        if (allocated(subDivisions)) deallocate(subDivisions)
        allocate( subDivisions(validCellCounter,3) )
        subDivisions(:,:) = 0

        ! If validCellCounter is not equal to the 
        ! size of the initial concentration array, 
        ! and somehow much smaller,  then it could 
        ! be of benefit to save the valid cell ids  
        ! and deallocate densityDistribution 
        saveValidCellNumbers = .false. 
        if ( validCellCounter .lt. 0.5*grid%CellCount ) then
          if ( allocated( validCellNumbers ) ) deallocate( validCellNumbers ) 
          allocate( validCellNumbers(validCellCounter) ) 
          saveValidCellNumbers = .true. 
        end if 

        ! Active cell total mass holder
        if ( allocated( activeCellTotalMass ) ) deallocate( activeCellTotalMass ) 
        allocate( activeCellTotalMass(validCellCounter) ) 
        activeCellTotalMass(:) = 0d0

        ! Loop over densityDistribution and compute:
        !  - totalDissolvedMass
        !  - totalAccumulatedMass: considers retardation factor
        !  - totalParticleCount
        totalDissolvedMass   = 0d0
        totalAccumulatedMass = 0d0
        totalParticleCount   = 0
        cellCounter          = 0
        select case(particleArrayFormat)
        ! for equispaced, quasi-random
        case(0,1)
          do nc = 1, grid%CellCount
            ! If no concentration, next
            if ( densityDistribution(nc) .eq. 0d0 ) cycle

            ! Increase cellCounter
            cellCounter = cellCounter + 1

            ! Compute cell volume
            cellVolume = 1d0
            do nd=1,3
              if( dimensionMask(nd).eq.0 ) cycle
              select case(nd)
              case(1)
                cellVolume = cellVolume*grid%DelX(nc)
              case(2)
                cellVolume = cellVolume*grid%DelY(nc)
              case(3)
                ! simple dZ
                cellVolume = cellVolume*(grid%Top(nc)-grid%Bottom(nc))
              end select
            end do
            cellDissolvedMass = 0d0
            cellTotalMass = 0d0
            ! Absolute value is required for the weird case that 
            ! densityDistribution contains negative values
            cellDissolvedMass = abs(densityDistribution(nc))*porosity(nc)*cellVolume
            totalDissolvedMass = totalDissolvedMass + cellDissolvedMass 
            cellTotalMass = cellDissolvedMass*this%Retardation(nc)
            totalAccumulatedMass = totalAccumulatedMass + cellTotalMass
           
            ! Save cell total mass
            activeCellTotalMass(cellCounter) = cellTotalMass

            ! nParticlesCell: estimate the number of particles 
            ! for the cell using the specified mass 
            nParticlesCell = cellTotalMass/particleMass

            ! If less than 0.5 particle, cycle to the next cell
            ! Note: might not be necessary for this initial condition format 
            !if ( nParticlesCell .lt. nParticlesCellMin ) cycle

            ! Compute shapeFactors only if dimension is active
            ! If not, will remain as zero
            sX = 0
            sY = 0
            sZ = 0
            do nd=1,3
              if( dimensionMask(nd).eq.0 ) cycle
              select case(nd)
              case(1)
                sX = grid%DelX(nc)/(cellVolume**(1d0/dble(nDim)))
              case(2)
                sY = grid%DelY(nc)/(cellVolume**(1d0/dble(nDim)))
              case(3)
                ! simple dZ
                sZ = (grid%Top(nc)-grid%Bottom(nc))/(cellVolume**(1d0/dble(nDim)))
              end select
            end do

            ! Estimate subdivisions
            nPX    = sX*( (nParticlesCell)**(1d0/dble(nDim)) ) 
            nPY    = sY*( (nParticlesCell)**(1d0/dble(nDim)) )
            nPZ    = sZ*( (nParticlesCell)**(1d0/dble(nDim)) )
            iNPX   = max(int( nPX + 0.5d0 ),1)
            iNPY   = max(int( nPY + 0.5d0 ),1)
            iNPZ   = max(int( nPZ + 0.5d0 ),1)
            NPCELL = iNPX*iNPY*iNPZ
            totalParticleCount = totalParticleCount + NPCELL

            ! Save in subdivisions
            subDivisions(cellCounter,1) = iNPX
            subDivisions(cellCounter,2) = iNPY
            subDivisions(cellCounter,3) = iNPZ

            ! Save valid cell numbers
            if ( saveValidCellNumbers ) then
              validCellNumbers(cellCounter) = nc
            end if 
          end do ! end loop over densityDistribution

        ! for random
        case(2)
          do nc = 1, grid%CellCount
            ! If no concentration, next
            if ( densityDistribution(nc) .eq. 0d0 ) cycle

            ! Increase cellCounter
            cellCounter = cellCounter + 1

            ! Compute cell volume
            cellVolume = 1d0
            do nd=1,3
              if( dimensionMask(nd).eq.0 ) cycle
              select case(nd)
              case(1)
                cellVolume = cellVolume*grid%DelX(nc)
              case(2)
                cellVolume = cellVolume*grid%DelY(nc)
              case(3)
                ! simple dZ
                cellVolume = cellVolume*(grid%Top(nc)-grid%Bottom(nc))
              end select
            end do
            cellDissolvedMass = 0d0
            cellTotalMass = 0d0
            ! Absolute value is required for the weird case that 
            ! densityDistribution contains negative values
            cellDissolvedMass = abs(densityDistribution(nc))*porosity(nc)*cellVolume
            totalDissolvedMass = totalDissolvedMass + cellDissolvedMass 
            cellTotalMass = cellDissolvedMass*this%Retardation(nc)
            totalAccumulatedMass = totalAccumulatedMass + cellTotalMass
           
            ! Save cell total mass
            activeCellTotalMass(cellCounter) = cellTotalMass

            ! nParticlesCell: estimate the number of particles 
            ! for the cell using the specified mass 
            nParticlesCell = cellTotalMass/particleMass

            ! If less than 0.5 particle, cycle to the next cell
            ! Note: might not be necessary for this initial condition format
            !if ( nParticlesCell .lt. nParticlesCellMin ) cycle

            ! For random distribution, does not need to 
            ! compute subdivisions, can simply take the integer
            NPCELL = int(nParticlesCell) + 1
            totalParticleCount = totalParticleCount + NPCELL

            ! Save in subdivisions
            ! NOTE: the first position contains
            ! all the particles
            subDivisions(cellCounter,1) = NPCELL
            !subDivisions(cellCounter,2) = 0 
            !subDivisions(cellCounter,3) = 0

            ! Save valid cell numbers
            if ( saveValidCellNumbers ) then
              validCellNumbers(cellCounter) = nc
            end if 
          end do ! end loop over densityDistribution

        end select

        ! Assign totalParticleCount and continue to the next IC if no particles
        particleGroups(nic)%TotalParticleCount = totalParticleCount
        if ( totalParticleCount .eq. 0 ) then 
          write(outUnit,*) ' Warning: initial condition ',&
              particleGroups(nic)%Name,' has zero particles, it will skip this group.'
          ! Process the next one
          cycle
        end if 

        ! The total number of particles is known

        ! Report
        write(outUnit,'(A,es18.9e3)') 'Original particle mass for initial condition = ', particleMass
        write(outUnit,'(A)') 'Particle mass is assigned on a per-cell basis.'
        ! With the number of particles and the total mass, the avg mass can be reported
        avgParticleMass = totalAccumulatedMass/dble(totalParticleCount)
        write(outUnit,'(A,es18.9e3)') 'Average particle mass for initial condition = ', avgParticleMass

        ! As the mass of particles is different, needs to aggregate
        totalAccumulatedParticlesMass = 0d0
        varParticleMass = 0d0
        if ( saveValidCellNumbers ) then
          ! If valid cell numbers were saved...

          ! Deallocate densityDistribution
          deallocate( densityDistribution ) 

          ! Allocate particles for this IC 
          if(allocated(particleGroups(nic)%Particles)) deallocate(particleGroups(nic)%Particles)
          allocate(particleGroups(nic)%Particles(totalParticleCount))

          ! Assign to the particle group 
          !particleGroups(nic)%Mass = effParticleMass 

          ! Create particles
          m = 0
          cellCounter = 0
          do nc=1,validCellCounter

            cellCounter = cellCounter + 1

            ! Skip this cell if all subDivisions remained as zero
            if ( all( subDivisions( cellCounter, : ) .eq. 0 ) ) cycle

            nParticlesCell = dble(product(subDivisions(cellCounter,:),mask=subDivisions(cellCounter,:)>0))
            effParticleMass = activeCellTotalMass(cellCounter)/nParticlesCell
            totalAccumulatedParticlesMass = totalAccumulatedParticlesMass + effParticleMass*nParticlesCell
            varParticleMass = varParticleMass + &
              (nParticlesCell*( effParticleMass - avgParticleMass )**2d0)/dble(totalParticleCount)
            cellNumber = validCellNumbers(cellCounter)

            ! 0: is for drape. TEMPORARY
            ! Drape = 0: particle placed in the cell. If dry, status to unreleased
            ! Drape = 1: particle placed in the uppermost active cell
            call CreateMassParticlesArray(& 
              particleGroups(nic), cellNumber, m, &
              subDivisions(cellCounter,1),&
              subDivisions(cellCounter,2),&
              subDivisions(cellCounter,3),& 
              0, effParticleMass, particleGroups(nic)%GetReleaseTime(1) )
          end do

        else
          ! Allocate particles without deallocating density
          ! distribution, might be memory demanding for large models

          ! Allocate particles for this IC 
          if(allocated(particleGroups(nic)%Particles)) deallocate(particleGroups(nic)%Particles)
          allocate(particleGroups(nic)%Particles(totalParticleCount))

          ! Assign to the particle group 
          !particleGroups(nic)%Mass = effParticleMass 

          ! Create particles
          m = 0
          cellCounter = 0
          do nc=1,grid%CellCount
            ! If no concentration, next
            if ( densityDistribution(nc) .eq. 0d0 ) cycle

            cellCounter = cellCounter + 1

            ! Skip this cell if all subDivisions remained as zero
            if ( all( subDivisions( cellCounter, : ) .eq. 0 ) ) cycle

            nParticlesCell = dble(product(subDivisions(cellCounter,:),mask=subDivisions(cellCounter,:)>0))
            effParticleMass = activeCellTotalMass(cellCounter)/nParticlesCell
            totalAccumulatedParticlesMass = totalAccumulatedParticlesMass + effParticleMass*nParticlesCell
            varParticleMass = varParticleMass + &
              (nParticlesCell*( effParticleMass - avgParticleMass )**2d0)/dble(totalParticleCount)

            ! 0: is for drape. TEMPORARY
            ! Drape = 0: particle placed in the cell. If dry, status to unreleased
            ! Drape = 1: particle placed in the uppermost active cell
            call CreateMassParticlesArray(& 
              particleGroups(nic), nc, m, &
              subDivisions(cellCounter,1),&
              subDivisions(cellCounter,2),&
              subDivisions(cellCounter,3),& 
              0, effParticleMass, particleGroups(nic)%GetReleaseTime(1) )
          end do

        end if 
        ! Report the standard deviation of the particles mass
        write(outUnit,'(A,es18.9e3)') 'Standard deviation of particle mass for initial condition = ', & 
                sqrt( varParticleMass )
        write(outUnit,'(A,es18.9e3)') 'Relative standard deviation of particle mass for initial condition = ', & 
                sqrt( varParticleMass )/avgParticleMass

        ! Assign for each particle of this group
        !  - id
        !  - group
        !  - sequence number
        !  - layer 
        idmax = 0
        do m = 1, totalParticleCount
          seqNumber = seqNumber + 1
          if(particleGroups(nic)%Particles(m)%Id .gt. idmax) idmax = particleGroups(nic)%Particles(m)%Id
          particleGroups(nic)%Particles(m)%Group = particleGroups(nic)%Group
          particleGroups(nic)%Particles(m)%SequenceNumber = seqNumber
          particleGroups(nic)%Particles(m)%InitialLayer = &
            grid%GetLayer(particleGroups(nic)%Particles(m)%InitialCellNumber)
          particleGroups(nic)%Particles(m)%Layer = &
            grid%GetLayer(particleGroups(nic)%Particles(m)%CellNumber)
        end do

        ! Done with this IC kind

      case default
        write(outUnit,*) 'Invalid initial condition kind ', initialConditionFormat, '. Stop.' 
        call ustop('Invalid initial condition kind')
      end select

      ! Report about total mass and number of particles
      ! They should be the same.
      if(this%RetardationFactorOption .gt. 1) then
        write(outUnit,'(A)') 'Retardation factor is considered in total accumulated mass.'
        write(outUnit,'(A,es18.9e3)') 'Total disolved mass for initial condition = ', totalDissolvedMass
        write(outUnit,'(A,es18.9e3)') 'Total accumulated mass for initial condition = ', totalAccumulatedMass
        write(outUnit,'(A,es18.9e3)') 'Total accumulated mass with particles = ', totalAccumulatedParticlesMass
      else 
        write(outUnit,'(A)') 'Retardation factor is unitary so total dissolved mass is the total mass.'
        write(outUnit,'(A,es18.9e3)') 'Total accumulated mass for initial condition = ', totalAccumulatedMass
        write(outUnit,'(A,es18.9e3)') 'Total accumulated mass with particles = ', totalAccumulatedParticlesMass
      end if 
      write(outUnit,'(A,I10)') 'Total number of particles for this initial condition = ', totalParticleCount 

      ! Increment valid counter
      nValidInitialConditions = nValidInitialConditions + 1 

      ! Assign the solute id 
      if ( this%ParticlesMassOption .eq. 2 ) then 
        particleGroups(nic)%Solute = soluteId
      end if 

      ! Incremenent particleCount
      particleCount = particleCount + particleGroups(nic)%TotalParticleCount

      ! Save the current sequence number
      this%currentSeqNumber = seqNumber

    end do ! loop over initial conditions 
    write(outUnit, '(a,i10)') 'Total number of particles on initial conditions = ', particleCount
    write(outUnit, *)

    ! Extend simulationdata to include these particle groups
    if ( nValidInitialConditions .gt. 0 ) then

      ! Update the pgroup count
      newParticleGroupCount = this%ParticleGroupCount + nValidInitialConditions

      ! If no previously defined pgroups, and all the new were 
      ! valid, then simply move_alloc
      if ( ( this%ParticleGroupCount .eq. 0 ) .and. &
           ( nInitialConditions .eq. nValidInitialConditions ) ) then 
        this%ParticleGroupCount = newParticleGroupCount
        this%TotalParticleCount = particleCount
        call move_alloc( particleGroups, this%ParticleGroups )
      else
        ! Needs something more flexible that doesn't 
        ! require copying the pgroups
        allocate(newParticleGroups(newParticleGroupCount))
        ! If some particle groups existed previously
        if( this%ParticleGroupCount .gt. 0 ) then 
          do n = 1, this%ParticleGroupCount
            newParticleGroups(n) = this%ParticleGroups(n)
          end do
        end if 
        pgCount = 0
        do n = 1, nInitialConditions
          if ( particleGroups(n)%TotalParticleCount .eq. 0 ) cycle
          pgCount = pgCount + 1 
          newParticleGroups(pgCount+this%ParticleGroupCount) = particleGroups(n)
        end do 
        if( this%ParticleGroupCount .gt. 0 ) then 
          call move_alloc( newParticleGroups, this%ParticleGroups )
          this%ParticleGroupCount = newParticleGroupCount
          this%TotalParticleCount = this%TotalParticleCount + particleCount
        else
          this%ParticleGroupCount = newParticleGroupCount
          this%TotalParticleCount = particleCount
          call move_alloc( newParticleGroups, this%ParticleGroups )
        end if
      end if 

    end if

    ! flush
    flush(outUnit)

    ! Close ic data file
    close( icUnit )


  end subroutine pr_ReadICData



  ! Read specific SRC data
  subroutine pr_ReadSRCData( this, srcFile, srcUnit, outUnit, grid, flowModelData )
    use UTL8MODULE,only : urword,ustop,u3dintmpusg,u3dintmp,ugetnode
    use linear_interpolation_module, only: linear_interp_1d
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    ! input 
    class(ModpathSimulationDataType), target     :: this
    character(len=200), intent(in)               :: srcFile
    integer, intent(in)                          :: srcUnit
    integer, intent(in)                          :: outUnit
    class(ModflowRectangularGridType),intent(in) :: grid
    class(FlowModelDataType),intent(in)          :: flowModelData
    ! local
    integer :: isThisFileOpen
    integer :: nSources, nValidSources
    integer :: particleCount, totalParticleCount, currentParticleCount
    integer :: nsrc, nSrcBudgets, nsb
    integer, dimension(:), pointer :: dimensionMask
    integer, pointer               :: nDim
    character(len=20) :: srcName
    character(len=20) :: srcSpecKind
    integer :: nAuxNames, naux, nTimes, nt, nCells, nc, nd, nr
    character(len=20),allocatable,dimension(:) :: srcPkgNames
    character(len=20),allocatable,dimension(:) :: auxNames
    integer,allocatable,dimension(:)           :: srcIFaceOpt
    logical :: validAuxNames = .false.
    doubleprecision, allocatable, dimension(:) :: auxMasses
    doubleprecision, allocatable, dimension(:) :: auxEffMasses
    integer, allocatable, dimension(:,:)       :: auxSubDivisions
    integer, dimension(12)                     :: auxFaceDivisions
    doubleprecision, allocatable, dimension(:,:)     :: flowTimeseries     ! nt x ncells
    doubleprecision, allocatable, dimension(:,:,:)   :: auxTimeseries      ! nt x ncells x nauxvars
    doubleprecision, allocatable, dimension(:)       :: timeIntervals      ! nt
    doubleprecision, allocatable, dimension(:)       :: times              ! nt + 1
    integer        , allocatable, dimension(:)       :: srcCellNumbers     ! nCells
    integer        , allocatable, dimension(:)       :: srcCellIFaces      ! nCells
    doubleprecision, allocatable, dimension(:,:,:)   :: cummMassTimeseries ! nt x ncells x nauxvars
    integer        , allocatable, dimension(:,:)     :: auxNPCell
    doubleprecision, allocatable, dimension(:)       :: totMass
    integer        , allocatable, dimension(:)       :: totMassLoc
    doubleprecision, allocatable, dimension(:)       :: nParticlesDbl
    integer        , allocatable, dimension(:)       :: nParticlesInt
    integer        , allocatable, dimension(:)       :: nReleases
    doubleprecision, allocatable, dimension(:)       :: releaseTimes
    doubleprecision, allocatable, dimension(:)       :: cummMassSeries
    doubleprecision :: initialTime, finalTime, minT
    doubleprecision :: effectiveMass, cummEffectiveMass, lastCummMass, totCummEffectiveMass
    doubleprecision :: npAuxDbl
    doubleprecision :: releaseTrackingTime
    integer :: npThisRelease, npNextRelease
    logical :: lastRelease = .false.
    integer :: firstnonzero, nti, nte  
    character(len=20) :: tempChar1
    character(len=20) :: tempChar2
    integer :: m, idmax, seqNumber, cellCounter, offset, cellNumber
    integer :: nValidPGroup = 0
    integer :: newParticleGroupCount, pgCount
    integer :: iFaceNumber, defaultIFaceNumber
    logical :: iFaceOption
    integer :: cellReadFormat
    integer :: layer, row, column
    integer :: nSpecies, ns
    integer :: nTimeIntervals,nMFTimeIntervals, nMFTimes
    integer :: mftCounter, tCounter, itCounter, doCounter, tcount
    doubleprecision, allocatable, dimension(:)     :: mergedTimes
    doubleprecision, allocatable, dimension(:)     :: inputTimes
    doubleprecision, allocatable, dimension(:)     :: mfTimes
    doubleprecision, allocatable, dimension(:)     :: particlesMass
    doubleprecision, allocatable, dimension(:,:)   :: allSpecData         ! nc+2 x nt (column major)
    doubleprecision, allocatable, dimension(:,:,:) :: concTimeseries      ! nt x nc x ns 
    doubleprecision, allocatable, dimension(:,:)   :: flowDataTimeseries  ! nt x ncells
    integer, allocatable, dimension(:)             :: intervalIndex
    integer, allocatable, dimension(:)             :: soluteIds
    integer :: readNTemplates
    integer, allocatable, dimension(:,:) :: nSubDivisions
    integer, allocatable, dimension(:)   :: cellsHolder
    integer, allocatable, dimension(:)   :: cellsPerLayer
    logical :: readCellsFromBudget   = .false.
    logical :: readAsUnstructured = .false.
    logical :: concPerCell = .false.
    logical :: isValid = .false.
    integer :: nColumns
    integer :: ktime, kinitial, kfinal, kdelta
    integer :: correctInterval
    type(ParticleGroupType),dimension(:),allocatable :: particleGroups
    type(ParticleGroupType),dimension(:),allocatable :: newParticleGroups
    ! urword
    character(len=200) :: line
    integer :: icol,istart,istop,n
    doubleprecision :: r
    ! linear interpolation 
    type( linear_interp_1d ) :: interp1d
    integer :: int1dstat
    ! error
    character(len=132) message
    ! u3d
    character(len=24)  :: aname(1)
    data aname(1) /'                   CELLS'/
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW SRC file data'
    write(outUnit, '(1x,a)') '------------------------'

    ! Verify if unit is open 
    isThisFileOpen = -1
    inquire( file=srcFile, number=isThisFileOpen )
    if ( isThisFileOpen .lt. 0 ) then 
      ! No bc file
      write(outUnit,'(A)') 'SRC package was not specified in name file.'
      ! And leave
      return
    end if

    ! Read number of specs
    read(srcUnit, *) nSources
    write(outUnit,'(A,I5)') 'Given number of source specifications = ', nSources

    if(nSources  .le. 0) then
      ! Report
      write(outUnit,'(A)') 'Number of given source specifications is .le. 0. Will not interpret the file .'
      ! And leave
      return
    end if


    ! RW dimensionality vars
    dimensionMask => this%TrackingOptions%dimensionMask
    nDim => this%TrackingOptions%nDim


    nValidSources = 0 
    ! Loop over source specs
    do nsrc = 1, nSources
      
      read(srcUnit, '(a)') line
      icol = 1
      call urword(line,icol,istart,istop,0,n,r,0,0)
      srcName = line(istart:istop)

      ! Report which SRC specification will be processed
      write(outUnit,'(A,A)') 'Processing SRC specification: ', trim(adjustl(srcName))

      ! Interpret specification kind
      read(srcUnit, '(a)') line
      icol = 1
      call urword(line,icol,istart,istop,1,n,r,0,0)
      srcSpecKind = line(istart:istop)

      select case (srcSpecKind) 
      ! Concentrations are read from AUX variables stored in the budget file.
      ! Flow-rates and times are extracted from the budget for the characteristic
      ! simulation times ( reftime, stoptime ).
      case ('AUX','AUXILIARY')

        ! Process reading auxiliary variables
        write(outUnit,'(A,I5)') 'SRC specification will be read from AUXILIARY variables.'

        read(srcUnit, '(a)') line
        icol = 1
        call urword(line,icol,istart,istop,2,n,r,0,0)
        nSrcBudgets = n
        
        if ( nSrcBudgets .lt. 1 ) then
          write(outUnit,'(A,A,A)') 'Number of source budgets for spec ', srcName ,' is .lt. 1. It should be at least 1.'
          call ustop('Number of source budgets is .lt. 1. It should be at least 1. Stop.')
        end if 

        ! Interpret source budgets
        if( allocated( srcPkgNames ) ) deallocate( srcPkgNames ) 
        allocate( srcPkgNames( nSrcBudgets ) ) 
        if( allocated( srcIFaceOpt ) ) deallocate( srcIFaceOpt ) 
        allocate( srcIFaceOpt( nSrcBudgets ) )
        srcIFaceOpt(:) = 0 
        do nsb=1,nSrcBudgets

          ! Read the pkg/budget header
          read(srcUnit, '(a)') line
          icol = 1
          call urword(line,icol,istart,istop,1,n,r,0,0)
          srcPkgNames(nsb) = line(istart:istop)

          ! Read iFaceOption
          call urword(line,icol,istart,istop,2,n,r,0,0)
          srcIFaceOpt(nsb) = n 

          ! Report
          write(outUnit,'(A,A)') 'Source budget name: ', trim(adjustl(srcPkgNames(nsb)))

          ! Activate/deactivate iFaceOption for this source
          iFaceOption = .false.
          if ( srcIFaceOpt(nsb) .gt. 0 ) then 
            iFaceOption = .true.
            write(outUnit,'(A,A)') 'Will interpret IFACE from the budget file.'
          end if  

          ! Number of aux variables and allocate
          read(srcUnit, '(a)') line
          icol = 1
          call urword(line,icol,istart,istop,2,n,r,0,0)
          nAuxNames = n 
          if ( nAuxNames .lt. 1 ) then
            write(outUnit,'(A,A,A)') 'Number of aux variables source ', trim(adjustl(srcPkgNames(nsb))) ,' should be at least 1.'
            call ustop('Number of aux variables is .lt. 1. It should be at least 1. Stop.')
          end if 
          if ( allocated( auxNames ) ) deallocate( auxNames ) 
          allocate( auxNames( nAuxNames ) )
          if ( allocated( auxMasses ) ) deallocate( auxMasses ) 
          allocate( auxMasses( nAuxNames ) )
          if ( allocated( auxSubDivisions ) ) deallocate( auxSubDivisions ) 
          allocate( auxSubDivisions( nAuxNames,3 ) )

          ! Read the solute ids if the simulation demands  
          if ( ( this%ParticlesMassOption .eq. 2 ) ) then 
            write(outUnit,'(A)') 'Will read species ids due to ParticlesMassOption.eq.2 '
            ! Read solute id
            if( allocated( soluteIds ) ) deallocate( soluteIds )
            allocate( soluteIds( nAuxNames ) )
          end if

          ! Loop over aux names and interpret data
          ! Requires some health checks
          do naux = 1, nAuxNames
            ! read the complete line and then extract specific params
            read(srcUnit, '(a)') line

            ! auxName
            icol = 1
            call urword(line,icol,istart,istop,1,n,r,0,0)
            auxNames(naux) = line(istart:istop)

            ! particles mass
            call urword(line,icol,istart,istop,3,n,r,0,0)
            auxMasses(naux) = r
            
            ! template
            call urword(line,icol,istart,istop,2,n,r,0,0)
            auxSubDivisions(naux,1) = n 
            call urword(line,icol,istart,istop,2,n,r,0,0)
            auxSubDivisions(naux,2) = n 
            call urword(line,icol,istart,istop,2,n,r,0,0)
            auxSubDivisions(naux,3) = n 
            ! Validate template

            ! Read solute id depending on solutes option
            if ( this%ParticlesMassOption .eq. 2 ) then 
              call urword(line,icol,istart,istop,2,n,r,0,0)
              ! If the soluteid is not given, stop, it is required by this option
              if ( n.lt. 1) then 
               write(outUnit,'(A)') 'SpeciesID for source is invalid. It was not given or is less than 1. Stop.'
               call ustop('SpeciesID for source is invalid. It was not given or is less than 1. Stop.')
              end if 
              soluteIds(naux) = n 
            end if

          end do ! naux = 1, nAuxNames 
        
          ! Until this point, necessary data for reading auxiliary 
          ! variables and transforming into particles is available
          ! for this source budget

          ! Validate given aux names
          validAuxNames = flowModelData%ValidateAuxVarNames( srcPkgNames( nsb ), auxNames, this%isMF6, iFaceOption )
          if ( .not. validAuxNames ) then 
            write(outUnit,'(A,A,A)') 'Not all aux variables were found in source ', trim(adjustl(srcPkgNames(nsb))), '.'
            call ustop('Not all aux variables were found in source or it does not support aux vars. Stop.')
          end if 

          ! While reading from AUX vars, uses simulation characteristic times
          ! The initial MODFLOW time is always ReferenceTime
          initialTime = this%ReferenceTime  
          ! The final MODFLOW time...
          ! StopTime was updated at MPathRW, so it already 
          ! contains the logic of StoppingTimeOption
          ! In cases that stoptime is infinite, force the MODFLOW limits 
          if ( this%TrackingOptions%BackwardTracking ) then 
            finalTime = this%ReferenceTime - this%StopTime ! Could it be less than zero ? 
            if ( finalTime .lt. 0d0 ) finalTime = 0d0
          else 
            finalTime = this%ReferenceTime + this%StopTime
            if ( finalTime.gt.this%tdisData%TotalTimes(size(this%tdisData%TotalTimes)) ) then 
              finalTime = this%tdisData%TotalTimes(size(this%tdisData%TotalTimes)) 
            end if
          end if 

          ! Obtain flow and aux vars timeseries.
          ! Function allocates and return necessary arrays
          ! Note: times is a vector built from tdisData%TotalTimes, 
          ! including intial and final times. timeIntervals is computed as diff(times)
          call flowModelData%LoadFlowAndAuxTimeseries( srcPkgNames( nsb ), auxNames,& 
                                  this%isMF6, initialTime, finalTime, this%tdisData,&
                                flowTimeseries, auxTimeseries, timeIntervals, times,&
                                srcCellNumbers, outUnit, iFaceOption, srcCellIFaces,&
                                              this%TrackingOptions%BackwardTracking )

          ! From here until the creation of particles the process 
          ! is the same than for the SPEC format so it could be considered 
          ! to wrap these ops on a common function 
          ! Now compute the cummulative mass function to obtain the total injected mass
          if ( allocated( cummMassTimeseries ) ) deallocate( cummMassTimeseries ) 
          allocate( cummMassTimeseries, mold=auxTimeseries ) 
          cummMassTimeseries(:,:,:) = 0d0

          ! Integrate for each aux var
          ! Considers STEPWISE quantities
          nTimes = size(timeIntervals)
          do naux=1, nAuxNames
            do nt=1,nTimes
              if ( nt.gt.1 ) then
                cummMassTimeseries(nt,:,naux) = &
                  cummMassTimeseries(nt-1,:,naux) + flowTimeseries(nt,:)*auxTimeseries(nt,:,naux)*timeIntervals(nt)
                cycle
              end if
              if ( nt.eq.1 ) cummMassTimeseries(nt,:,naux) = flowTimeseries(nt,:)*auxTimeseries(nt,:,naux)*timeIntervals(nt)
            end do
          end do

          ! It needs to transform the cummulative mass into particles
          nCells = size(srcCellNumbers)
          if ( allocated(totMass) ) deallocate(totMass)
          allocate( totMass( nCells ) )
          if ( allocated(totMassLoc) ) deallocate(totMassLoc)
          allocate( totMassLoc( nCells ) )
          if ( allocated(nParticlesDbl) ) deallocate(nParticlesDbl)
          allocate( nParticlesDbl( nCells ) )
          if ( allocated(nParticlesInt) ) deallocate(nParticlesInt)
          allocate( nParticlesInt( nCells ) )
          if ( allocated(nReleases) ) deallocate(nReleases)
          allocate( nReleases( nCells ) )
          if ( allocated( auxEffMasses ) ) deallocate( auxEffMasses ) 
          allocate( auxEffMasses, mold=auxMasses )
          if ( allocated( auxNPCell ) ) deallocate( auxNPCell ) 
          allocate( auxNPCell( nAuxNames, nCells ) )
          ! Allocate srcCellIFaces as a placeholder with zeros
          ! if it was not allocated at the extraction function
          if ( .not. allocated(srcCellIFaces) ) then
            allocate( srcCellIFaces(nCells) )
            srcCellIFaces(:) = 0
          end if 

          ! For this source budget, it should create 
          ! particlegroups for what specifically ? 

          ! One for each aux var  ? MAYBE
          ! One per source budget ? COULD ALSO BE
          ! One for each cell     ? NO

          ! For multidispersion or multispecies simulations,
          ! the current program structure relates a solute with 
          ! particle group ids. 

          ! In the most likely scenario that each aux var represents 
          ! a different solute, then it makes sense to consider 
          ! a particle group per aux var. Also, the package structure
          ! will request a soluteid for each auxvar in case the simulation 
          ! is multispecies.

          ! Carrier for candidate particle groups 
          if ( allocated(particleGroups) ) deallocate( particleGroups )
          allocate(particleGroups(nAuxNames))

          ! Intialize some counters
          nValidPGroup = 0
          particleCount = 0

          ! Do it for each aux var
          do naux=1,nAuxNames

            ! Report
            write(outUnit,'(A,A)') 'Processing auxiliary variable: ', trim(adjustl(auxNames(naux)))

            ! Give a name to the particle group
            ! SRCname_idsrcbudget_idauxvar
            tempChar1 = ''
            tempChar2 = ''
            write( unit=tempChar1, fmt=* )nsb 
            write( unit=tempChar2, fmt=* )naux 
            particleGroups(naux)%Name = trim(adjustl(srcName))//'_'//trim(adjustl(tempChar1))//'_'//trim(adjustl(tempChar2))

            ! Increase pgroup counter
            particleGroups(naux)%Group = this%ParticleGroupCount + naux

            ! Assign soluteIds if simulation demands
            if ( this%ParticlesMassOption .eq. 2 ) then 
              particleGroups(naux)%Solute = soluteIds(naux)
            end if

            ! Correction to the particle template based on dimension mask
            ! If dimension is active, leave the given value. 
            ! If not, switch it to one
            do nd = 1, 3
              if ( dimensionMask(nd) .eq. 1 ) cycle 
              auxSubDivisions(naux,nd) = 1
            end do
                
            ! Particles per cell
            auxNPCell(naux,:) = product(auxSubDivisions(naux,:))

            ! Interpret srcCellIFaces and apply
            ! correction to particles per cell.
            if ( iFaceOption ) then
              do nc = 1, nCells
                auxFaceDivisions(:) = 0
                select case(srcCellIFaces(nc))
                case(1,2)
                  auxFaceDivisions(2*srcCellIFaces(nc)-1) = auxSubDivisions(naux,3) ! nver
                  auxFaceDivisions(2*srcCellIFaces(nc))   = auxSubDivisions(naux,2) ! nrow
                case(3,4)
                  auxFaceDivisions(2*srcCellIFaces(nc)-1) = auxSubDivisions(naux,3) ! nver
                  auxFaceDivisions(2*srcCellIFaces(nc))   = auxSubDivisions(naux,1) ! ncol
                case(5,6)
                  auxFaceDivisions(2*srcCellIFaces(nc)-1) = auxSubDivisions(naux,2) ! nrow
                  auxFaceDivisions(2*srcCellIFaces(nc))   = auxSubDivisions(naux,1) ! ncol
                case default
                  cycle
                end select
                ! Update NP cell
                auxNPCell(naux,nc) = auxFaceDivisions(2*srcCellIFaces(nc)-1)*auxFaceDivisions(2*srcCellIFaces(nc))
              end do
            end if

            ! Max cumm mass stats
            totMass    = maxval(cummMassTimeseries(:,:,naux), dim=1) ! Max cumm mass in time
            totMassLoc = maxloc(cummMassTimeseries(:,:,naux), dim=1) ! Loc of max cumm mass in time
            ! Note: the location of max mass corresponds to the end of the interval at index loc

            ! With the given particles mass, estimate the number of particles
            ! necessary for achieving the cummulative mass 
            nParticlesDbl = totMass/auxMasses(naux)
            ! ... and using the number of particles per cell, estimate
            ! the number of releases per cell.
            nReleases     = int(nParticlesDbl/auxNPCell(naux,:)+0.5)

            ! Up to this point, the number of particles 
            ! needed can be computed as nReleases*NPCELL

            ! Cycle to the next aux var if the total 
            ! number of particles is zero.
            totalParticleCount = sum(nReleases*auxNPCell(naux,:))
            if ( totalParticleCount .lt. 1 ) then 
              write(outUnit,*) ' Warning: pgroup ',&
                trim(adjustl(particleGroups(naux)%Name)),' has zero particles, it will skip this aux var.'
              ! Process the next aux variable
              cycle
            end if 

            ! Allocate particles 
            particleGroups(naux)%TotalParticleCount = totalParticleCount
            if(allocated(particleGroups(naux)%Particles)) deallocate(particleGroups(naux)%Particles)
            allocate(particleGroups(naux)%Particles(totalParticleCount))

            ! Variables for defining the particles of this group
            idmax = 0
            offset = 0
            seqNumber = this%currentSeqNumber
            cellCounter = 0
            currentParticleCount = 0 

            ! Given the number of releases, interpolate release times
            ! for each cell, using as source the cummulative mass function
            totCummEffectiveMass = 0d0
            do nc=1, nCells

              if ( nReleases(nc) .lt. 1 ) cycle ! At least one release

              ! DEPRECATE
              cellCounter = cellCounter + 1
              if ( allocated( releaseTimes ) ) deallocate( releaseTimes ) 
              allocate( releaseTimes( nReleases(nc) ) )
              ! END DEPRECATE

              ! The effective mass of particles for the aux var of this cell
              effectiveMass = totMass(nc)/(nReleases(nc)*auxNPCell(naux,nc))

              ! Create particle template for this cell
              ! 0: is for drape. (TEMP)
              ! Drape = 0: particle placed in the cell. If dry, status to unreleased
              ! Drape = 1: particle placed in the uppermost active cell
              ! -999: is a placeholder for release time, it will assigned later

              ! Interpret srcCellIFaces
              ! This check requires allocated srcCellIFaces, it could be improved. 
              if ( iFaceOption .and. (srcCellIFaces(nc).gt.0)  ) then
                ! Create particles at cell face
                auxFaceDivisions(:) = 0
                select case(srcCellIFaces(nc))
                case(1,2)
                  auxFaceDivisions(2*srcCellIFaces(nc)-1) = auxSubDivisions(naux,3) ! nver
                  auxFaceDivisions(2*srcCellIFaces(nc))   = auxSubDivisions(naux,2) ! nrow
                case(3,4)
                  auxFaceDivisions(2*srcCellIFaces(nc)-1) = auxSubDivisions(naux,3) ! nver
                  auxFaceDivisions(2*srcCellIFaces(nc))   = auxSubDivisions(naux,1) ! ncol
                case(5,6)
                  auxFaceDivisions(2*srcCellIFaces(nc)-1) = auxSubDivisions(naux,2) ! nrow
                  auxFaceDivisions(2*srcCellIFaces(nc))   = auxSubDivisions(naux,1) ! ncol
                end select
                call CreateMassParticlesOnFaces(&
                  particleGroups(naux),    &
                  srcCellNumbers(nc),      &
                  currentParticleCount,    &
                  auxFaceDivisions,        &
                  0, effectiveMass, -999d0 )
              else
                ! Normal internal array of particles
                call CreateMassParticlesAsInternalArray(& 
                  particleGroups(naux),     &
                  srcCellNumbers(nc),       &
                  currentParticleCount,     &
                  auxSubDivisions(naux,1),  &
                  auxSubDivisions(naux,2),  &
                  auxSubDivisions(naux,3),  & 
                  0, effectiveMass, -999d0  )
              end if

              ! Assign values for particle template
              !   - layer
              !   - group
              do m = 1, currentParticleCount 
                particleGroups(naux)%Particles(m)%Group = particleGroups(naux)%Group
                particleGroups(naux)%Particles(m)%InitialLayer =                       &
                  grid%GetLayer(particleGroups(naux)%Particles(m)%InitialCellNumber)
                particleGroups(naux)%Particles(m)%Layer =                              &
                  grid%GetLayer(particleGroups(naux)%Particles(m)%CellNumber)
              end do

              ! Note: 
              ! Interpolator requires strictly increasing x and 
              ! the cummulative mass function may have flat sections.

              ! Array for the cumm series considering a zero at the beginning
              if ( allocated( cummMassSeries ) ) deallocate( cummMassSeries ) 
              allocate( cummMassSeries(totMassLoc(nc)+1) )
              cummMassSeries(:)  = 0d0
              cummMassSeries(2:) = cummMassTimeseries(1:totMassLoc(nc),nc,naux)

              ! nti, nte: time indexes to be passed to interpolator as source
              nti = 0
              nte = 0
              firstnonzero = 0
              lastCummMass = 0d0
              cummEffectiveMass = 0d0
              npNextRelease = 0

              ! Loop over cummulative mass series
              do nt = 1, totMassLoc(nc)+1

                ! Advance loop until finding the first non-zero
                if (&
                  ( cummMassSeries( nt ) .eq. 0d0 ) .and. &
                  ( firstnonzero .eq. 0 ) ) cycle

                ! If found a nonzero mass, save the first nonzero
                if ( firstnonzero .eq. 0 ) then 
                  nti = nt-1
                  if( nt .eq. 1 ) nti = 1 ! just in case, but it shouldn't
                  firstnonzero = nt
                  lastCummMass = cummMassSeries(nt)
                  if( nt .ne. (totMassLoc(nc)+1) ) cycle 
                end if

                ! If the current mass is higher than the last found, 
                ! and the initial index is defined, then everything is ok, continue
                ! to the next index if it is not the last.
                if ( (cummMassSeries(nt) .gt. lastCummMass) .and. (nti.gt.0) ) then 
                  lastCummMass = cummMassSeries(nt)
                  if( nt.ne.(totMassLoc(nc)+1) ) cycle
                else if ( (cummMassSeries(nt) .eq. lastCummMass) .and. (nti.eq.0) ) then
                  ! Still on a flat area
                  cycle
                else if ( (cummMassSeries(nt) .gt. lastCummMass) .and. (nti.eq.0) ) then
                  ! Define the new starting point and update lastCummMass
                  nti = nt - 1
                  lastCummMass = cummMassSeries(nt)
                  if( nt.ne.(totMassLoc(nc)+1) ) cycle
                end if

                ! If it is equal to the last, and the first 
                ! index is defined, then we are on a flat zone
                if ( (cummMassSeries(nt) .eq. lastCummMass) .and. (nti.gt.0) ) then
                  ! This is the last index for interpolation 
                  nte = nt-1 
                  ! If it is the last loop, set at the last index in array
                  if ( nt.eq.(totMassLoc(nc)+1) ) nte = nt

                  ! Can it occur ?
                  if ( nte .eq. nti ) then 
                    write(outUnit,'(A)') 'Range of time indexes for 1d interpolation is inconsistent.'
                    call ustop('Range of time indexes for 1d interpolation is inconsistent. Stop.')
                  end if

                  ! Initialize interpolator
                  call interp1d%initialize(cummMassSeries(nti:nte),times(nti:nte), int1dstat)
                  if ( int1dstat .gt. 0 ) then 
                    write(outUnit,'(A)') 'There was a problem while initializing 1d interpolator.'
                    call ustop('There was a problem while initializing 1d interpolator. Stop.')
                  end if

                  ! Initialize release variables
                  lastRelease = .false.
                  npThisRelease = auxNPCell(naux,nc)
                  if ( npNextRelease .gt. 0 ) then 
                    npThisRelease = npNextRelease
                    npNextRelease = 0
                  end if 

                  ! And loop until a maximum of nReleases
                  ! It will break based on cummulative mass
                  do nr=1, nReleases(nc)

                    ! Increase the cummulative mass counter
                    cummEffectiveMass = cummEffectiveMass + npThisRelease*effectiveMass

                    if ( cummEffectiveMass .ge. lastCummMass ) then
                      ! If above the current limit, add some particles
                      ! as long as the limit remains below lastCummMass 
                      ! and correct cummEffectiveMass
                      cummEffectiveMass = cummEffectiveMass - npThisRelease*effectiveMass
                      npAuxDbl = (lastCummMass - cummEffectiveMass)/effectiveMass
                      npThisRelease = int(npAuxDbl+0.5)

                      ! If no particles, done
                      if ( npThisRelease .lt. 1 ) exit

                      ! This variable keeps track of the particles 
                      ! needed to begin releases next time and keep 
                      ! consistency with total mass for the computed
                      ! effectiveMass
                      npNextRelease = auxNPCell(naux,nc) - npThisRelease
                      cummEffectiveMass = cummEffectiveMass + npThisRelease*effectiveMass
                      lastRelease = .true.
                    end if

                    ! Interpolate the release time        
                    call interp1d%evaluate( cummEffectiveMass, releaseTimes(nr) ) ! Replace array of release times by a single value

                    ! Fix the release time considering backward/forward tracking
                    if ( this%TrackingOptions%BackwardTracking ) then
                      ! The interpolated value is the time in the context of the MODFLOW model, 
                      ! if backward tracking, this time is decreasing with respect to reference time.
                      releaseTrackingTime = this%ReferenceTime - releaseTimes(nr)
                    else
                      ! In forward tracking, the interpolated time is increasing with respect to 
                      ! reference time 
                      releaseTrackingTime = releaseTimes(nr) - this%ReferenceTime
                    end if 

                    ! Assign to particles
                    do m = 1, npThisRelease
                     idmax = idmax + 1
                     seqNumber = seqNumber + 1
                     particleGroups(naux)%Particles(offset+m)%Id = idmax
                     particleGroups(naux)%Particles(offset+m)%SequenceNumber = seqNumber
                     particleGroups(naux)%Particles(offset+m)%Group = particleGroups(naux)%Group
                     particleGroups(naux)%Particles(offset+m)%Drape = particleGroups(naux)%Particles(m)%Drape
                     particleGroups(naux)%Particles(offset+m)%Status = particleGroups(naux)%Particles(m)%Status
                     particleGroups(naux)%Particles(offset+m)%InitialCellNumber = & 
                             particleGroups(naux)%Particles(m)%InitialCellNumber
                     particleGroups(naux)%Particles(offset+m)%InitialLayer = particleGroups(naux)%Particles(m)%InitialLayer
                     particleGroups(naux)%Particles(offset+m)%InitialFace = particleGroups(naux)%Particles(m)%InitialFace
                     particleGroups(naux)%Particles(offset+m)%InitialLocalX = particleGroups(naux)%Particles(m)%InitialLocalX
                     particleGroups(naux)%Particles(offset+m)%InitialLocalY = particleGroups(naux)%Particles(m)%InitialLocalY
                     particleGroups(naux)%Particles(offset+m)%InitialLocalZ = particleGroups(naux)%Particles(m)%InitialLocalZ
                     particleGroups(naux)%Particles(offset+m)%InitialTrackingTime = releaseTrackingTime
                     particleGroups(naux)%Particles(offset+m)%TrackingTime = & 
                             particleGroups(naux)%Particles(offset+m)%InitialTrackingTime
                     particleGroups(naux)%Particles(offset+m)%CellNumber = particleGroups(naux)%Particles(m)%CellNumber
                     particleGroups(naux)%Particles(offset+m)%Layer = particleGroups(naux)%Particles(m)%Layer
                     particleGroups(naux)%Particles(offset+m)%Face = particleGroups(naux)%Particles(m)%Face
                     particleGroups(naux)%Particles(offset+m)%LocalX = particleGroups(naux)%Particles(m)%LocalX
                     particleGroups(naux)%Particles(offset+m)%LocalY = particleGroups(naux)%Particles(m)%LocalY
                     particleGroups(naux)%Particles(offset+m)%LocalZ = particleGroups(naux)%Particles(m)%LocalZ
                     ! Mass ( potentially different effectiveMass for different cells )
                     particleGroups(naux)%Particles(offset+m)%Mass   = effectiveMass
                    end do

                    ! Increase offset
                    offset = offset + npThisRelease
                    
                    ! Restart npThisRelease
                    npThisRelease = auxNPCell(naux,nc)

                    if ( lastRelease ) then
                      ! And exit loop over releases
                      exit
                    end if 

                  end do ! nr=1, nReleases(nc)

                  ! After interpolating release times, destroy interpolator
                  call interp1d%destroy()

                  ! It should reset nti, nte
                  nti = 0
                  nte = 0

                end if

              end do ! nt = 1, totMassLoc(nc)+1

              ! Save the current particle count, it will add more
              ! corresponding to the next cell
              currentParticleCount = offset

              ! Accumulate over all cells for report 
              totCummEffectiveMass = totCummEffectiveMass + cummEffectiveMass

            end do ! nc=1, nCells

            ! Save the current sequence number once finished the 
            ! processing of this aux var
            this%currentSeqNumber = seqNumber

            ! Increase counters
            if ( particleGroups(naux)%TotalParticleCount .gt. 0 ) then 
              nValidPGroup = nValidPGroup + 1
              particleCount =  particleCount + particleGroups(naux)%TotalParticleCount 
              write(outUnit,'(A,es18.9e3)') 'Total released mass related to aux var  = ', totCummEffectiveMass
              write(outUnit,'(A,I10)') 'Total number of particles related to aux var = ', particleGroups(naux)%TotalParticleCount 
            end if


          end do ! naux = 1, nAuxNames
       
          ! It needs to add the newly created groups to simulationData%ParticleGroups

          ! Extend simulationdata to include these particle groups
          if ( nValidPGroup .gt. 0 ) then 
            newParticleGroupCount = this%ParticleGroupCount + nValidPGroup
            allocate(newParticleGroups(newParticleGroupCount))
            ! If some particle groups existed previously
            if( this%ParticleGroupCount .gt. 0 ) then 
              do n = 1, this%ParticleGroupCount
                newParticleGroups(n) = this%ParticleGroups(n)
              end do
            end if 
            pgCount = 0
            do n = 1, nAuxNames
              if ( particleGroups(n)%TotalParticleCount .eq. 0 ) cycle
              pgCount = pgCount + 1 
              newParticleGroups(pgCount+this%ParticleGroupCount) = particleGroups(n)
            end do 
            if( this%ParticleGroupCount .gt. 0 ) then 
              call move_alloc( newParticleGroups, this%ParticleGroups )
              this%ParticleGroupCount = newParticleGroupCount
              this%TotalParticleCount = this%TotalParticleCount + particleCount
            else
              this%ParticleGroupCount = newParticleGroupCount
              this%TotalParticleCount = particleCount
              call move_alloc( newParticleGroups, this%ParticleGroups )
            end if
          end if


        end do ! nsb=1,nSrcBudgets


      ! Concentrations, injection times and cells are specified by the user.
      ! Flow-rates are extracted from a given BUDGET header
      case ('SPEC','SPECIFIED')

        ! Report 
        write(outUnit,'(A,I5)') 'SRC specification will be read with the SPECIFIED format.'

        read(srcUnit, '(a)') line
        icol = 1
        call urword(line,icol,istart,istop,2,n,r,0,0)
        nSrcBudgets = n
        
        if ( nSrcBudgets .lt. 1 ) then
          write(outUnit,'(A,A,A)') 'Number of source budgets for spec ', srcName ,' is .lt. 1. It should be at least 1.'
          call ustop('Number of source budgets is .lt. 1. It should be at least 1. Stop.')
        end if 

        ! Interpret source budgets
        if( allocated( srcPkgNames ) ) deallocate( srcPkgNames ) 
        allocate( srcPkgNames( nSrcBudgets ) )
        if( allocated( srcIFaceOpt ) ) deallocate( srcIFaceOpt ) 
        allocate( srcIFaceOpt( nSrcBudgets ) )
        srcIFaceOpt(:) = 0 
        do nsb=1,nSrcBudgets

          ! Read the budget header
          read(srcUnit, '(a)') line
          icol = 1
          call urword(line,icol,istart,istop,1,n,r,0,0)
          srcPkgNames(nsb) = line(istart:istop)

          ! Report
          write(outUnit,'(A,A)') 'Source budget header: ', trim(adjustl(srcPkgNames(nsb)))

          ! Read iFaceOption
          call urword(line,icol,istart,istop,2,n,r,0,0)
          srcIFaceOpt(nsb) = n 

          ! Read defaultIFaceNumber
          call urword(line,icol,istart,istop,2,n,r,0,0)
          defaultIFaceNumber = n 

          ! Activate/deactivate iFaceOption for this source
          iFaceOption = .false.
          if ( srcIFaceOpt(nsb) .gt. 0 ) then 
            iFaceOption = .true.
          end if  

          ! Read cell format 
          read(srcUnit, '(a)') line
          icol = 1
          call urword(line,icol,istart,istop,2,n,r,0,0)
          cellReadFormat = n

          ! Read cells 
          readCellsFromBudget = .false.
          select case(cellReadFormat)
          ! From budget file
          case(0)
            write(outUnit,'(A)') 'Cells will be taken from the budget header.'
            readCellsFromBudget = .true.
            ! nCells is determined after extracting flow-rates
            if ( iFaceOption ) then 
              write(outUnit,'(A)') 'IFACE option is activated for this source.'
            end if 
          ! As list of given cells
          case(1)
            write(outUnit,'(A)') 'Cells are expected to be specified as a list.'
            if ( iFaceOption ) then 
              write(outUnit,'(A)') 'Will read IFACE next to the cell ids.'
            end if 
            ! Number of cells
            nCells = 0
            call urword(line,icol,istart,istop,2,n,r,0,0)
            nCells = n
            if ( nCells .lt. 1 ) then 
              write(outUnit,'(A)') 'Number of specified cells is invalid. It should be at least 1.'
              call ustop('Number of specified cells is invalid. It should be at least 1. Stop.')
            end if 
            ! Cell number format
            readAsUnstructured = .false.
            call urword(line,icol,istart,istop,2,n,r,0,0)
            if ( n.gt.0 ) then 
              readAsUnstructured = .true.
            end if
            ! Different concentration per cell
            concPerCell = .false.
            call urword(line,icol,istart,istop,2,n,r,0,0)
            if ( n.gt.0 ) then
              concPerCell = .true. 
              write(outUnit,'(A)') 'Will consider different concentration per cell. Expects nSpecies*nCells columns.'
            end if  
            ! Allocate array for cell ids
            if ( allocated( srcCellNumbers ) ) deallocate( srcCellNumbers )
            allocate( srcCellNumbers(nCells) )
            ! iFaceOption
            if ( iFaceOption ) then 
              if ( allocated(srcCellIFaces) ) deallocate( srcCellIFaces )
              allocate( srcCellIFaces(nCells) )
            end if
            ! Read the cells
            if( readAsUnstructured ) then
              if ( .not. iFaceOption ) then 
                do nc = 1, nCells
                  read(srcUnit,*) cellNumber
                  srcCellNumbers(nc) = cellNumber
                end do
              else
                ! Read as cellNumber, iFaceNumber
                ! It uses urword to allow for undefined iFaceNumber
                do nc = 1, nCells
                  iFaceNumber = 0
                  read(srcUnit, '(a)') line
                  icol = 1
                  call urword(line,icol,istart,istop,2,n,r,0,0)
                  cellNumber = n 
                  call urword(line,icol,istart,istop,2,n,r,0,0)
                  iFaceNumber = n 
                  if ( n.lt.1 ) iFaceNumber = defaultIFaceNumber
                  srcCellNumbers(nc) = cellNumber
                  srcCellIFaces(nc)  = iFaceNumber
                end do
              end if  
            else
              ! Validate grid type for reading as layer, row, column
              ! GridType .eq. 1 = RectangularGridDis       (MF2005)
              ! GridType .eq. 2 = RectangularGridDisuMfusg (MFUSG)
              ! GridType .eq. 3 = RectangularGridDisMf6    (MF6DIS)
              ! GridType .eq. 4 = RectangularGridDisvMf6   (MF6DISV)
              ! GridType .eq. 5 = RectangularGridDisuMf6   (MF6DISU)
              select case(grid%GridType)
              case(2,4,5) ! mfusg, disv, disu
              write(outUnit,'(A)') 'Error: reading cells as (lay,row,col) only allowed for structured grids (DIS).'
              call ustop('Error: reading cells as (lay,row,col) only allowed for structured grids (DIS).')
              end select

              if ( .not. iFaceOption ) then 
                ! Read as layer, row, column
                do nc = 1, nCells
                  read(srcUnit, *) layer, row, column
                  call ugetnode(&
                    grid%LayerCount,   &
                    grid%RowCount,     &
                    grid%ColumnCount,  &
                    layer, row, column,&
                    cellNumber )
                  srcCellNumbers(nc) = cellNumber
                end do 
              else
                ! Read as layer, row, column iFaceNumber
                ! It uses urword to allow for undefined iFaceNumber
                do nc = 1, nCells
                  iFaceNumber = 0
                  read(srcUnit, '(a)') line
                  icol = 1
                  call urword(line,icol,istart,istop,2,n,r,0,0)
                  layer = n 
                  call urword(line,icol,istart,istop,2,n,r,0,0)
                  row = n 
                  call urword(line,icol,istart,istop,2,n,r,0,0)
                  column = n 
                  call urword(line,icol,istart,istop,2,n,r,0,0)
                  iFaceNumber = n 
                  if ( n.lt.1 ) iFaceNumber = defaultIFaceNumber
                  call ugetnode(&
                    grid%LayerCount,   &
                    grid%RowCount,     &
                    grid%ColumnCount,  &
                    layer, row, column,&
                    cellNumber )
                  srcCellNumbers(nc) = cellNumber
                  srcCellIFaces(nc)  = iFaceNumber
                end do 
              end if
            end if
          ! As 3d array
          case(2)
            write(outUnit,'(A)') 'Cells are expected to be specified as array.'

            if ( iFaceOption ) then 
              write(outUnit,'(A)') 'IFACE option is activated for this source.'
            end if 

            ! Required for u3d
            if(allocated(cellsHolder)) deallocate(cellsHolder)
            allocate(cellsHolder(grid%CellCount))
            cellsHolder(:) = 0

            ! Read cells
            if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
              call u3dintmp(srcUnit, outUnit, grid%LayerCount, grid%RowCount, &
                grid%ColumnCount, grid%CellCount, cellsHolder, aname(1)) 
            else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
              if ( .not. allocated( cellsPerLayer ) ) then
                allocate(cellsPerLayer(grid%LayerCount))
                do n = 1, grid%LayerCount
                  cellsPerLayer(n) = grid%GetLayerCellCount(n)
                end do
              end if 
              call u3dintmpusg(srcUnit, outUnit, grid%CellCount, grid%LayerCount, &
                cellsHolder, aname(1), cellsPerLayer)
            else
              write(outUnit,*) 'Invalid grid type specified when reading CELLS array data.'
              write(outUnit,*) 'Stop.'
              call ustop(' ')          
            end if

            ! Count how many obs cells specified 
            nCells = count(cellsHolder/=0)
            if ( nCells .eq. 0 ) then 
              write(outUnit,*) 'No cells found in array for source.' 
              call ustop('No cells found in array for source. Stop.')
            end if

            ! Depending on the number of cells 
            ! allocate array for cell ids
            if ( allocated( srcCellNumbers ) ) deallocate( srcCellNumbers )
            allocate( srcCellNumbers(nCells) )

            ! Fill srcCellNumbers
            cellCounter = 0
            do nc =1,grid%CellCount
              if(cellsHolder(nc).eq.0) cycle
              cellCounter = cellCounter + 1
              srcCellNumbers(cellCounter) = nc
            end do

            ! Apply default iface if the iface option was given 
            if ( iFaceOption ) then 
              if ( allocated(srcCellIFaces) ) deallocate( srcCellIFaces )
              allocate( srcCellIFaces(nCells) )
              srcCellIFaces(:) = defaultIFaceNumber
            end if

            if ( allocated( cellsHolder ) ) deallocate( cellsHolder ) 
          ! Invalid
          case default
            write(outUnit,'(A,I6)') 'Cells reading format not available. Given ', cellReadFormat
            call ustop('Cells reading format not available. Stop.')
          end select

          ! Read the number of concentrations/species
          read(srcUnit, '(a)') line
          icol = 1
          call urword(line,icol,istart,istop,2,n,r,0,0)
          nSpecies = n 
          if ( nSpecies .lt. 1 ) then 
            write(outUnit,'(A)') 'Number of species/concentration columns should be at least 1.'
            call ustop('Number of species/concentration columns should be at least 1. Stop.')
          end if 
          if ( nSpecies .gt. 50 ) then 
            write(outUnit,'(A)') 'Number of concentration columns is large. An arbitrary limit of 50 is established.'
            call ustop('Number of concentration columns is large. An arbitrary limit of 50 is established. Stop.')
          end if 

          ! Read the option to determine how to read 
          ! the release template
          ! If the option:
          ! - not given or 0  : read 1 template and apply the same for all species
          ! - given and .gt. 0: read nSpecies templates
          readNTemplates = 1
          call urword(line,icol,istart,istop,2,n,r,0,0)
          if ( n.gt.0 ) then 
            readNTemplates = nSpecies
          end if  
          if( allocated( nSubDivisions ) ) deallocate( nSubDivisions ) 
          allocate( nSubDivisions(nSpecies,3) )
          nSubDivisions(:,:) = 1
          do nr=1,readNTemplates
            read(srcUnit,*) (nSubDivisions(nr,ns),ns=1,3)
          end do
          if ( (readNTemplates .eq. 1) .and. (nSpecies.gt.1) ) then 
            do ns=2,nSpecies
              nSubDivisions(ns,:) = nSubDivisions(1,:)
            end do
          end if
          ! Verify validity  
          if ( any(nSubDivisions.lt.1) ) then 
            write(outUnit,'(A)') 'Error: some invalid subdivisions. Verify that all are .gt. 0.'
            call ustop('Error: some invalid subdivisions. Verify that all are .gt. 0. Stop.')
          end if 

          ! It needs to read particles mass, for each species
          if( allocated( particlesMass ) ) deallocate( particlesMass )
          allocate( particlesMass( nSpecies ) )
          particlesMass(:) = 0d0 
          ! Read masses old school
          read(srcUnit,*) (particlesMass(ns),ns=1,nSpecies)

          ! And some validation
          ! Invalid if:
          ! - A value of mass is .le. 0
          if ( any(particlesMass.le.0) ) then 
            write(outUnit,'(A)') 'Error: invalid particles mass. Verify that all are .gt. 0.'
            call ustop('Error: invalid particles mass. Verify that all are .gt. 0. Stop.')
          end if 

          ! Read the solute ids if the simulation demands  
          if ( ( this%ParticlesMassOption .eq. 2 ) ) then 
            write(outUnit,'(A)') 'Will read species ids due to ParticlesMassOption.eq.2 '
            ! Read solute id
            if( allocated( soluteIds ) ) deallocate( soluteIds )
            allocate( soluteIds( nSpecies ) )
            read(srcUnit,*) (soluteIds(ns),ns=1,nSpecies)
            do ns=1,nSpecies
             if ( soluteIds(ns) .lt. 1 ) then 
              write(outUnit,'(A)') 'SpeciesID for source is invalid. It was not given or is less than 1. Stop.'
              call ustop('SpeciesID for source is invalid. It was not given or is less than 1. Stop.')
             end if 
            end do
          end if
   
          ! Read the number of time intervals
          read(srcUnit, '(a)') line
          icol = 1
          call urword(line,icol,istart,istop,2,n,r,0,0)
          nTimeIntervals = n 
          if ( nTimeIntervals .lt. 1 ) then 
            write(outUnit,'(A)') 'Number of time intervals should be at least 1.'
            call ustop('Number of time intervals should be at least 1. Stop.')
          end if 
          if ( nTimeIntervals .gt. 1000000 ) then 
            write(outUnit,'(A)') 'Number of time intervals is large. An arbitrary limit of 1e6 is established.'
            call ustop('Number of time intervals is large. An arbitrary limit of 1e6 is established. Stop.')
          end if 

          ! Allocate allSpecData
          if ( allocated( allSpecData ) ) deallocate( allSpecData ) 
          ! NOTICE THE NCELLS,NTIMES TO FOLLOW FORTRAN MAJOR FORMAT
          if ( .not. concPerCell ) then 
            nColumns = nSpecies + 2 
            allocate( allSpecData(nColumns,nTimeIntervals) )
          else
            nColumns = nSpecies*nCells + 2 
            allocate( allSpecData(nColumns,nTimeIntervals) )
          end if 
        
          ! Load everything in allSpecData. 
          ! "columns" 1 and 2 (rows in fact, but seen as columns in input file), represent the time intervals
          ! and the rest are the values for species/concentrations
          do nt=1,nTimeIntervals           
           read(srcUnit,*) (allSpecData(nd,nt),nd=1,nColumns)
          end do

          ! Validate the time values
          ! Invalid intervals when:
          ! Any tend .le. tstart
          ! tstart(nt+1) .lt. tend(nt)
          isValid = .true.
          do nt=1,nTimeIntervals
            ! If end .le. start
            if( allSpecData(2,nt) .le. allSpecData(1,nt) ) then 
              isValid = .false.
              exit
            end if
            ! If the last loop and didn't fail
            ! in previous check, ready
            if ( nt.eq.nTimeIntervals ) then 
              exit
            end if 
            ! If start+1 .lt. end
            if( allSpecData(1,nt+1) .lt. allSpecData(2,nt) ) then 
              isValid = .false.
              exit
            end if  
          end do
          if( .not. isValid ) then 
          write(outUnit,'(A)') 'Error: invalid time intervals. Verify that do not overlap and that start is .le. end.'
          call ustop('Error: invalid time intervals. Verify that do not overlap and that start is .le. end. Stop.')
          end if 

          ! It needs to build a time vector considering characteristic simulation times.
          ! It was already validated that input times were increasing. 

          ! Define the times for extracting information from MODFLOW.
          ! The initial time is always ReferenceTime.
          initialTime = this%ReferenceTime 
          ! The final MODFLOW time...
          ! StopTime was updated at MPathRW, so it already contains the logic of StoppingTimeOption.
          ! In cases that stoptime is infinite, force the MODFLOW limits 
          if ( this%TrackingOptions%BackwardTracking ) then 
            finalTime = this%ReferenceTime - this%StopTime
            if ( finalTime .lt. 0d0 ) finalTime = 0d0
          else 
            finalTime = this%ReferenceTime + this%StopTime
            if ( finalTime.gt.this%tdisData%TotalTimes(size(this%tdisData%TotalTimes)) ) then 
              finalTime = this%tdisData%TotalTimes(size(this%tdisData%TotalTimes)) 
            end if
          end if 


          ! Given initial and final times, 
          ! compute the initial and final time step indexes.
          kinitial = this%tdisData%FindContainingTimeStep(initialTime)
          kfinal   = this%tdisData%FindContainingTimeStep(finalTime)
          if ( (kfinal .eq. 0) ) then ! Could it be one ?
            write(outUnit,'(a)') 'Warning: kfinal is assumed to be CumulativeTimeStepCount'
            write(outUnit,'(a,e15.7)') 'finalTime is ', finalTime
            kfinal = this%tdisData%CumulativeTimeStepCount
          end if

          kdelta = 1
          correctInterval = 1
          ! Modify values for backward tracking
          if ( this%TrackingOptions%BackwardTracking ) then 
            !sign   = -1d0
            kdelta = -1 
            ! This verification avoids creating an additional unnecessary interval.
            ! Taking the previous example, if initialTime 1.5dt, and finalTime is dt, 
            ! FindContainingTimeStep returns 2 and 1 respectively, hence nTimeIntervals 
            ! is 2 if computed as abs(kfinal-kinitial)+1, when in reality is only 1 interval.
            if ( finalTime.eq.this%tdisData%TotalTimes(kfinal) ) correctInterval = 0
          end if

          ! The number of intervals
          nMFTimeIntervals = abs(kfinal - kinitial) + correctInterval
          nMFTimes = nMFTimeIntervals + 1 
          ! Something wrong with times 
          if ( nMFTimeIntervals .lt. 1 ) then 
            write(message,'(A)') 'Error: the number of time intervals is .lt. 1. Stop.'
            message = trim(message)
            call ustop(message)
          end if

          ! mfTimes: obtained from determined kinitial and kfinal
          ! including initial and final times 
          if ( allocated( mfTimes ) ) deallocate( mfTimes )
          allocate( mfTimes(nMFTimes) )

          ! Fill mfTimes
          !do ktime=1,nMFTimeIntervals
          !  if ( ktime .eq. 1 ) mfTimes(1) = initialTime
          !  if ( ktime .lt. nMFTimeIntervals ) mftimes(ktime+1)  = this%tdisData%TotalTimes(ktime+kinitial-1)
          !  if ( ktime .eq. nMFTimeIntervals ) mftimes(nMFTimes) = finalTime 
          !end do
          ! Fill times and time intervals
          if ( this%TrackingOptions%BackwardTracking ) then 
            do ktime=1,nMFTimeIntervals
              if ( ktime .eq. 1 ) mftimes(1) = initialTime
              if ( ktime .lt. nMFTimeIntervals ) mftimes(ktime+1) = this%tdisData%TotalTimes(kinitial-ktime)
              if ( ktime .eq. nMFTimeIntervals ) mftimes(nMFTimes) = finalTime 
              !timeIntervals(ktime) = abs(times(ktime+1) - times(ktime))
            end do
          else
            do ktime=1,nMFTimeIntervals
              if ( ktime .eq. 1 ) mftimes(1) = initialTime
              if ( ktime .lt. nMFTimeIntervals ) mftimes(ktime+1) = this%tdisData%TotalTimes(ktime+kinitial-1)
              if ( ktime .eq. nMFTimeIntervals ) mftimes(nMFTimes) = finalTime 
              !timeIntervals(ktime) = times(ktime+1) - times(ktime)
            end do
          end if 

          ! Allocate input times, will have at most nTimeIntervals*2 elements
          if ( allocated( inputTimes ) ) deallocate( inputTimes ) 
          allocate( inputTimes(2*nTimeIntervals) )
          inputTimes (:) = 0d0

          ! Flatten the input times
          itCounter = 0
          if ( this%TrackingOptions%BackwardTracking ) then 
            ! Decreasing order
            do nt = nTimeIntervals, 1, -1
              if (nt.eq.nTimeIntervals) then
                inputTimes(1) = allSpecData(2,nt)
                inputTimes(2) = allSpecData(1,nt)
                itCounter = 2
                cycle
              end if
              if ( nt.ge.1 ) then 
                if ( allSpecData(2,nt) .lt. allSpecData(1,nt+1) ) then
                 itCounter = itCounter + 1
                 inputTimes(itCounter) = allSpecData(2,nt)
                end if
                itCounter = itCounter + 1
                inputTimes(itCounter) = allSpecData(1,nt)
              end if 
            end do
          else
            ! Increasing order
            do nt = 1, nTimeIntervals
              if (nt.eq.1) then
                inputTimes(1) = allSpecData(1,nt)
                inputTimes(2) = allSpecData(2,nt)
                itCounter = 2
                cycle
              end if
              if ( nt.le.nTimeIntervals ) then
                if ( allSpecData(1,nt) .gt. allSpecData(2,nt-1) ) then
                 itCounter = itCounter + 1
                 inputTimes(itCounter) = allSpecData(1,nt)
                end if
                itCounter = itCounter + 1
                inputTimes(itCounter) = allSpecData(2,nt)
              end if
            end do
          end if 

          ! Allocate mergedTimes, will have at most 
          ! tCounter + nMFTimes elements
          if ( allocated( mergedTimes ) ) deallocate( mergedTimes )
          allocate( mergedTimes(itCounter+nMFTimes) )
          mergedTimes(:) = 0d0
          tCounter = 1
          mftCounter = 1

          if ( this%TrackingOptions%BackwardTracking ) then
            ! Merge in decreasing order 
            do nt = 1, (itCounter+nMFTimes)
             if( (tCounter.le.itCounter).and.(mftCounter.le.nMFTimes) ) then 
               minT = maxval((/inputTimes(tCounter),mfTimes(mftCounter)/))
             else
               if ( tCounter.gt.itCounter ) minT = mfTimes(mftCounter) 
               if ( mftCounter.gt.nMFTimes ) minT = inputTimes(tCounter)
             end if 
             if ( minT.eq.inputTimes(tCounter)) tCounter = tCounter + 1
             if ( minT.eq.mfTimes(mftCounter) ) mftCounter = mftCounter + 1
             mergedTimes(nt) = minT
             if( (tCounter.gt.itCounter).and.(mftCounter.gt.nMFTimes) ) exit
            end do 
            ! Of the merged vector, how many are within the valid range 
            ! The range (1:nt) is to avoid the potential zeros at the end of the array
            tCounter =  count((mergedTimes(1:nt).ge.finalTime).and.(mergedTimes(1:nt).le.initialTime))
            ! Filter times and fill the time vector
            if( allocated( times ) ) deallocate( times )
            allocate( times(tCounter) )
            times(:) = 0d0
            tcount = 0
            do nt = 1, (itCounter+nMFTimes)
              if( mergedTimes(nt) .gt. initialTime ) cycle
              if( mergedTimes(nt) .lt. finalTime ) exit
              tcount = tcount + 1 
              if ( tcount .gt. tCounter ) exit
              times(tcount) = mergedTimes(nt)
            end do
            ! For the vector of times, determine the input data interval.
            ! Will be used to assign concentrations. 
            if ( allocated( intervalIndex ) ) deallocate( intervalIndex ) 
            allocate( intervalIndex( size(times) ) ) 
            intervalIndex(:) = 0
            itCounter = nTimeIntervals
            do nt=2,size(times)
              doCounter = 0
              do
                if( itCounter .lt. 1 ) exit
                if( times(nt).ge.allSpecData(2,itCounter) ) then
                  exit
                end if
                if(&
                  (times(nt-1).le.allSpecData(2,itCounter)).and.& 
                  (times(nt).ge.allSpecData(1,itCounter)) ) then
                  intervalIndex(nt) = itCounter
                  exit
                end if
                if(&
                  (times(nt-1).le.allSpecData(1,itCounter)).and.& 
                  (times(nt).ge.allSpecData(2,itCounter-1)) ) then
                  itCounter = itCounter - 1
                  exit
                end if
                doCounter = doCounter + 1
                if( doCounter .gt. 1e5 ) exit ! just in case
              end do
              if( itCounter .lt. 1 ) exit
              if( doCounter .gt. 1e5 ) exit
            end do 
          else
            ! Merge in increasing order
            do nt = 1, (itCounter+nMFTimes)
             if( (tCounter.le.itCounter).and.(mftCounter.le.nMFTimes) ) then 
               minT = minval((/inputTimes(tCounter),mfTimes(mftCounter)/))
             else
               if ( tCounter.gt.itCounter ) minT = mfTimes(mftCounter) 
               if ( mftCounter.gt.nMFTimes ) minT = inputTimes(tCounter)
             end if 
             if ( minT.eq.inputTimes(tCounter)) tCounter = tCounter + 1
             if ( minT.eq.mfTimes(mftCounter) ) mftCounter = mftCounter + 1
             mergedTimes(nt) = minT
             if( (tCounter.gt.itCounter).and.(mftCounter.gt.nMFTimes) ) exit
            end do 
            ! Of the merged vector, how many are within the valid range 
            ! The range (1:nt) is to avoid the potential zeros at the end of the array
            tCounter =  count((mergedTimes(1:nt).ge.initialTime).and.(mergedTimes(1:nt).le.finalTime))
            ! Filter times and fill the time vector
            if( allocated( times ) ) deallocate( times )
            allocate( times(tCounter) )
            times(:) = 0d0
            tcount = 0
            do nt = 1, (itCounter+nMFTimes)
              if( mergedTimes(nt) .lt. initialTime ) cycle
              if( mergedTimes(nt) .gt. finalTime ) exit
              tcount = tcount + 1 
              if ( tcount .gt. tCounter ) exit
              times(tcount) = mergedTimes(nt)
            end do
            ! For the vector of times, determine the input data interval.
            ! Will be used to assign concentrations. 
            if ( allocated( intervalIndex ) ) deallocate( intervalIndex ) 
            allocate( intervalIndex( size(times) ) ) 
            intervalIndex(:) = 0
            itCounter = 1
            do nt=2,size(times)
              doCounter = 0
              do
                if( itCounter .gt. nTimeIntervals ) exit
                if( times(nt).le.allSpecData(1,itCounter) ) then
                  exit
                end if
                if(&
                  (times(nt-1).ge.allSpecData(1,itCounter)).and.& 
                  (times(nt).le.allSpecData(2,itCounter)) ) then
                  intervalIndex(nt) = itCounter
                  exit
                end if
                if(&
                  (times(nt-1).ge.allSpecData(2,itCounter)).and.& 
                  (times(nt).le.allSpecData(1,itCounter+1)) ) then
                  itCounter = itCounter + 1
                  exit
                end if
                doCounter = doCounter + 1
                if( doCounter .gt. 1e5 ) exit ! just in case
              end do
              if( itCounter .gt. nTimeIntervals ) exit
              if( doCounter .gt. 1e5 ) exit
            end do 
          end if 
          if ( doCounter.gt.1e5 ) then 
          write(outUnit,'(A)') 'Error: something went wrong while analyzing time intervals for assigning concentrations.'
          call ustop('Error: something went wrong while analyzing time intervals for assigning concentrations. Stop.')
          end if 

          ! Some cleaning
          deallocate( mergedTimes ) 
          deallocate( inputTimes  )

          ! Once the times are known, it can validate the existence of the budget header using the 
          ! range of stress periods within the initial and final times. Validate given aux names.
          isValid = .false.
          isValid = flowModelData%ValidateBudgetHeader(srcPkgNames(nsb),&
                         initialTime, finalTime, this%tdisData, outUnit,&
                      this%TrackingOptions%BackwardTracking, this%isMF6 )
          if ( .not. isValid ) then 
            write(message,'(A,A,A)') 'Given header ', trim(adjustl(srcPkgNames(nsb))),' was not found in budget file. Stop.'
            call ustop(message)
          end if

          ! If the header exists, and it was specified 
          ! to extract the cells from the budget, do it.
          ! Generate the flow-rates timeseries
          call flowModelData%LoadFlowTimeseries( srcPkgNames(nsb), &
            initialTime, finalTime, this%tdisData, srcCellNumbers, &
                 flowDataTimeseries, readCellsFromBudget, outUnit, & 
                            this%TrackingOptions%BackwardTracking, &
                                                        this%isMF6 )

          ! The allocation of concTimeseries needs to be after loadflowtimeseries
          ! in case nCells is determined after reading cells from the budget.
          nCells = size(srcCellNumbers)
          ! Assign concentrations and fill timeIntervals
          nTimeIntervals = size(times)-1
          if ( allocated( concTimeseries ) ) deallocate( concTimeseries ) 
          allocate( concTimeseries( nTimeIntervals, nCells, nSpecies ) )  ! conc during the interval
          concTimeseries(:,:,:) = 0d0

          if ( allocated( timeIntervals ) ) deallocate( timeIntervals ) 
          allocate( timeIntervals( nTimeIntervals ) )  ! interval length
          timeIntervals(:) = 0d0

          do nt=1,size(times)-1
            if( intervalIndex(nt+1).gt.0 ) then 
              ! Assume same concentration for all cells
              if ( .not. concPerCell ) then 
                do nc=1,nCells
                  concTimeseries(nt,nc,:) = allSpecData(3:,intervalIndex(nt+1))
                end do
              else
                do nc=1,nCells
                  concTimeseries(nt,nc,:) = & 
                    allSpecData(2+(nc-1)*nSpecies+1:2+nc*nSpecies,intervalIndex(nt+1))
                end do 
              end if 
            end if
            timeIntervals(nt) = abs(times(nt+1)-times(nt))
          end do 

          ! For the vector of times, determine the corresponding MODFLOW data interval.
          ! Will be used to assign flow-rates. 
          ! Note: mfTimes, by definition, is without blank-jumps in between intervals. 
          ! Meaning that the end of one interval is the beggining of the next. 
          ! mfTimes is built with the consideration that begins at ReferenceTime and 
          ! after intersecting with the information provided by the user, the times vector
          ! also begins at ReferenceTime.
          intervalIndex(:) = 0
          itCounter = 1
          if ( this%TrackingOptions%BackwardTracking ) then 
            do nt=2,size(times)
              doCounter = 0
              do
                if( itCounter .gt. nMFTimeIntervals ) exit
                if(&
                  (times(nt-1).le.mfTimes(itCounter)).and.& 
                  (times(nt).ge.mfTimes(itCounter+1)) ) then
                  intervalIndex(nt) = itCounter
                  if( (times(nt).eq.mfTimes(itCounter+1)) ) itCounter = itCounter + 1
                  exit
                end if
                doCounter = doCounter + 1
                if( doCounter .gt. 1e5 ) exit ! just in case
              end do
              if( itCounter .gt. nMFTimeIntervals ) exit
              if( doCounter .gt. 1e5 ) exit
            end do 
          else
            do nt=2,size(times)
              doCounter = 0
              do
                if( itCounter .gt. nMFTimeIntervals ) exit
                if(&
                  (times(nt-1).ge.mfTimes(itCounter)).and.& 
                  (times(nt).le.mfTimes(itCounter+1)) ) then
                  intervalIndex(nt) = itCounter
                  if( (times(nt).eq.mfTimes(itCounter+1)) ) itCounter = itCounter + 1
                  exit
                end if
                doCounter = doCounter + 1
                if( doCounter .gt. 1e5 ) exit ! just in case
              end do
              if( itCounter .gt. nMFTimeIntervals ) exit
              if( doCounter .gt. 1e5 ) exit
            end do 
          end if
          if ( doCounter.gt.1e5 ) then 
          write(outUnit,'(A)') 'Error: something went wrong while analyzing time intervals for assigning flow-rates.'
          call ustop('Error: something went wrong while analyzing time intervals for assigning flow-rates. Stop.')
          end if

          ! Assign flow-rates 
          if ( allocated( flowTimeseries ) ) deallocate( flowTimeseries ) 
          allocate( flowTimeseries( nTimeIntervals, nCells ) )
          flowTimeseries(:,:) = 0d0
          do nt=1,nTimeIntervals
            if( intervalIndex(nt+1).gt.0 ) then 
              flowTimeseries(nt,:) = flowDataTimeseries(intervalIndex(nt+1),:)
            end if
          end do 
          ! Apply default iface if cells were read from budget
          if ( iFaceOption.and.readCellsFromBudget ) then 
            if ( allocated(srcCellIFaces) ) deallocate( srcCellIFaces )
            allocate( srcCellIFaces(nCells) )
            srcCellIFaces(:) = defaultIFaceNumber
          end if

          ! From here until the creation of particles the process 
          ! is the same than for the AUX format so it could be considered 
          ! to wrap these ops on a common function 

          ! Now compute the cummulative mass function to obtain the total injected mass
          if ( allocated( cummMassTimeseries ) ) deallocate( cummMassTimeseries ) 
          !allocate( cummMassTimeseries, mold=auxTimeseries ) ! auxTimeseries  nt,nc,ns
          allocate( cummMassTimeseries(nTimeIntervals,nCells,nSpecies) ) ! try to be consistent with the aux format 
          cummMassTimeseries(:,:,:) = 0d0

          ! Integrate for each concentration
          ! Considers STEPWISE quantities
          do ns=1, nSpecies
            do nt=1,nTimeIntervals
              if ( nt.gt.1 ) then
                cummMassTimeseries(nt,:,ns) = &
                  cummMassTimeseries(nt-1,:,ns) + flowTimeseries(nt,:)*concTimeseries(nt,:,ns)*timeIntervals(nt)
                cycle
              end if
              if ( nt.eq.1 ) cummMassTimeseries(nt,:,ns) = flowTimeseries(nt,:)*concTimeseries(nt,:,ns)*timeIntervals(nt)
            end do
          end do


          ! NOTICE THE FOLLOWING                   !
          ! TRICKS FOR CONSISTENCY BETWEEN FORMATS !
          nAuxNames = nSpecies
          if( allocated( auxSubDivisions ) ) deallocate( auxSubDivisions ) 
          auxSubDivisions = nSubDivisions
          if ( (.not. iFaceOption).and.allocated(srcCellIFaces) ) then
            deallocate( srcCellIFaces ) 
          end if
          if ( allocated(auxMasses) ) deallocate( auxMasses ) 
          auxMasses = particlesMass
          ! TRICKS FOR CONSISTENCY BETWEEN FORMATS !
          ! NOTICE THE FOLLOWING                   !


          if ( allocated(totMass) ) deallocate(totMass)
          allocate( totMass( nCells ) )
          if ( allocated(totMassLoc) ) deallocate(totMassLoc)
          allocate( totMassLoc( nCells ) )
          if ( allocated(nParticlesDbl) ) deallocate(nParticlesDbl)
          allocate( nParticlesDbl( nCells ) )
          if ( allocated(nParticlesInt) ) deallocate(nParticlesInt)
          allocate( nParticlesInt( nCells ) )
          if ( allocated(nReleases) ) deallocate(nReleases)
          allocate( nReleases( nCells ) )
          if ( allocated( auxEffMasses ) ) deallocate( auxEffMasses ) 
          allocate( auxEffMasses, mold=auxMasses )
          if ( allocated( auxNPCell ) ) deallocate( auxNPCell ) 
          allocate( auxNPCell( nAuxNames, nCells ) )
          ! Allocate srcCellIFaces as a placeholder with zeros
          ! if it was not allocated at the extraction function
          if ( .not. allocated(srcCellIFaces) ) then
            allocate( srcCellIFaces(nCells) )
            srcCellIFaces(:) = 0
          end if


          ! Carrier for candidate particle groups 
          if ( allocated(particleGroups) ) deallocate( particleGroups )
          allocate(particleGroups(nAuxNames))

          ! Intialize some counters
          nValidPGroup = 0
          particleCount = 0

          ! Do it for each aux var
          do naux=1,nAuxNames

            ! Report
            !write(outUnit,'(A,A)') 'Processing auxiliary variable: ', trim(adjustl(auxNames(naux)))
            write(outUnit,'(A,I3)') 'Processing specie: ', naux

            ! Give a name to the particle group
            ! SRCname_idsrcbudget_idauxvar
            tempChar1 = ''
            tempChar2 = ''
            write( unit=tempChar1, fmt=* )nsb 
            write( unit=tempChar2, fmt=* )naux 
            particleGroups(naux)%Name = trim(adjustl(srcName))//'_'//trim(adjustl(tempChar1))//'_'//trim(adjustl(tempChar2))

            ! Increase pgroup counter
            particleGroups(naux)%Group = this%ParticleGroupCount + naux

            ! Assign soluteIds if simulation demands
            if ( this%ParticlesMassOption .eq. 2 ) then 
              particleGroups(naux)%Solute = soluteIds(naux)
            end if

            ! Correction to the particle template based on dimension mask
            ! If dimension is active, leave the given value. 
            ! If not, switch it to one
            do nd = 1, 3
              if ( dimensionMask(nd) .eq. 1 ) cycle 
              auxSubDivisions(naux,nd) = 1
            end do
                
            ! Particles per cell
            auxNPCell(naux,:) = product(auxSubDivisions(naux,:))

            ! Interpret srcCellIFaces and apply
            ! correction to particles per cell.
            if ( iFaceOption ) then
              do nc = 1, nCells
                auxFaceDivisions(:) = 0
                select case(srcCellIFaces(nc))
                case(1,2)
                  auxFaceDivisions(2*srcCellIFaces(nc)-1) = auxSubDivisions(naux,3) ! nver
                  auxFaceDivisions(2*srcCellIFaces(nc))   = auxSubDivisions(naux,2) ! nrow
                case(3,4)
                  auxFaceDivisions(2*srcCellIFaces(nc)-1) = auxSubDivisions(naux,3) ! nver
                  auxFaceDivisions(2*srcCellIFaces(nc))   = auxSubDivisions(naux,1) ! ncol
                case(5,6)
                  auxFaceDivisions(2*srcCellIFaces(nc)-1) = auxSubDivisions(naux,2) ! nrow
                  auxFaceDivisions(2*srcCellIFaces(nc))   = auxSubDivisions(naux,1) ! ncol
                case default
                  cycle
                end select
                ! Update NP cell
                auxNPCell(naux,nc) = auxFaceDivisions(2*srcCellIFaces(nc)-1)*auxFaceDivisions(2*srcCellIFaces(nc))
              end do
            end if

            ! Max cumm mass stats
            totMass    = maxval(cummMassTimeseries(:,:,naux), dim=1) ! Max cumm mass in time
            totMassLoc = maxloc(cummMassTimeseries(:,:,naux), dim=1) ! Loc of max cumm mass in time
            ! Note: the location of max mass corresponds to the end of the interval at index loc

            ! With the given particles mass, estimate the number of particles
            ! necessary for achieving the cummulative mass 
            nParticlesDbl = totMass/auxMasses(naux)
            ! ... and using the number of particles per cell, estimate
            ! the number of releases per cell.
            nReleases     = int(nParticlesDbl/auxNPCell(naux,:)+0.5)

            ! Up to this point, the number of particles 
            ! needed can be computed as nReleases*NPCELL

            ! Cycle to the next aux var if the total 
            ! number of particles is zero.
            totalParticleCount = sum(nReleases*auxNPCell(naux,:))
            if ( totalParticleCount .lt. 1 ) then 
              write(outUnit,*) ' Warning: pgroup ',&
                trim(adjustl(particleGroups(naux)%Name)),' has zero particles, it will skip this aux var.'
              ! Process the next aux variable
              cycle
            end if 

            ! Allocate particles 
            particleGroups(naux)%TotalParticleCount = totalParticleCount
            if(allocated(particleGroups(naux)%Particles)) deallocate(particleGroups(naux)%Particles)
            allocate(particleGroups(naux)%Particles(totalParticleCount))

            ! Variables for defining the particles of this group
            idmax = 0
            offset = 0
            seqNumber = this%currentSeqNumber 
            cellCounter = 0
            currentParticleCount = 0 

            ! Given the number of releases, interpolate release times
            ! for each cell, using as source the cummulative mass function
            totCummEffectiveMass = 0d0
            do nc=1, nCells

              if ( nReleases(nc) .lt. 1 ) cycle ! At least one release

              ! DEPRECATE
              cellCounter = cellCounter + 1
              if ( allocated( releaseTimes ) ) deallocate( releaseTimes ) 
              allocate( releaseTimes( nReleases(nc) ) )
              ! END DEPRECATE

              ! The effective mass of particles for the aux var of this cell
              effectiveMass = totMass(nc)/(nReleases(nc)*auxNPCell(naux,nc))

              ! Create particle template for this cell
              ! 0: is for drape. (TEMP)
              ! Drape = 0: particle placed in the cell. If dry, status to unreleased
              ! Drape = 1: particle placed in the uppermost active cell
              ! -999: is a placeholder for release time, it will assigned later

              ! Interpret srcCellIFaces
              ! This check requires allocated srcCellIFaces, it could be improved. 
              if ( iFaceOption .and. (srcCellIFaces(nc).gt.0)  ) then
                ! Create particles at cell face
                auxFaceDivisions(:) = 0
                select case(srcCellIFaces(nc))
                case(1,2)
                  auxFaceDivisions(2*srcCellIFaces(nc)-1) = auxSubDivisions(naux,3) ! nver
                  auxFaceDivisions(2*srcCellIFaces(nc))   = auxSubDivisions(naux,2) ! nrow
                case(3,4)
                  auxFaceDivisions(2*srcCellIFaces(nc)-1) = auxSubDivisions(naux,3) ! nver
                  auxFaceDivisions(2*srcCellIFaces(nc))   = auxSubDivisions(naux,1) ! ncol
                case(5,6)
                  auxFaceDivisions(2*srcCellIFaces(nc)-1) = auxSubDivisions(naux,2) ! nrow
                  auxFaceDivisions(2*srcCellIFaces(nc))   = auxSubDivisions(naux,1) ! ncol
                end select
                call CreateMassParticlesOnFaces(&
                  particleGroups(naux),    &
                  srcCellNumbers(nc),      &
                  currentParticleCount,    &
                  auxFaceDivisions,        &
                  0, effectiveMass, -999d0 )
              else
                ! Normal internal array of particles
                call CreateMassParticlesAsInternalArray(& 
                  particleGroups(naux),     &
                  srcCellNumbers(nc),       &
                  currentParticleCount,     &
                  auxSubDivisions(naux,1),  &
                  auxSubDivisions(naux,2),  &
                  auxSubDivisions(naux,3),  & 
                  0, effectiveMass, -999d0  )
              end if
              ! Assign values for particle template
              !   - layer
              !   - group
              do m = 1, currentParticleCount 
                particleGroups(naux)%Particles(m)%Group = particleGroups(naux)%Group
                particleGroups(naux)%Particles(m)%InitialLayer =                       &
                  grid%GetLayer(particleGroups(naux)%Particles(m)%InitialCellNumber)
                particleGroups(naux)%Particles(m)%Layer =                              &
                  grid%GetLayer(particleGroups(naux)%Particles(m)%CellNumber)
              end do

              ! Note: 
              ! Interpolator requires strictly increasing x and 
              ! the cummulative mass function may have flat sections.

              ! Array for the cumm series considering a zero at the beginning
              if ( allocated( cummMassSeries ) ) deallocate( cummMassSeries ) 
              allocate( cummMassSeries(totMassLoc(nc)+1) )
              cummMassSeries(:)  = 0d0
              cummMassSeries(2:) = cummMassTimeseries(1:totMassLoc(nc),nc,naux)

              ! nti, nte: time indexes to be passed to interpolator as source
              nti = 0
              nte = 0
              firstnonzero = 0
              lastCummMass = 0d0
              cummEffectiveMass = 0d0
              npNextRelease = 0
              ! Loop over cummulative mass series
              do nt = 1, totMassLoc(nc)+1
                ! Advance loop until finding the first non-zero
                if (&
                  ( cummMassSeries( nt ) .eq. 0d0 ) .and. &
                  ( firstnonzero .eq. 0 ) ) cycle

                ! If found a nonzero mass, save the first nonzero
                if ( firstnonzero .eq. 0 ) then 
                  nti = nt-1
                  if( nt .eq. 1 ) nti = 1 ! just in case, but it shouldn't
                  firstnonzero = nt
                  lastCummMass = cummMassSeries(nt)
                  if( nt .ne. (totMassLoc(nc)+1) ) cycle 
                end if
                ! If the current mass is higher than the last found, 
                ! and the initial index is defined, then everything is ok, continue
                ! to the next index if it is not the last.
                if ( (cummMassSeries(nt) .gt. lastCummMass) .and. (nti.gt.0) ) then 
                  lastCummMass = cummMassSeries(nt)
                  if( nt.ne.(totMassLoc(nc)+1) ) cycle
                else if ( (cummMassSeries(nt) .eq. lastCummMass) .and. (nti.eq.0) ) then
                  ! Still on a flat area
                  cycle
                else if ( (cummMassSeries(nt) .gt. lastCummMass) .and. (nti.eq.0) ) then
                  ! Define the new starting point and update lastCummMass
                  nti = nt - 1
                  lastCummMass = cummMassSeries(nt)
                  if( nt.ne.(totMassLoc(nc)+1) ) cycle
                end if

                ! If it is equal to the last, and the first 
                ! index is defined, then we are on a flat zone
                if ( (cummMassSeries(nt) .eq. lastCummMass) .and. (nti.gt.0) ) then
                  ! This is the last index for interpolation 
                  nte = nt-1 
                  ! If it is the last loop, set at the last index in array
                  if ( nt.eq.(totMassLoc(nc)+1) ) nte = nt

                  ! Can it occur ?
                  if ( nte .eq. nti ) then 
                    write(outUnit,'(A)') 'Range of time indexes for 1d interpolation is inconsistent.'
                    call ustop('Range of time indexes for 1d interpolation is inconsistent. Stop.')
                  end if

                  ! Initialize interpolator
                  call interp1d%initialize(cummMassSeries(nti:nte),times(nti:nte), int1dstat)
                  if ( int1dstat .gt. 0 ) then 
                    write(outUnit,'(A)') 'There was a problem while initializing 1d interpolator.'
                    call ustop('There was a problem while initializing 1d interpolator. Stop.')
                  end if

                  ! Initialize release variables
                  lastRelease = .false.
                  npThisRelease = auxNPCell(naux,nc)
                  if ( npNextRelease .gt. 0 ) then 
                    npThisRelease = npNextRelease
                    npNextRelease = 0
                  end if 

                  ! And loop until a maximum of nReleases
                  ! It will break based on cummulative mass
                  do nr=1, nReleases(nc)
                    ! Increase the cummulative mass counter
                    cummEffectiveMass = cummEffectiveMass + npThisRelease*effectiveMass

                    if ( cummEffectiveMass .ge. lastCummMass ) then
                      ! If above the current limit, add some particles
                      ! as long as the limit remains below lastCummMass 
                      ! and correct cummEffectiveMass
                      cummEffectiveMass = cummEffectiveMass - npThisRelease*effectiveMass
                      npAuxDbl = (lastCummMass - cummEffectiveMass)/effectiveMass
                      npThisRelease = int(npAuxDbl+0.5)

                      ! If no particles, done
                      if ( npThisRelease .lt. 1 ) exit

                      ! This variable keeps track of the particles 
                      ! needed to begin releases next time and keep 
                      ! consistency with total mass for the computed
                      ! effectiveMass
                      npNextRelease = auxNPCell(naux,nc) - npThisRelease
                      cummEffectiveMass = cummEffectiveMass + npThisRelease*effectiveMass
                      lastRelease = .true.
                    end if

                    ! Interpolate the release time        
                    call interp1d%evaluate( cummEffectiveMass, releaseTimes(nr) ) ! THIS ARRAY RELEASE TIMES BYE BYE

                    ! Fix the release time considering backward/forward tracking
                    if ( this%TrackingOptions%BackwardTracking ) then
                      ! The interpolated value is the time in the context of the MODFLOW model, 
                      ! if backward tracking, this time is decreasing with respect to reference time.
                      releaseTrackingTime = this%ReferenceTime - releaseTimes(nr)
                    else
                      ! In forward tracking, the interpolated time is increasing with respect to 
                      ! reference time 
                      releaseTrackingTime = releaseTimes(nr) - this%ReferenceTime
                    end if 

                    ! Assign to particles
                    do m = 1, npThisRelease
                     idmax = idmax + 1
                     seqNumber = seqNumber + 1
                     particleGroups(naux)%Particles(offset+m)%Id = idmax
                     particleGroups(naux)%Particles(offset+m)%SequenceNumber = seqNumber
                     particleGroups(naux)%Particles(offset+m)%Group = particleGroups(naux)%Group
                     particleGroups(naux)%Particles(offset+m)%Drape = particleGroups(naux)%Particles(m)%Drape
                     particleGroups(naux)%Particles(offset+m)%Status = particleGroups(naux)%Particles(m)%Status
                     particleGroups(naux)%Particles(offset+m)%InitialCellNumber = & 
                             particleGroups(naux)%Particles(m)%InitialCellNumber
                     particleGroups(naux)%Particles(offset+m)%InitialLayer = particleGroups(naux)%Particles(m)%InitialLayer
                     particleGroups(naux)%Particles(offset+m)%InitialFace = particleGroups(naux)%Particles(m)%InitialFace
                     particleGroups(naux)%Particles(offset+m)%InitialLocalX = particleGroups(naux)%Particles(m)%InitialLocalX
                     particleGroups(naux)%Particles(offset+m)%InitialLocalY = particleGroups(naux)%Particles(m)%InitialLocalY
                     particleGroups(naux)%Particles(offset+m)%InitialLocalZ = particleGroups(naux)%Particles(m)%InitialLocalZ
                     particleGroups(naux)%Particles(offset+m)%InitialTrackingTime = releaseTrackingTime
                     particleGroups(naux)%Particles(offset+m)%TrackingTime = & 
                             particleGroups(naux)%Particles(offset+m)%InitialTrackingTime
                     particleGroups(naux)%Particles(offset+m)%CellNumber = particleGroups(naux)%Particles(m)%CellNumber
                     particleGroups(naux)%Particles(offset+m)%Layer = particleGroups(naux)%Particles(m)%Layer
                     particleGroups(naux)%Particles(offset+m)%Face = particleGroups(naux)%Particles(m)%Face
                     particleGroups(naux)%Particles(offset+m)%LocalX = particleGroups(naux)%Particles(m)%LocalX
                     particleGroups(naux)%Particles(offset+m)%LocalY = particleGroups(naux)%Particles(m)%LocalY
                     particleGroups(naux)%Particles(offset+m)%LocalZ = particleGroups(naux)%Particles(m)%LocalZ
                     ! Mass ( potentially different effectiveMass for different cells )
                     particleGroups(naux)%Particles(offset+m)%Mass   = effectiveMass
                    end do

                    ! Increase offset
                    offset = offset + npThisRelease
                    
                    ! Restart npThisRelease
                    npThisRelease = auxNPCell(naux,nc)

                    if ( lastRelease ) then
                      ! And exit loop over releases
                      exit
                    end if 

                  end do ! nr=1, nReleases(nc)

                  ! After interpolating release times, destroy interpolator
                  call interp1d%destroy()

                  ! It should reset nti, nte
                  nti = 0
                  nte = 0

                end if

              end do ! nt = 1, totMassLoc(nc)+1

              ! Save the current particle count, it will add more
              ! corresponding to the next cell
              currentParticleCount = offset

              ! Accumulate over all cells for report 
              totCummEffectiveMass = totCummEffectiveMass + cummEffectiveMass

            end do ! nc=1, nCells

            ! Save the current sequence number once finished the 
            ! processing of this source 
            this%currentSeqNumber = seqNumber

            ! Increase counters
            if ( particleGroups(naux)%TotalParticleCount .gt. 0 ) then 
              nValidPGroup = nValidPGroup + 1
              particleCount =  particleCount + particleGroups(naux)%TotalParticleCount 
              write(outUnit,'(A,es18.9e3)') 'Total released mass related to aux var  = ', totCummEffectiveMass
              write(outUnit,'(A,I10)') 'Total number of particles related to aux var = ', particleGroups(naux)%TotalParticleCount 
            end if


          end do ! naux = 1, nAuxNames
       
          ! It needs to add the newly created groups to simulationData%ParticleGroups

          ! Extend simulationdata to include these particle groups
          if ( nValidPGroup .gt. 0 ) then 
            newParticleGroupCount = this%ParticleGroupCount + nValidPGroup
            allocate(newParticleGroups(newParticleGroupCount))
            ! If some particle groups existed previously
            if( this%ParticleGroupCount .gt. 0 ) then 
              do n = 1, this%ParticleGroupCount
                newParticleGroups(n) = this%ParticleGroups(n)
              end do
            end if 
            pgCount = 0
            do n = 1, nAuxNames
              if ( particleGroups(n)%TotalParticleCount .eq. 0 ) cycle
              pgCount = pgCount + 1 
              newParticleGroups(pgCount+this%ParticleGroupCount) = particleGroups(n)
            end do 
            if( this%ParticleGroupCount .gt. 0 ) then 
              call move_alloc( newParticleGroups, this%ParticleGroups )
              this%ParticleGroupCount = newParticleGroupCount
              this%TotalParticleCount = this%TotalParticleCount + particleCount
            else
              this%ParticleGroupCount = newParticleGroupCount
              this%TotalParticleCount = particleCount
              !allocate(this%ParticleGroups(this%ParticleGroupCount))
              call move_alloc( newParticleGroups, this%ParticleGroups )
            end if
          end if


        end do ! nsb=1,nSrcBudgets


      case default
        ! Not implemented
        write(message,'(A,A,A)') 'The SRC specification kind ',trim(adjustl(srcSpecKind)),' has not been implemented. Stop.'
        message = trim(message)
        write(outUnit,'(A)') message
        call ustop(message)
      end select
      write(outUnit,'(A,A,A,I10)')'Total number of particles in ',trim(adjustl(srcName)),'       = ', particleCount
      write(outUnit,'(1x,a)') '------------------------'
      write(outUnit, *)

    end do ! nSources/specs
       
    ! flush
    flush(outUnit)

    ! Close data file
    close( srcUnit )

  end subroutine pr_ReadSRCData


end module ModpathSimulationDataModule
