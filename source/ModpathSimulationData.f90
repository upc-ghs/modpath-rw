module ModpathSimulationDataModule
  use ParticleTrackingOptionsModule,only : ParticleTrackingOptionsType
  use ParticleGroupModule,only : ParticleGroupType
  use ModflowRectangularGridModule,only : ModflowRectangularGridType
  use StartingLocationReaderModule,only : ReadAndPrepareLocations
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
    character(len=200) :: DispersionFile          ! RWPT
    logical            :: isRWSimulation =.false. ! RWPT
    logical            :: anyObservation =.false. ! RWPT
    integer            :: ParticlesMassOption     ! RWPT
    integer            :: SolutesOption           ! RWPT
    integer,dimension(:),allocatable :: BudgetCells
    integer,dimension(:),allocatable :: Zones
    doubleprecision,dimension(:),allocatable :: Retardation
    doubleprecision,dimension(:),allocatable :: TimePoints
    type(ParticleGroupType),dimension(:),allocatable :: ParticleGroups
    type(ParticleTrackingOptionsType),allocatable :: TrackingOptions
  contains
    procedure :: ReadFileHeaders=>pr_ReadFileHeaders
    procedure :: ReadData=>pr_ReadData
    procedure :: ReadGPKDEData=>pr_ReadGPKDEData ! RWPT
    procedure :: ReadOBSData=>pr_ReadOBSData     ! RWPT
  end type


contains


  subroutine pr_ReadFileHeaders(this, inUnit)
  use UTL8MODULE,only : u8rdcom
  implicit none
  class(ModpathSimulationDataType) :: this
  integer,intent(in) :: inUnit
  integer :: outUnit, errorCode
  character(len=200) line
  
  outUnit = 0
  call u8rdcom(inUnit, outUnit, line, errorCode)
  
  ! Assign the name file
  this%NameFile = line
  
  ! Read MODPATH listing file filename
  read(inUnit, '(a)') this%ListingFile
  
  end subroutine pr_ReadFileHeaders
 

  ! Read simulation data 
  subroutine pr_ReadData(this, inUnit, outUnit, ibound, timeDiscretization, grid)
  use UTL8MODULE,only : urword, ustop, u1dint, u1drel, u1ddbl, u8rdcom, u3ddblmpusg, u3dintmp, u3dintmpusg, u3ddblmp, ugetnode
  use TimeDiscretizationDataModule,only : TimeDiscretizationDataType
  use ObservationModule, only: ObservationType
  implicit none
  class(ModpathSimulationDataType), target :: this
  class(ModflowRectangularGridType),intent(in) :: grid
  integer,intent(in) :: inUnit, outUnit
  integer,dimension(:),allocatable :: cellsPerLayer
  integer,dimension(grid%CellCount),intent(in) :: ibound
  type(TimeDiscretizationDataType),intent(in) :: timeDiscretization
  integer :: icol, istart, istop, n, nc, kper, kstp, seqNumber, particleCount, nn, slocUnit, errorCode
  integer :: releaseOption, releaseTimeCount
  doubleprecision :: initialReleaseTime, releaseInterval
  doubleprecision,dimension(:),allocatable :: releaseTimes
  doubleprecision :: frac, r, tinc
  character*24 aname(3)
  character(len=200) line
  DATA aname(1) /'              ZONE ARRAY'/
  DATA aname(2) /'                 RFACTOR'/
  DATA aname(3) /'                OBSCELLS'/

  ! RWPT
  character(len=200) :: dispersionFile
  integer :: iodispersion   = 0
  integer :: ioInUnit       = 0

  ! OBS
  integer :: nObservations  = 0
  character(len=100) :: tempChar
  integer :: layerCount, rowCount, columnCount
  type( ObservationType ), pointer :: obs => null()
  integer :: readStyle, no, cellNumber, layer, row, column, cellCount, ocount
  integer,dimension(:),allocatable :: obsCells

  ! GPKDE
  character(len=200) :: gpkdeFile
  integer :: gpkdeUnit = 89
  integer :: iogpkde   = 0


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
          IF(this%AdvectiveObservationsOption.EQ.2) then
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
          IF(this%AdvectiveObservationsOption.EQ.2) then
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
          IF(this%AdvectiveObservationsOption.EQ.2) then
            read(inUnit, '(a)') this%AdvectiveObservationsFile
            icol = 1
            call urword(this%AdvectiveObservationsFile, icol, istart, istop, 0, n, r,0,0)
            this%AdvectiveObservationsFile = this%AdvectiveObservationsFile(istart:istop)
          end if
          this%isRWSimulation = .true.
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
          IF(this%AdvectiveObservationsOption.EQ.2) then
            read(inUnit, '(a)') this%AdvectiveObservationsFile
            icol = 1
            call urword(this%AdvectiveObservationsFile, icol, istart, istop, 0, n, r,0,0)
            this%AdvectiveObservationsFile = this%AdvectiveObservationsFile(istart:istop)
          end if
          this%isRWSimulation = .true.
      case(7)
          write(outUnit,'(A,I2,A)') 'RWPT Endpoint Analysis (Simulation type =', this%SimulationType, ')'
          read(inUnit, '(a)') this%EndpointFile
          icol = 1
          call urword(this%EndpointFile,icol,istart,istop,0,n,r,0,0)
          this%EndpointFile = this%EndpointFile(istart:istop)
          this%isRWSimulation = .true.
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
  if ( ( ( this%SimulationType .eq. 5 ) .or.  &
         ( this%SimulationType .eq. 6 ) .or.  & 
         ( this%SimulationType .eq. 7 )    )  &
      .and. ( this%TrackingDirection .eq. 2 ) ) then 
      call ustop('Random Walk Particle Tracking only accepts Forward tracking. Stop.')
  end if

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
        write(outUnit,'(/A)') 'The retardation factor array will be read.'
        if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, this%Retardation, ANAME(2))            
        else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,            &
              this%Retardation, aname(2), cellsPerLayer)
        else
            write(outUnit,*) 'Invalid grid type specified when reading retardation array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')            
        end if
    else
        write(outUnit,'(/A)') 'The retardation factor for all cells = 1'
        do n = 1, grid%CellCount
            this%Retardation(n) = 1.0d0
        end do
    end if
      
    ! Particle data
    read(inUnit, *) this%ParticleGroupCount
    write(outUnit,'(/A,I5)') 'Number of particle groups = ', this%ParticleGroupCount
  
    seqNumber = 0
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
          ! Read group mass, is a proxy for concentrations
          ! when mass is uniform for a pgroup
          read(inUnit, *) this%ParticleGroups(n)%Mass
          this%ParticleGroups(n)%Particles(:)%Mass = this%ParticleGroups(n)%Mass
          ! Read the solute id for this group 
          if ( this%ParticlesMassOption .eq. 2 ) then 
            read(inUnit, *) this%ParticleGroups(n)%Solute
          end if
        end if

      end do

      this%TotalParticleCount = particleCount
      write(outUnit, '(a,i10)') 'Total number of particles = ', this%TotalParticleCount
      write(outUnit, *)
    end if

    ! TrackingOptions data
    allocate(this%TrackingOptions)
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

    ! RWPT
    ! Assign specific RWPT options 
    if ( (this%SimulationType .eq. 5) .or. & 
         (this%SimulationType .eq. 6) .or. & 
         (this%SimulationType .eq. 7) ) then

        ! Identify the simulation, used in trackingEngine and probably
        ! in more places. Is not the same as isRWSimulation for 
        ! reading mass from pgroups reading 
        this%TrackingOptions%RandomWalkParticleTracking = .true.

        !! Open the dispersion data file
        !! Requires error handling
        !read(inUnit, '(a)') this%DispersionFile
        !icol = 1
        !call urword(this%DispersionFile,icol,istart,istop,0,n,r,0,0)
        !this%DispersionFile = this%DispersionFile(istart:istop)

        ! Initialization and reading of dispersion file data is handled in TransportModelData.f90        

    end if


  end subroutine pr_ReadData


  ! Read specific GPKDE data
  subroutine pr_ReadGPKDEData( this, gpkdeFile, gpkdeUnit, outUnit )
    use UTL8MODULE,only : urword
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    class(ModpathSimulationDataType), target :: this
    character(len=200), intent(in)           :: gpkdeFile
    integer, intent(in)                      :: gpkdeUnit
    integer, intent(in)                      :: outUnit
    ! local
    integer :: isThisFileOpen = -1
    integer :: icol,istart,istop,n
    doubleprecision    :: r
    character(len=200) :: line
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW GPKDE file data'
    write(outUnit, '(1x,a)') '--------------------------'

    ! Verify if GPKDE unit is open 
    inquire( file=gpkdeFile, number=isThisFileOpen )
    if ( isThisFileOpen .lt. 0 ) then 
      ! No gpkde 
      write(outUnit,'(A)') 'GPKDE reconstruction is disabled'
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
      read(gpkdeUnit, '(a)') this%TrackingOptions%gpkdeOutputFile
      icol = 1
      call urword(this%TrackingOptions%gpkdeOutputFile,icol,istart,istop,0,n,r,0,0)
      this%TrackingOptions%gpkdeOutputFile = this%TrackingOptions%gpkdeOutputFile(istart:istop)
    
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
    
      ! Read binSize
      read(gpkdeUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      this%TrackingOptions%gpkdeBinSize(1) = r
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      this%TrackingOptions%gpkdeBinSize(2) = r
      call urword(line, icol, istart, istop, 3, n, r, 0, 0)
      this%TrackingOptions%gpkdeBinSize(3) = r
    
      ! Read nOptimizationLoops
      read(gpkdeUnit, '(a)') line
      icol = 1
      call urword(line, icol, istart, istop, 2, n, r, 0, 0)
      this%TrackingOptions%gpkdeNOptLoops = n
    
      ! Read reconstruction method
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
        ! Read kernel params
        ! - min   h/lambda
        ! - max   h/lambda
        read(gpkdeUnit, '(a)') line
        icol = 1
        call urword(line, icol, istart, istop, 3, n, r, 0, 0)
        this%TrackingOptions%gpkdeKDBParams(1) = r
        this%TrackingOptions%gpkdeKDBParams(2) = 0d0 ! NOT USED
        call urword(line, icol, istart, istop, 3, n, r, 0, 0)
        this%TrackingOptions%gpkdeKDBParams(3) = r
      end if 
    
    else

      ! If simulation is not timeseries
      write(outUnit,'(A)') 'GPKDE reconstruction requires a timeseries. Will remain disabled.'

    end if

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
    class(ModpathSimulationDataType), target :: this
    character(len=200), intent(in)           :: obsFile
    integer, intent(in)                      :: obsUnit
    integer, intent(in)                      :: outUnit
    class(ModflowRectangularGridType),intent(in) :: grid
    ! local
    integer :: nObservations  = 0
    character(len=100) :: tempChar
    integer :: layerCount, rowCount, columnCount, cellCount
    type( ObservationType ), pointer :: obs => null()
    integer :: readStyle, no, cellNumber, layer, row, column, ocount
    integer,dimension(:),allocatable :: obsCells
    integer,dimension(:),allocatable :: cellsPerLayer
    integer :: isThisFileOpen = -1
    integer :: icol,istart,istop,n,nc
    integer :: ioInUnit = 0
    doubleprecision    :: r
    character(len=200) :: line
    character(len=24)  :: aname(1)
    DATA aname(1) /'                OBSCELLS'/
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW OBS file data'
    write(outUnit, '(1x,a)') '--------------------------'

    ! Verify if OBS unit is open 
    inquire( file=obsFile, number=isThisFileOpen )
    if ( isThisFileOpen .lt. 0 ) then 
      ! No obs
      write(outUnit,'(A)') 'No observations were specified'
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
      this%anyObservation = .true.
      write(outUnit,'(1X,A,I6)') 'Given number of observations: ', nObservations

      ! Allocate observation arrays
      call this%TrackingOptions%InitializeObservations( nObservations )

      ! It might be needed downstream
      layerCount  = grid%LayerCount
      rowCount    = grid%RowCount
      columnCount = grid%ColumnCount
      cellCount   = grid%CellCount
      
      ! Allocate id arrays in tracking options
      if(allocated(this%TrackingOptions%isObservation)) & 
          deallocate(this%TrackingOptions%isObservation)
      allocate(this%TrackingOptions%isObservation(cellCount))
      if(allocated(this%TrackingOptions%idObservation)) & 
          deallocate(this%TrackingOptions%idObservation)
      allocate(this%TrackingOptions%idObservation(cellCount))
      this%TrackingOptions%isObservation(:) = .false.
      this%TrackingOptions%idObservation(:) = -999

      ! Read observation cells and assign 
      ! proper variables
      do nc = 1, nObservations

        ! A pointer
        obs => this%TrackingOptions%Observations(nc) 

        ! Read observation id
        read(obsUnit, '(a)', iostat=ioInUnit) line
        icol = 1
        call urword(line, icol, istart, istop, 2, n, r, 0, 0)
        obs%id = n 
       
        ! Read observation filename and assign an output unit 
        read(obsUnit, '(a)') obs%outputFileName
        icol = 1
        call urword(obs%outputFileName,icol,istart,istop,0,n,r,0,0)
        obs%outputFileName = obs%outputFileName(istart:istop)
        obs%outputUnit     = 5500 + nc
        obs%auxOutputUnit  = 7700 + nc
        tempChar           = 'temp'
        write( unit=obs%auxOutputFileName, fmt='(a)')&
            trim(adjustl(tempChar))//'_'//trim(adjustl(obs%outputFileName))

        ! Read observation style (sink obs, normal count of particles obs)
        read(obsUnit, '(a)', iostat=ioInUnit) line
        icol = 1
        call urword(line, icol, istart, istop, 2, n, r, 0, 0)
        obs%style = n 

        ! Is the style that requires flow-rates ?
        if ( obs%style .eq. 2 ) then 
          this%TrackingOptions%anySinkObservation = .true.
        end if 

        ! Read observation cell option
        ! Determine how to read cells
        read(obsUnit, '(a)', iostat=ioInUnit) line
        icol = 1
        call urword(line, icol, istart, istop, 2, n, r, 0, 0)
        obs%cellOption = n 

        ! Load observation cells
        select case( obs%cellOption )
          ! In case 1, a list of cell ids is specified, that 
          ! compose the observation.  
          case (1)
            ! Read number of observation cells 
            read(obsUnit, '(a)', iostat=ioInUnit) line
            icol = 1
            call urword(line, icol, istart, istop, 2, n, r, 0, 0)
            obs%nCells = n 
            
            ! Depending on the number of cells 
            ! allocate array for cell ids
            if ( allocated( obs%cells ) ) deallocate( obs%cells )
            allocate( obs%cells(obs%nCells) )
            if ( allocated( obs%nRecordsCell ) ) deallocate( obs%nRecordsCell )
            allocate( obs%nRecordsCell(obs%nCells) )
            obs%nRecordsCell(:) = 0

            ! Are these ids as (lay,row,col) or (cellid) ?
            read(obsUnit, '(a)', iostat=ioInUnit) line
            icol = 1
            call urword(line, icol, istart, istop, 2, n, r, 0, 0)
            readStyle = n

            ! Load the observation cells
            if( readStyle .eq. 1) then
              ! Read as layer, row, column
              do no = 1, obs%nCells
                read(obsUnit, *) layer, row, column
                call ugetnode(layerCount, rowCount, columnCount, layer, row, column,cellNumber)
                obs%cells(no) = cellNumber
              end do 
            else if ( readStyle .eq. 2 ) then 
              do no = 1, obs%nCells
                read(obsUnit,*)  cellNumber
                obs%cells(no) = cellNumber
              end do 
            else
              call ustop('Invalid observation kind. Stop.')
            end if

          case (2)
            ! In case 2, observation cells are given by specifying a 3D array
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
            if ( allocated( obs%nRecordsCell ) ) deallocate( obs%nRecordsCell )
            allocate( obs%nRecordsCell(obs%nCells) )
            obs%nRecordsCell(:) = 0

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
          this%TrackingOptions%idObservation(obs%cells(no)) = nc
        end do


        ! Read observation cell time option
        ! Determine how to reconstruct timeseries
        read(obsUnit, '(a)', iostat=ioInUnit) line
        icol = 1
        call urword(line, icol, istart, istop, 2, n, r, 0, 0)
        obs%timeOption = n 

        ! Timeoption determine from where 
        ! to obtain the timeseries considered for 
        ! reconstruction
        select case(obs%timeOption)
          case(1)
            ! Get it from the timeseries run 
            ! Is there any gpkde config for this ?
            continue
          case (2)
            ! Create it by reading input 
            ! params like for example the 
            ! number of datapoints

            ! Needs reading

            continue
          case default
            ! Get it from the timeseries run 
            continue
        end select

        ! Depending on parameters, initialize observation file as 
        ! binary or plain-text 

      end do 

      ! Close the OBS unit
      close( obsUnit ) 

    end if  


  end subroutine pr_ReadOBSData



    
end module ModpathSimulationDataModule
