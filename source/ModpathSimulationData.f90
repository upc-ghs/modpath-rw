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
    procedure :: ReadGPKDEData=>pr_ReadGPKDEData   ! RWPT
    procedure :: ReadOBSData=>pr_ReadOBSData       ! RWPT
    procedure :: ReadRWOPTSData=>pr_ReadRWOPTSData ! RWPT
    procedure :: ReadICData=>pr_ReadICData         ! RWPT
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
 

  ! Read simulation data 
  subroutine pr_ReadData(this, inUnit, outUnit, ibound, timeDiscretization, grid)
    use UTL8MODULE,only : urword, ustop, u1dint, u1drel, u1ddbl, u8rdcom, &
                          u3ddblmpusg, u3dintmp, u3dintmpusg, u3ddblmp, ugetnode
    use TimeDiscretizationDataModule,only : TimeDiscretizationDataType
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
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
    if ((( this%SimulationType .eq. 5 ) .or.  &
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



  ! Read specific RWOPTS data
  subroutine pr_ReadRWOPTSData( this, rwoptsFile, rwoptsUnit, outUnit )
    use UTL8MODULE,only : urword,ustop
    !--------------------------------------------------------------
    ! Specifications
    !--------------------------------------------------------------
    implicit none
    class(ModpathSimulationDataType), target :: this
    character(len=200), intent(in)           :: rwoptsFile
    integer, intent(in)                      :: rwoptsUnit
    integer, intent(in)                      :: outUnit
    ! local
    type(ParticleTrackingOptionsType), pointer :: trackingOptions
    integer :: isThisFileOpen = -1
    integer :: icol,istart,istop,n,nd,currentDim
    doubleprecision    :: r
    character(len=200) :: line
    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW RWOPTS file data'
    write(outUnit, '(1x,a)') '---------------------------'

    ! Verify if unit is open 
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
    if ( trackingOptions%nDim .eq. 0 ) then
      ! No dimensions
      write(outUnit,'(A)') 'No dimensions were given for RW displacements at RWOPTS, nDim .eq. 0.'
      call ustop('No dimensions were given for RW displacements at RWOPTS, nDim .eq. 0. Stop.')
    end if 

    ! Detect idDim and report dimensions
    ! where displacements will be applied
    select case(trackingOptions%nDim)
      ! 1D
      case(1)
        trackingOptions%twoDimensions = .true. ! TEMP
        write(outUnit,'(A)') 'RW displacements for 1 dimension.'
        ! Relate x,y,z dimensions to 1 dimensions
        do nd = 1,3
          if ( trackingOptions%dimensionMask( nd ) .eq. 0 ) cycle
          select case(nd) 
            case (1)
              trackingOptions%idDim1 = nd
              write(outUnit,'(A)') 'RW displacements for X dimension.'
            case (2)
              trackingOptions%idDim1 = nd
              write(outUnit,'(A)') 'RW displacements for Y dimension.'
            case (3)
              trackingOptions%idDim1 = nd
              write(outUnit,'(A)') 'RW displacements for Z dimension.'
          end select   
          ! Use the first found
          exit
        end do
      ! 2D
      case(2)
        trackingOptions%twoDimensions = .true. ! TEMP
        write(outUnit,'(A)') 'RW displacements for 2 dimensions.'
        ! Relate x,y,z dimensions to 1,2 dimensions
        do nd = 1,3
          if ( trackingOptions%dimensionMask( nd ) .eq. 0 ) cycle
          currentDim = sum( trackingOptions%dimensionMask(1:nd) )
          select case(nd) 
            case (1)
              trackingOptions%idDim1 = n
              write(outUnit,'(A)') 'RW displacements for X dimension.'
            case (2)
              if ( currentDim .eq. 1 ) then 
                trackingOptions%idDim1 = nd
              else if ( currentDim .eq. 2 ) then
                trackingOptions%idDim2 = nd
              end if
              write(outUnit,'(A)') 'RW displacements for Y dimension.'
            case (3)
              trackingOptions%idDim2 = nd
              write(outUnit,'(A)') 'RW displacements for Z dimension.'
          end select   
        end do
      ! 3D
      case(3)
        trackingOptions%twoDimensions = .false.! TEMP
        write(outUnit,'(A)') 'RW displacements for 3 dimensions.'
    end select


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
    class(ModpathSimulationDataType), target     :: this
    character(len=200), intent(in)               :: icFile
    integer, intent(in)                          :: icUnit
    integer, intent(in)                          :: outUnit
    class(ModflowRectangularGridType),intent(in) :: grid
    doubleprecision, dimension(:)                :: porosity
    ! local
    integer :: isThisFileOpen = -1
    integer :: icol,istart,istop,n
    integer :: nc, nic, nd
    doubleprecision    :: r
    character(len=200) :: line
    
    integer :: nInitialConditions, nValidInitialConditions 
    integer :: initialConditionFormat
    integer, dimension(:), pointer :: dimensionMask
    integer, pointer               :: nDim

    type(ParticleGroupType),dimension(:),allocatable :: particleGroups
    type(ParticleGroupType),dimension(:),allocatable :: newParticleGroups

    doubleprecision, dimension(:), allocatable :: densityDistribution
    integer, dimension(:), allocatable :: cellsPerLayer
    doubleprecision :: particleMass, totalAccumulatedMass
    doubleprecision :: cellVolume,dX,dY,dZ,sX,sY,sZ
    integer :: soluteId
    doubleprecision :: initialReleaseTime

    character(len=24),dimension(1) :: aname
    data aname(1) /'            IC'/



  !! CLEAN AND VERIFY
  !character(len=24),dimension(5) :: aname
  !data aname(1) /'          BOUNDARY ARRAY'/
  !data aname(2) /'            DISPERSIVITY'/
  !data aname(3) /'            ICBOUND'/
  !data aname(4) /'            IC'/
  !data aname(5) /'          MEDIUMDISTANCE'/
  !integer :: tempAlphaUnit = 666
  !character(len=200) :: tempAlphaFile

  !! Initial conditions
  !integer :: nInitialConditions, particleCount, seqNumber, slocUnit
  !integer :: nValidInitialConditions
  !integer :: releaseOption, releaseTimeCount
  !doubleprecision :: initialReleaseTime, releaseInterval
  !doubleprecision,dimension(:),allocatable :: releaseTimes
  !doubleprecision :: frac, tinc
  !type(ParticleGroupType),dimension(:),allocatable :: particleGroups
  !type(ParticleGroupType),dimension(:),allocatable :: newParticleGroups
  !integer :: initialConditionFormat, massProcessingFormat
  !doubleprecision, dimension(:), allocatable :: densityDistribution
  !doubleprecision, dimension(:), allocatable :: rawNParticles
  !doubleprecision, dimension(:), allocatable :: nParticles
  !doubleprecision, dimension(:), allocatable :: cellVolumes
  !doubleprecision, dimension(:), allocatable :: delZ
  !! FOR DETERMINATION OF SUBDIVISIONS
  !doubleprecision, dimension(:), allocatable :: shapeFactorX
  !doubleprecision, dimension(:), allocatable :: shapeFactorY
  !doubleprecision, dimension(:), allocatable :: shapeFactorZ
  !doubleprecision, dimension(:), allocatable :: nParticlesX
  !doubleprecision, dimension(:), allocatable :: nParticlesY
  !doubleprecision, dimension(:), allocatable :: nParticlesZ

  !integer :: nic
  !integer :: newParticleGroupCount
  !doubleprecision :: minOneParticleDensity

  !! FROM READLOCATIONS3
  !type(ParticleGroupType) :: pGroup
  !integer :: totalParticleCount,templateCount,templateCellCount,nc,nr,nl,row,column,idmax
  !integer :: count,np,face,i,j,k,layerCount,rowCount,columnCount,       &    
  !           subCellCount,cell,offset,npcell
  !doubleprecision :: dx,dy,dz,x,y,z,faceCoord,rowCoord,columnCoord,dr,dc 
  !integer,dimension(:),allocatable :: templateSubDivisionTypes,         &
  !  templateCellNumbers,templateCellCounts, drape
  !integer,dimension(:,:),allocatable :: subDiv
  !integer,dimension(12) :: sdiv



    


    !--------------------------------------------------------------

    write(outUnit, *)
    write(outUnit, '(1x,a)') 'MODPATH-RW IC file data'
    write(outUnit, '(1x,a)') '-----------------------'

    ! Verify if unit is open 
    inquire( file=icFile, number=isThisFileOpen )
    if ( isThisFileOpen .lt. 0 ) then 
      ! No ic 
      write(outUnit,'(A)') 'No IC package for the RW simulation was specified.'
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

      ! Increase pgroup counter
      particleGroups(nic)%Group = this%ParticleGroupCount + nic

      ! Set release time for initial condition.
      ! It is an initial condition, then
      ! assumes release at referencetime
      initialReleaseTime = this%ReferenceTime
      call particleGroups(nic)%SetReleaseOption1(initialReleaseTime)

      ! Read id 
      read(icUnit, '(a)') particleGroups(nic)%Name

      ! Initial condition format
      ! 1: concentration
      read(icUnit, *) initialConditionFormat

      select case ( initialConditionFormat )
      ! Read initial condition as resident concentration (ML^-3)
      case (1) 
        ! Given a value for the mass of particles, 
        ! use flowModelData to compute cellvolume
        ! and a shape factor from which the number 
        ! of particles per cell is estimated

        if(allocated(densityDistribution)) deallocate(densityDistribution)
        allocate(densityDistribution(grid%CellCount))

        ! Read particles mass
        read(icUnit, *) particleMass

        if ( ( this%ParticlesMassOption .eq. 2 ) .or. & 
             ( this%SolutesOption .eq. 1 ) ) then 
          ! Read solute id
          ! Some validation
          ! Maybe read solutes before the IC's and 
          ! and validate against that information 
          read(icUnit, *) soluteId
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

        ! Validity of initial density distribution 
        if( all(abs(densityDistribution).eq.0d0) ) then
        write(outUnit,'(a,a,a)') 'Warning: initial condition ',&
                trim(particleGroups(nic)%Name),' has a distribution only with zeros'
          write(outUnit,'(a)') 'It will not create a particle group. Continue to the next.'
          ! Process the next one
          cycle
        end if 


        ! Loop over densityDistribution
        !  - Compute totalAccumulatedMass
        totalAccumulatedMass = 0d0
        do nc = 1, grid%CellCount
          if ( densityDistribution(nd) .eq. 0d0 ) cycle
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
          totalAccumulatedMass = densityDistribution(nc)*porosity(nc)*cellVolume
        end do 

        write(outUnit,'(/A,es18.9e3)') 'Total accumulated mass for initial condition = ', totalAccumulatedMass

      case default
        write(outUnit,*) 'Invalid initial condition kind ', initialConditionFormat, '. Stop.' 
        call ustop('Invalid initial condition kind')
      end select

    end do


      !!!!!!!!!!!!!!!!!!!!!!!!!!1
      !! CLEAN AND VERIFY 
      !!!!!!!!!!!!!!!!!!!!!!!!!!1

      !! Required for u3d
      !allocate(cellsPerLayer(grid%LayerCount))
      !do n = 1, grid%LayerCount
      !    cellsPerLayer(n) = grid%GetLayerCellCount(n)
      !end do

      !! Cell size
      !rowCount    = grid%RowCount
      !columnCount = grid%ColumnCount
      !layerCount  = grid%LayerCount


      !! Determine model dimensions (READIT!)
      !dimensionMask = 0
      !if ( grid%ColumnCount .gt. 1 ) dimensionMask(1) = 1 ! this dimension is active
      !if ( grid%RowCount .gt. 1 )    dimensionMask(2) = 1 ! this dimension is active
      !if ( grid%LayerCount .gt. 1 )  dimensionMask(3) = 1 ! this dimension is active
      !nDim = sum(dimensionMask)


      !! Write header to the listing file
      !write(outUnit, *)
      !write(outUnit, '(1x,a)') '------------------------------------'
      !write(outUnit, '(1x,a)') ' MODPATH-RW configuration file      '
      !write(outUnit, '(1x,a)') '------------------------------------'

      !! Process initial conditions 
      !read(inUnit, *) nInitialConditions
      !write(outUnit,'(/A,I5)') ' Number of initial conditions = ', nInitialConditions
      !nValidInitialConditions = 0

      !! What is this ?
      !particleCount = 0
      !slocUnit      = 0
      !seqNumber     = 0

      !! Extend particle groups
      !if(nInitialConditions .gt. 0) then

      !  ! Carrier of particle groups
      !  allocate(particleGroups(nInitialConditions))

      !  ! Loop over initial conditions
      !  do nic = 1, nInitialConditions

      !    ! Increase pgroup counter
      !    particleGroups(nic)%Group = simulationData%ParticleGroupCount + nic

      !    ! Set release time for initial condition.
      !    ! It is an initial condition, then
      !    ! assumes release at referencetime
      !    initialReleaseTime = simulationData%ReferenceTime
      !    call particleGroups(nic)%SetReleaseOption1(initialReleaseTime)

      !    ! Read id 
      !    read(inUnit, '(a)') particleGroups(nic)%Name

      !    ! Initial condition format
      !    ! 1: concentration
      !    read(inUnit, *) initialConditionFormat

      !    select case ( initialConditionFormat )

      !    ! Read initial condition as concentration times 
      !    ! porosity  (ML^-3)
      !    case (1) 

      !      ! Given a value for the mass of particles, 
      !      ! use flowModelData to compute cellvolume
      !      ! and a shape factor from which the number 
      !      ! of particles per cell is estimated

      !      ! Density/concentration arrays are expected 
      !      ! to be consistent with flow model grid.
      !      if(allocated(densityDistribution)) deallocate(densityDistribution)
      !      allocate(densityDistribution(grid%CellCount))
      !      if(allocated(rawNParticles)) deallocate(rawNParticles)
      !      allocate(rawNParticles(grid%CellCount))
      !      if(allocated(nParticles)) deallocate(nParticles)
      !      allocate(nParticles(grid%CellCount))
      !      if(allocated(cellVolumes)) deallocate(cellVolumes)
      !      allocate(cellVolumes(grid%CellCount))
      !      if(allocated(delZ)) deallocate(delZ)
      !      allocate(delZ(grid%CellCount))
      !      if(allocated(shapeFactorX)) deallocate(shapeFactorX)
      !      allocate(shapeFactorX(grid%CellCount))
      !      if(allocated(shapeFactorY)) deallocate(shapeFactorY)
      !      allocate(shapeFactorY(grid%CellCount))
      !      if(allocated(shapeFactorZ)) deallocate(shapeFactorZ)
      !      allocate(shapeFactorZ(grid%CellCount))
      !      if(allocated(nParticlesX)) deallocate(nParticlesX)
      !      allocate(nParticlesX(grid%CellCount))
      !      if(allocated(nParticlesY)) deallocate(nParticlesY)
      !      allocate(nParticlesY(grid%CellCount))
      !      if(allocated(nParticlesZ)) deallocate(nParticlesZ)
      !      allocate(nParticlesZ(grid%CellCount))

      !      ! Read particles mass
      !      read(inUnit, *) particleMass
      !
      !      if ( ( simulationData%ParticlesMassOption .eq. 2 ) .or. & 
      !           ( simulationData%SolutesOption .eq. 1 ) ) then 
      !        ! Read solute id
      !        ! Some validation 
      !        read(inUnit, *) soluteId
      !      end if 

      !      ! Read as density/concentration
      !      if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
      !        call u3ddblmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
      !          grid%ColumnCount, grid%CellCount, densityDistribution, ANAME(4))
      !      else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
      !        call u3ddblmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,  &
      !          densityDistribution, aname(4), cellsPerLayer)
      !      else
      !        write(outUnit,*) 'Invalid grid type specified when reading IC array ', & 
      !            particleGroups(nic)%Name, ' name.'
      !        write(outUnit,*) 'Stopping.'
      !        call ustop(' ')          
      !      end if
      !   
      !      ! Validity of initial density distribution 
      !      if( all(abs(densityDistribution).eq.0d0) ) then
      !        write(outUnit,*) 'Warning: initial condition ',&
      !            particleGroups(nic)%Name,' has a densityDistribution only with zeros'
      !        write(outUnit,*) 'It will not create a particle group. Continue to the next.'

      !        ! Process the next one
      !        cycle

      !      end if 

      !      ! Initialize
      !      nParticles   = 0
      !      cellVolumes  = 0d0
      !      shapeFactorX = 0
      !      shapeFactorY = 0
      !      shapeFactorZ = 0
      !      nParticlesX  = 0
      !      nParticlesY  = 0
      !      nParticlesZ  = 0

      !      ! Simple delZ 
      !      delZ = grid%Top-grid%Bottom
      !      ! LayerType if .eq. 1 convertible
      !      where( this%Grid%CellType .eq. 1 )
      !          ! ONLY IF HEADS .lt. TOP 
      !          delZ = flowModelData%Heads-grid%Bottom
      !      end where
      !      where ( delZ .le. 0d0 )
      !          delZ = 0d0 
      !      end where

      !      ! Compute cell volumes
      !      cellVolumes = 1d0
      !      do n =1, 3
      !        if ( dimensionMask(n) .eq. 1 ) then 
      !          select case( n ) 
      !            case(1)
      !              cellVolumes = cellVolumes*grid%DelX
      !            case(2)
      !              cellVolumes = cellVolumes*grid%DelY
      !            case(3)
      !              cellVolumes = cellVolumes*delZ
      !          end select
      !        end if 
      !      end do 

      !      ! Particles density 
      !      ! absolute value is required for the case that 
      !      ! densityDsitribution contains negative values
      !      rawNParticles = abs(densityDistribution)/particleMass

      !      ! nParticles
      !      ! notice rawNParticles = mass/volume/mass = 1/volume: particle density
      !      nParticles  = rawNParticles*cellVolumes

      !      ! The minimum measurable absolute concentration/density 
      !      minOneParticleDensity = minval(particleMass/cellVolumes)

      !      ! Allocate subdivision arrays
      !      templateCount     = 1 
      !      templateCellCount = grid%CellCount 
      !      allocate(subDiv(templateCount,12))
      !      allocate(templateSubDivisionTypes(templateCount))
      !      allocate(templateCellCounts(templateCount))
      !      allocate(drape(templateCount))
      !      allocate(templateCellNumbers(templateCellCount))

      !      subDiv(:,:) = 0
      !      drape       = 0 ! Should come from somewhere
      !      templateSubDivisionTypes(1) = 1
      !      templateCellCounts(1)       = grid%CellCount

      !      np     = 0
      !      npcell = 0
      !      offset = 0
      !      ! There is only one templateCount
      !      do n = 1, templateCount
      !        ! Determine subdivisions
      !        if ( dimensionMask(1) .eq. 1 ) then 
      !          ! Shape factors
      !          do m =1, grid%CellCount
      !              if ( cellVolumes(m) .le. 0d0 ) cycle
      !              shapeFactorX(m) = grid%DelX(m)/(cellVolumes(m)**(1d0/nDim))
      !          end do 
      !        end if 

      !        if ( dimensionMask(2) .eq. 1 ) then 
      !          ! Shape factors
      !          do m =1, grid%CellCount
      !              if ( cellVolumes(m) .le. 0d0 ) cycle
      !              shapeFactorY(m) = grid%DelY(m)/(cellVolumes(m)**(1d0/nDim))
      !          end do 
      !        end if 

      !        if ( dimensionMask(3) .eq. 1 ) then 
      !          ! Shape factors
      !          do m =1, grid%CellCount
      !              if ( cellVolumes(m) .le. 0d0 ) cycle
      !              shapeFactorZ(m) = delZ(m)/(cellVolumes(m)**(1d0/nDim))
      !          end do 
      !        end if 

      !        where ( nParticles .lt. 1d0 )
      !           nParticles = 0d0 
      !        end where 
      !        nParticlesX = shapeFactorX*( (nParticles)**(1d0/nDim) ) 
      !        nParticlesY = shapeFactorY*( (nParticles)**(1d0/nDim) )
      !        nParticlesZ = shapeFactorZ*( (nParticles)**(1d0/nDim) )

      !        ! Loop through cells and count the number of particles
      !        do cell = 1, templateCellCounts(n)
      !          if(ibound(cell) .ne. 0) then
      !            ! Could be a user defined threshold ?
      !            if ( nParticles(cell) .lt. 1d0 ) cycle
      !            ! Cell subdivisions
      !            subDiv(n,1) = int( nParticlesX(cell) ) + 1
      !            subDiv(n,2) = int( nParticlesY(cell) ) + 1
      !            subDiv(n,3) = int( nParticlesZ(cell) ) + 1
      !            npcell      = subDiv(n,1)*subDiv(n,2)*subDiv(n,3)
      !            np = np + npcell
      !          end if
      !        end do
      !        ! Increment the offset
      !        offset = offset + templateCellCounts(n)
      !      end do

      !      ! Calculate the total number of particles for all release time points.
      !      totalParticleCount = 0
      !      totalParticleCount = np*particleGroups(nic)%GetReleaseTimeCount()
      !      particleGroups(nic)%TotalParticleCount = totalParticleCount
      !      if ( totalParticleCount .eq. 0 ) then 
      !        write(outUnit,*) ' Warning: initial condition ',&
      !            particleGroups(nic)%Name,' has zero particles, it will skip this group.'

      !        ! Process the next one
      !        cycle

      !      end if 

      !      if(allocated(particleGroups(nic)%Particles)) deallocate(particleGroups(nic)%Particles)
      !      allocate(particleGroups(nic)%Particles(totalParticleCount))

      !      ! Update particles mass
      !      ! If for a given cell the value of density or concentration 
      !      ! is negative, then assign the number of particles
      !      ! and modify the sign of particles mass for that cell 
      !      ! with a negative sign 
      !      particleMass = sum( abs(densityDistribution)*cellVolumes )/totalParticleCount 
      !      write(outUnit,'(/A,es18.9e3)') ' Effective particleMass for initial condition = ', particleMass

      !      ! Assign to the particle group 
      !      particleGroups(nic)%Mass = particleMass

      !      ! Set the data for particles at the first release time point
      !      m = 0
      !      offset = 0
      !      if(templateCount .gt. 0) then
      !        do n = 1, templateCount
      !          do cell = 1, templateCellCounts(n)
      !            if(ibound(cell) .ne. 0) then
      !              if ( nParticles(cell) .lt. 1d0 ) cycle
      !              sdiv(1) = int( nParticlesX(cell) ) + 1
      !              sdiv(2) = int( nParticlesY(cell) ) + 1
      !              sdiv(3) = int( nParticlesZ(cell) ) + 1
      !              ! For the weird requirement where density
      !              ! might be negative...
      !              if ( densityDistribution(cell) .gt. 0 ) then 
      !                particleMass = abs( particleMass ) 
      !              else 
      !                particleMass = -1*abs( particleMass )
      !              end if 
      !              call CreateMassParticlesAsInternalArray(& 
      !                particleGroups(nic), cell, m, sdiv(1), sdiv(2), sdiv(3), & 
      !                drape(n), particleMass, particleGroups(nic)%GetReleaseTime(1) )
      !            end if
      !          end do
      !          ! Increment the offset
      !          offset = offset + templateCellCounts(n)
      !        end do
      !      end if
      !    
      !      ! Assign layer value to each particle
      !      idmax = 0
      !      seqNumber = 0
      !      do m = 1, totalParticleCount
      !          seqNumber = seqNumber + 1
      !          if(particleGroups(nic)%Particles(m)%Id .gt. idmax) idmax = particleGroups(nic)%Particles(m)%Id
      !          particleGroups(nic)%Particles(m)%Group = particleGroups(nic)%Group
      !          particleGroups(nic)%Particles(m)%SequenceNumber = seqNumber
      !          particleGroups(nic)%Particles(m)%InitialLayer =                                   &
      !            grid%GetLayer(particleGroups(nic)%Particles(m)%InitialCellNumber)
      !          particleGroups(nic)%Particles(m)%Layer =                                          &
      !            grid%GetLayer(particleGroups(nic)%Particles(m)%CellNumber)
      !      end do
  

      !      ! Deallocate temporary arrays
      !      deallocate(subDiv)
      !      deallocate(templateSubDivisionTypes)
      !      deallocate(drape)
      !      deallocate(templateCellCounts)
      !      deallocate(templateCellNumbers)

      !    case default
      !      write(outUnit,*) ' Invalid initial condition kind ', initialConditionFormat 
      !      write(outUnit,*) ' Stop.'
      !      call ustop(' Error: Invalid initial condition kind ')
      !    end select

      !    ! Increment valid couter
      !    nValidInitialConditions = nValidInitialConditions + 1 

      !    if ( simulationData%ParticlesMassOption .eq. 2 ) then 
      !      ! Assign the solute id 
      !      particleGroups(nic)%Solute = soluteId
      !    end if 

      !    ! Report number of particles
      !    write(outUnit, '(a,i4,a,i10,a)') ' Initial condition ', nic, ' contains ',   &
      !      particleGroups(nic)%TotalParticleCount, ' particles.'
      !    particleCount = particleCount + particleGroups(nic)%TotalParticleCount

      !  end do ! loop over initial conditions 
      !  write(outUnit, '(a,i10)') ' Total number of particles on initial conditions = ', particleCount
      !  write(outUnit, *)


      !  ! Extend simulationdata to include these particle groups
      !  if ( nValidInitialConditions .gt. 0 ) then 
      !    newParticleGroupCount = simulationData%ParticleGroupCount + nValidInitialConditions
      !    allocate(newParticleGroups(newParticleGroupCount))
      !    ! If some particle groups existed previously
      !    if( simulationData%ParticleGroupCount .gt. 0 ) then 
      !      do n = 1, simulationData%ParticleGroupCount
      !        newParticleGroups(n) = simulationData%ParticleGroups(n)
      !      end do
      !    end if 
      !    ncount = 0
      !    do n = 1, nInitialConditions
      !      if ( particleGroups(n)%TotalParticleCount .eq. 0 ) cycle
      !      ncount = ncount + 1 
      !      newParticleGroups(ncount+simulationData%ParticleGroupCount) = particleGroups(n)
      !    end do 
      !    if( simulationData%ParticleGroupCount .gt. 0 ) then 
      !      call move_alloc( newParticleGroups, simulationData%ParticleGroups )
      !      simulationData%ParticleGroupCount = newParticleGroupCount
      !      simulationData%TotalParticleCount = simulationData%TotalParticleCount + particleCount
      !    else
      !      simulationData%ParticleGroupCount = newParticleGroupCount
      !      simulationData%TotalParticleCount = particleCount
      !      allocate(simulationData%ParticleGroups(simulationData%ParticleGroupCount))
      !      call move_alloc( newParticleGroups, simulationData%ParticleGroups )
      !    end if
      !  end if

      !  ! Will process the next at some point
      !  deallocate( particleGroups )
     
      !end if ! if nInitialConditions .gt. 0

















    ! Close rwopts data file
    close( icUnit )


  end subroutine pr_ReadICData








end module ModpathSimulationDataModule
