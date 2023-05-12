module ParticleTrackingOptionsModule
  use PrecisionModule, only : fp 
  use ObservationModule, only : ObservationType
  implicit none
  
  ! Set default access status to private
  private

  type,public :: ParticleTrackingOptionsType
    logical  :: DebugMode = .false.
    logical  :: BackwardTracking = .false.
    logical  :: CreateTrackingLog = .false.
    logical  :: StopAtWeakSinks = .false.
    logical  :: StopAtWeakSources = .false.
    logical  :: ExtendSteadyState = .true.
    logical  :: SpecifyStoppingTime = .false.
    logical  :: SpecifyStoppingZone = .false.
    real(fp) :: Stoptime = 1.0e+30_fp
    integer  :: StopZone = 0

    ! RWPT
    logical                            :: RandomWalkParticleTracking = .false.
    integer                            :: timeStepKind
    real(fp), dimension(2)             :: timeStepParameters = 0.0_fp
    integer                            :: advectionKind
    integer, dimension(3)              :: dimensionMask = 1
    integer                            :: nDim = 3
    integer                            :: idDim1, idDim2
    integer, dimension(:), allocatable :: dimensions
    integer                            :: randomGenFunction

    ! GPKDE
    logical                :: GPKDEReconstruction = .false.
    real(fp), dimension(3) :: gpkdeDomainSize
    real(fp), dimension(3) :: gpkdeBinSize
    real(fp), dimension(3) :: gpkdeDomainOrigin
    integer                :: gpkdeNOptLoops 
    character(len=200)     :: gpkdeOutputFile
    integer                :: gpkdeOutputUnit = 125
    logical                :: gpkdeKernelDatabase
    logical                :: gpkdeAsConcentration = .false.
    real(fp)               :: gpkdeBinVolume
    real(fp)               :: gpkdeScalingFactor
    real(fp), dimension(3) :: gpkdeKDBParams  ! minHLambda, deltaHLambda, maxHLambda

    ! TimeseriesOutputOption
    logical                :: skipTimeseriesWriter = .false.

    ! Observation Cells 
    integer :: nObservations
    logical :: anyObservation = .false.
    logical :: anyResObservation = .false.
    logical :: anySinkObservation = .false.
    integer, allocatable, dimension(:) :: observationCells
    integer, allocatable, dimension(:) :: observationUnits
    character(len=200), allocatable, dimension(:) :: observationFiles
    integer, allocatable, dimension(:) :: obsRecordCounts
    logical, allocatable, dimension(:) :: isObservation
    integer, allocatable, dimension(:) :: idObservation
    type( ObservationType ), allocatable, dimension(:) :: Observations

    ! DEPRECATION WARNNING
    ! NONLINEAR DISPERSION RWPT (TEMP)
    real(fp) :: betaTrans, betaLong
    real(fp) :: mediumDistance, mediumDelta
    real(fp) :: Dmol = 0.0_fp
    ! DEPRECATION WARNING

  contains
     procedure :: Reset=>pr_Reset
     ! Observations
     procedure :: InitializeObservations=>pr_InitializeObservations 
     procedure :: IdObservationCell=>pr_IdObservationCell
  end type
  
contains


  subroutine pr_Reset(this)
    !-----------------------------------------------------
    ! 
    !-----------------------------------------------------
    ! Specifications
    !-----------------------------------------------------
    implicit none
    class(ParticleTrackingOptionsType) :: this
    !-----------------------------------------------------
    ! Deallocate observation cells 
    this%nObservations         = 0
    this%anyObservation        = .false.
    this%anySinkObservation    = .false.
    this%anyResObservation     = .false.
    if(allocated(this%observationCells)) deallocate(this%observationCells)
    if(allocated(this%observationUnits)) deallocate(this%observationUnits)
    if(allocated(this%observationFiles)) deallocate(this%observationFiles)
    if(allocated(this%obsRecordCounts)) deallocate(this%obsRecordCounts)
    if(allocated(this%Observations)) deallocate(this%Observations)
    if(allocated(this%isObservation)) deallocate(this%isObservation)
    if(allocated(this%idObservation)) deallocate(this%idObservation)
  end subroutine pr_Reset


  subroutine pr_InitializeObservations( this, nObservations )
    !-----------------------------------------------------
    ! Initialize observation cells data
    !-----------------------------------------------------
    ! Specifications
    !-----------------------------------------------------
    implicit none
    class(ParticleTrackingOptionsType) :: this
    ! input
    integer, intent(in) :: nObservations
    !-------------------------

    ! Allocate observation files arrays
    this%nObservations = nObservations
    allocate(this%observationCells(nObservations))
    allocate(this%observationUnits(nObservations))
    allocate(this%observationFiles(nObservations))
    allocate(this%obsRecordCounts(nObservations))
    allocate(this%Observations(nObservations))
  end subroutine pr_InitializeObservations

  
  function pr_IdObservationCell( this, cellNumber ) result( output )
    !-------------------------------------------------------------
    ! Determine if cellNumber in the list of observation cells
    ! 
    ! return the index
    !-------------------------------------------------------------
    ! Specifications
    !-------------------------------------------------------------
    implicit none
    class(ParticleTrackingOptionsType) :: this
    ! input
    integer, intent(in) :: cellNumber
    ! output
    !logical :: output = .false.
    integer :: output 
    ! local
    integer :: n
    !-------------------------

    output = -999

    ! Check if cellNumber is on the list of observation cells     
    do n=1, this%nObservations
      if ( this%observationCells( n ) .eq. cellNumber ) then 
        output = n
        return
      end if 
    end do 

    return 
  end function pr_IdObservationCell 


end module ParticleTrackingOptionsModule
