module ParticleTrackingOptionsModule
  use ObservationModule, only : ObservationType
  use BoundaryConditionsModule, only : PrescribedType
  implicit none
  
! Set default access status to private
  private

  type,public :: ParticleTrackingOptionsType
    logical :: DebugMode = .false.
    logical :: BackwardTracking = .false.
    logical :: CreateTrackingLog = .false.
    logical :: StopAtWeakSinks = .false.
    logical :: StopAtWeakSources = .false.
    logical :: ExtendSteadyState = .true.
    logical :: SpecifyStoppingTime = .false.
    logical :: SpecifyStoppingZone = .false.
    doubleprecision :: Stoptime = 1.0d+30
    integer :: StopZone = 0

    ! RWPT
    logical                            :: RandomWalkParticleTracking = .false.
    integer                            :: timeStepKind
    doubleprecision, dimension(2)      :: timeStepParameters = 0
    integer                            :: advectionKind
    integer, dimension(3)              :: dimensionMask = 1
    integer                            :: nDim = 3
    integer                            :: idDim1, idDim2
    integer, dimension(:), allocatable :: dimensions
    integer                            :: dispersionModel ! Needed to indicate dispersion function
    logical                            :: twoDimensions = .false. ! DEPRECATE

    ! GPKDE
    logical                       :: GPKDEReconstruction = .false.
    doubleprecision, dimension(3) :: gpkdeDomainSize
    doubleprecision, dimension(3) :: gpkdeBinSize
    doubleprecision, dimension(3) :: gpkdeDomainOrigin
    integer                       :: gpkdeNOptLoops 
    character(len=200)            :: gpkdeOutputFile
    integer                       :: gpkdeOutputUnit = 125
    logical                       :: gpkdeKernelDatabase
    logical                       :: gpkdeAsConcentration = .false.
    doubleprecision               :: gpkdeBinVolume
    doubleprecision               :: gpkdeScalingFactor
    doubleprecision, dimension(3) :: gpkdeKDBParams  ! minHLambda, deltaHLambda, maxHLambda

    ! TimeseriesOutputOption
    logical                       :: skipTimeseriesWriter = .false.

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

    !! Prescribed concentration boundaries
    !type( PrescribedType ), allocatable, dimension(:) :: PrescribedBoundaries
    !logical, allocatable, dimension(:) :: isPrescribed
    !integer, allocatable, dimension(:) :: idPrescribed
    !integer :: nPrescribed
    !integer :: npcs

    ! DEPRECATION WARNNING
    ! NONLINEAR DISPERSION RWPT (TEMP)
    doubleprecision :: betaTrans, betaLong
    doubleprecision :: mediumDistance, mediumDelta
    doubleprecision :: Dmol = 0d0
    ! DEPRECATION WARNING

  contains

     procedure :: Reset=>pr_Reset

     ! Observations
     procedure :: InitializeObservations=>pr_InitializeObservations 
     procedure :: IdObservationCell=>pr_IdObservationCell
     ! Prescribed boundaries
     !procedure :: InitializePrescribedBoundaries=>pr_InitializePrescribedBoundaries 
  end type
  
contains


  subroutine pr_Reset(this)
      !-----------------------------------------------------
      ! Specifications
      !-----------------------------------------------------
      implicit none
      class(ParticleTrackingOptionsType) :: this
      !-----------------------------------------------------

      ! This function is incomplete

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
      !
      !-----------------------------------------------------
      ! Specifications
      !-----------------------------------------------------
      implicit none
      class(ParticleTrackingOptionsType) :: this
      ! input
      integer, intent(in) :: nObservations
      integer :: n 
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


  !! TEMPORARY, find a better place
  !subroutine pr_InitializePrescribedBoundaries( this, nPrescribed )
  !  !-----------------------------------------------------
  !  ! Initialize observation cells data
  !  !-----------------------------------------------------
  !  ! Specifications
  !  !-----------------------------------------------------
  !  implicit none
  !  class(ParticleTrackingOptionsType) :: this
  !  ! input
  !  integer, intent(in) :: nPrescribed
  !  integer :: n 
  !  !-------------------------

  !  ! Allocate observation files arrays
  !  this%nPrescribed = nPrescribed

  !  allocate(this%PrescribedBoundaries(nPrescribed))

  !end subroutine pr_InitializePrescribedBoundaries


end module ParticleTrackingOptionsModule
