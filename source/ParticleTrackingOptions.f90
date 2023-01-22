module ParticleTrackingOptionsModule
  use ObservationModule, only : ObservationType
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
    logical                       :: RandomWalkParticleTracking = .false.
    integer                       :: timeStepKind
    doubleprecision, dimension(2) :: timeStepParameters = 0
    integer                       :: advectionKind

    ! RW displacements, defaults to 3D
    integer, dimension(3)         :: dimensionMask = 1
    integer                       :: nDim = 3
    logical                       :: twoDimensions = .false.
    integer                       :: idDim1, idDim2


    integer                       :: dispersionModel
    doubleprecision               :: Dmol = 0



    ! NONLINEAR DISPERSION RWPT (TEMP)
    doubleprecision :: betaTrans, betaLong
    doubleprecision :: mediumDistance, mediumDelta



    ! GPKDE
    logical                       :: GPKDEReconstruction = .false.
    doubleprecision, dimension(3) :: gpkdeDomainSize
    doubleprecision, dimension(3) :: gpkdeBinSize
    doubleprecision, dimension(3) :: gpkdeDomainOrigin
    integer                       :: gpkdeNOptLoops 
    character(len=200)            :: gpkdeOutputFile
    integer                       :: gpkdeOutputUnit = 125
    logical                       :: gpkdeKernelDatabase
    doubleprecision, dimension(3) :: gpkdeKDBParams  ! minHLambda, deltaHLambda, maxHLambda

    ! TIMESERIES
    logical                       :: skipTimeseriesWriter = .false.

    ! OBSERVATION CELLS
    integer :: nObservations
    logical :: observationSimulation = .false.
    integer, allocatable, dimension(:) :: observationCells
    integer, allocatable, dimension(:) :: observationUnits
    character(len=200), allocatable, dimension(:) :: observationFiles
    integer, allocatable, dimension(:) :: obsRecordCounts
    logical, allocatable, dimension(:) :: isObservation
    integer, allocatable, dimension(:) :: idObservation
    logical :: anySinkObservation = .false.
    type( ObservationType ), allocatable, dimension(:) :: Observations


  contains

     procedure :: Reset=>pr_Reset
     procedure :: InitializeObservations=>pr_InitializeObservations 
     procedure :: IdObservationCell=>pr_IdObservationCell

  end type
  
contains


  subroutine pr_Reset(this)
      !-----------------------------------------------------
      ! Specifications
      !-----------------------------------------------------
      implicit none
      class(ParticleTrackingOptionsType) :: this
      !-----------------------------------------------------

      ! Deallocate observation cells 
      this%nObservations         = 0
      this%observationSimulation = .false.
      this%anySinkObservation    = .false.
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
      this%observationSimulation = .true.
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
