module ParticleTrackingOptionsModule
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
    doubleprecision               :: Dmol = 0
    !doubleprecision               :: alphaL, alphaT, Dmol = 0
    integer                       :: timeStepKind
    doubleprecision, dimension(2) :: timeStepParameters = 0
    integer                       :: advectionKind
    logical                       :: twoDimensions = .false.

    ! GPKDE
    logical                       :: GPKDEReconstruction = .false.
    doubleprecision, dimension(3) :: gpkdeDomainSize
    doubleprecision, dimension(3) :: gpkdeBinSize
    doubleprecision, dimension(3) :: gpkdeDomainOrigin
    integer                       :: gpkdeNOptLoops 
    character(len=200)            :: gpkdeOutputFile
    integer                       :: gpkdeOutputUnit = 125
    logical                       :: gpkdeKernelDatabase
    doubleprecision, dimension(3) :: gpkdeKDBParams


    ! TEMPORARY NONLINEAR RWPT
    integer         :: dispersionModel
    doubleprecision :: betaTrans, betaLong
    doubleprecision :: mediumDistance, mediumDelta


    ! OBS
    integer :: nObservations
    logical :: observationSimulation = .false.
    integer, allocatable, dimension(:) :: observationCells
    integer, allocatable, dimension(:) :: observationUnits
    character(len=200), allocatable, dimension(:) :: observationFiles

  contains
     procedure :: Reset=>pr_Reset
     procedure :: InitializeObservations=>pr_InitializeObservations 
     !procedure :: IsObservationCell=>pr_IsObservationCell
     procedure :: IdObservationCell=>pr_IdObservationCell
  end type
  
contains


  subroutine pr_Reset(this)

      !-----------------------
      ! Deallocate observation cells data
      !
      !--------------
      ! Specifications
      !--------
      implicit none
      class(ParticleTrackingOptionsType) :: this
      !------------------------------------

      this%nObservations         = 0
      this%observationSimulation = .false.

      if(allocated(this%observationCells)) deallocate(this%observationCells)
      if(allocated(this%observationUnits)) deallocate(this%observationUnits)
      if(allocated(this%observationFiles)) deallocate(this%observationFiles)

      return 

  end subroutine pr_Reset


  subroutine pr_InitializeObservations( this, nObservations )
      !-----------------------
      ! Initialize observation cells data
      !
      !--------------
      ! Specifications
      !--------
      implicit none
      class(ParticleTrackingOptionsType) :: this
      ! input
      integer, intent(in) :: nObservations
      integer :: n 
      !-------------------------

      ! Assign observation cells data      
      this%nObservations = nObservations
      this%observationSimulation = .true.
      allocate(this%observationCells(nObservations))
      allocate(this%observationUnits(nObservations))
      allocate(this%observationFiles(nObservations))

      !do n = 1, nObservations
      !    this%observationCells( n ) = 
      !end do


      return 

  end subroutine pr_InitializeObservations

  
  function pr_IdObservationCell( this, cellNumber ) result( output )
      !-----------------------
      ! Determine if cellNumber in the list of observation cells 
      ! 
      ! return the index
      !--------------
      ! Specifications
      !--------
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

      ! See if cellNumber is on the list of observation cells     
      do n=1, this%nObservations
          if ( this%observationCells( n ) .eq. cellNumber ) then 
              output = n
              return
          end if 
      end do 

      return 

  end function pr_IdObservationCell 


end module ParticleTrackingOptionsModule
