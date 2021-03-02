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

    !RWPT
    logical                       :: RandomWalkParticleTracking = .false.
    doubleprecision               :: alphaL, alphaT, Dmol = 0
    integer                       :: timeStepKind
    doubleprecision, dimension(2) :: timeStepParameters = 0
    integer                       :: advectionKind
    logical                       :: twoDimensions = .false.
  end type
  
contains

end module ParticleTrackingOptionsModule
