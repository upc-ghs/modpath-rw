module ObservationModule
  use ParticleModule,only : ParticleType
  use ParticleCoordinateModule,only : ParticleCoordinateType
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none
  
  ! Set default access status to private
  private

  ! Observation type
  type, public :: ObservationType
    integer, allocatable, dimension(:) :: cells
    integer                            :: nCells
    integer, allocatable, dimension(:) :: nRecordsCell
    integer                            :: nAuxRecords = 0
    integer                            :: id
    character(len=20)                  :: stringid
    integer                            :: style
    character(len=20)                  :: stylestringid
    integer                            :: cellOption
    integer                            :: timeOption
    integer                            :: outputOption
    integer                            :: postprocessOption
    logical                            :: doPostprocess
    integer                            :: outputUnit
    character(len=200)                 :: outputFileName
    integer                            :: recOutputUnit
    character(len=200)                 :: recOutputFileName
    integer                            :: auxOutputUnit
    character(len=200)                 :: auxOutputFileName
    doubleprecision                    :: cummSinkFlow
    doubleprecision, dimension(:), allocatable :: cummSinkFlowSeries
    ! Reconstruction params
    integer :: effectiveWeightFormat    = 0
    doubleprecision :: binSizeFactor    = 3d0
    integer :: initialSmoothingFormat   = 0
    integer :: nOptLoops                = 10
    doubleprecision :: errorConvergence = 0.01d0 
    integer :: histogramBinFormat       = 1
    doubleprecision :: histogramBin     = 0d0
    integer :: histogramOptions         = 0
    integer :: reconstructionOptions    = 0
    doubleprecision :: timeStepOut      = 0d0
    logical :: adaptGridToCoords        = .false.
    logical :: doInterp                 = .false.
    doubleprecision :: binSizeFraction  = 0.75d0 
    doubleprecision, dimension(:), allocatable :: obsSeries
    doubleprecision, dimension(:), allocatable :: timeSeries
    integer :: timePointCount
    type(SeriesType), allocatable, dimension(:) :: series

    ! Record writers
    procedure(SinkObsWriter), pass, pointer :: WriteSinkObs=>null()
    procedure(ResidentObsWriter), pass, pointer :: WriteResidentObs=>null()

  end type

  ! Series type
  type, public :: SeriesType
    doubleprecision, dimension(:)  , allocatable :: timeSeries
    doubleprecision, dimension(:,:), allocatable :: dataSeries ! for the same time may have multiple data columns
  end type


  ! Interfaces
  abstract interface
    subroutine SinkObsWriter( this, timeStep, timePointIndex, particle, & 
                                            soluteID, flowRate, outUnit )
      !-------------------------------------------------------------
      import ObservationType
      import ParticleType
      !-------------------------------------------------------------
      class(ObservationType) :: this
      ! input
      integer,intent(in)              :: timeStep, timePointIndex
      type(ParticleType),intent(in)   :: particle
      integer, intent(in)             :: soluteID
      doubleprecision, intent(in)     :: flowRate
      integer, intent(in)             :: outUnit
      !----------------------------------------------
    end subroutine SinkObsWriter

    subroutine ResidentObsWriter( this, timeStep, timePointIndex, particle, &
                      pCoord, soluteID, rFactor, waterVolume, outUnit )
      !-------------------------------------------------------------
      import ObservationType
      import ParticleType
      import ParticleCoordinateType
      !-------------------------------------------------------------
      class(ObservationType) :: this
      integer,intent(in)                      :: timeStep, timePointIndex
      type(ParticleType),intent(in)           :: particle
      type(ParticleCoordinateType),intent(in) :: pCoord
      integer, intent(in)                     :: soluteID
      doubleprecision, intent(in)             :: rFactor, waterVolume 
      integer, intent(in)                     :: outUnit
    end subroutine ResidentObsWriter

  end interface

  
  ! Open some methods
  public WriteSinkObsRecord
  public WriteSinkObsRecordBinary
  public WriteResidentObsRecord
  public WriteResidentObsRecordBinary


contains


  subroutine WriteSinkObsRecord( this, timeStep, timePointIndex, particle, &
                                               soluteID, flowRate, outUnit )
  !---------------------------------------------------------------------------------
  ! Write observation sink cell record
  !---------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  class(ObservationType) :: this
  ! input
  integer,intent(in)              :: timeStep, timePointIndex
  type(ParticleType),intent(in)   :: particle
  integer, intent(in)             :: soluteID
  doubleprecision, intent(in)     :: flowRate
  integer, intent(in)             :: outUnit
  !---------------------------------------------------------------------------------

    ! These could make use of the properties in the obs object

    write(outUnit, '(2I8,es18.9e3,i10,es18.9e3,2i5,2i10,es18.9e3)')                &
      timePointIndex, timeStep, particle%TrackingTime, particle%ID, particle%Mass, & 
        particle%Group, soluteID, particle%CellNumber, particle%Layer, flowRate


  end subroutine WriteSinkObsRecord


  subroutine WriteSinkObsRecordBinary( this, timeStep, timePointIndex, particle, & 
                                                     soluteID, flowRate, outUnit )
  !---------------------------------------------------------------------------------
  ! Write observation sink cell record
  !---------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  class(ObservationType) :: this
  ! input
  integer,intent(in)              :: timeStep, timePointIndex
  type(ParticleType),intent(in)   :: particle
  integer, intent(in)             :: soluteID
  doubleprecision, intent(in)     :: flowRate
  integer, intent(in)             :: outUnit
  !---------------------------------------------------------------------------------

    ! These could make use of the properties in the obs object

    write(outUnit) &
      timePointIndex, timeStep, particle%TrackingTime, particle%ID, particle%Mass, & 
        particle%Group, soluteID, particle%CellNumber, particle%Layer, flowRate


  end subroutine WriteSinkObsRecordBinary


  subroutine WriteResidentObsRecord( this, timeStep, timePointIndex, particle, pCoord, &
                                               soluteID, rFactor, waterVolume, outUnit )
  !--------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  class(ObservationType) :: this
  integer,intent(in)                      :: timeStep, timePointIndex
  type(ParticleType),intent(in)           :: particle
  type(ParticleCoordinateType),intent(in) :: pCoord
  integer, intent(in)                     :: soluteID
  doubleprecision, intent(in)             :: rFactor, waterVolume 
  integer, intent(in)                     :: outUnit
  doubleprecision :: modelX, modelY
  !---------------------------------------------------------------------------------
  
    modelX = pCoord%GlobalX
    modelY = pCoord%GlobalY
    
    write(outUnit, '(2I8,es18.9e3,i10,es18.9e3,2i5,2i10,5es18.9e3)')             &
      timePointIndex, timeStep, pCoord%TrackingTime, particle%ID, particle%Mass, & 
                      particle%Group, soluteID, pCoord%CellNumber, pCoord%Layer, &
                           rFactor, waterVolume, modelX, modelY, pCoord%GlobalZ


  end subroutine WriteResidentObsRecord


  subroutine WriteResidentObsRecordBinary( this, timeStep, timePointIndex, particle, pCoord, &
                                                     soluteID, rFactor, waterVolume, outUnit )
  !--------------------------------------------------------------------------------
  ! Specifications
  !---------------------------------------------------------------------------------
  implicit none
  class(ObservationType) :: this
  integer,intent(in)                      :: timeStep, timePointIndex
  type(ParticleType),intent(in)           :: particle
  type(ParticleCoordinateType),intent(in) :: pCoord
  integer, intent(in)                     :: soluteID
  doubleprecision, intent(in)             :: rFactor, waterVolume 
  integer, intent(in)                     :: outUnit
  doubleprecision :: modelX, modelY
  !---------------------------------------------------------------------------------
  
    modelX = pCoord%GlobalX
    modelY = pCoord%GlobalY
    
    write(outUnit) &
      timePointIndex, timeStep, pCoord%TrackingTime, particle%ID, particle%Mass, & 
                      particle%Group, soluteID, pCoord%CellNumber, pCoord%Layer, &
                           rFactor, waterVolume, modelX, modelY, pCoord%GlobalZ


  end subroutine WriteResidentObsRecordBinary



end module ObservationModule
