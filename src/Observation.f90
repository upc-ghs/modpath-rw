module ObservationModule
  implicit none
  
  ! Set default access status to private
  private

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
  end type


  type, public :: SeriesType
    doubleprecision, dimension(:)  , allocatable :: timeSeries
    doubleprecision, dimension(:,:), allocatable :: dataSeries ! for the same time may have multiple data columns
  end type


contains

  ! Some subroutines


end module ObservationModule
