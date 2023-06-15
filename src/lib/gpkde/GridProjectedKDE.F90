module GridProjectedKDEModule
  !----------------------------------------------------------------------------
  ! Main module providing Grid Projected Kernel Density Estimation
  !----------------------------------------------------------------------------
  use PrecisionModule, only : fp 
  use ConstantsModule, only : fZERO, fONE, fTWO,    & 
                              fTHREE, fFOUR, fFIVE, & 
                              fSIX, fEIGHT, pi, sqrtEightPi
  use HistogramModule, only : HistogramType
  use KernelMultiGaussianModule, only : InitializeKernelDimensions, &
                                           KernelMultiGaussianType, &
                                       KernelSecondDerivativeXType, &
                                       KernelSecondDerivativeYType, &
                                       KernelSecondDerivativeZType, &
                                                        KernelType
  use GridCellModule, only : GridCellType
  implicit none
  !----------------------------------------------------------------------------

  ! Default parameters
  integer  , parameter :: defaultKernelRange                     = 3
  integer  , parameter :: defaultKernelSDRange                   = 4
  integer  , parameter :: defaultNOptLoops                       = 10
  logical  , parameter :: defaultDatabaseOptimization            = .true.
  logical  , parameter :: defaultLogKernelDatabase               = .true.
  real(fp) , parameter :: defaultMaxHOverDelta                   = 15.0_fp
  real(fp) , parameter :: defaultMinHOverDelta                   = 0.3_fp
  real(fp) , parameter :: defaultDeltaHOverDelta                 = 0.3_fp
  real(fp) , parameter :: defaultRelativeErrorConvergence        = 0.02_fp
  real(fp) , parameter :: defaultRelaxedRelativeErrorConvergence = 0.05_fp
  integer  , parameter :: defaultInitialSmoothingSelection       = 0
  real(fp) , parameter :: defaultInitialSmoothingFactor          = 5.0_fp
  real(fp) , parameter :: defaultInitialSigmaFactor              = 2.0_fp
  logical  , parameter :: defaultAdaptGridToCoords               = .false. 
  real(fp) , parameter :: defaultMinRelativeRoughness            = 1.0e-3_fp
  real(fp) , parameter :: defaultMinLimitRoughness               = 1.0e-3_fp
  integer  , parameter :: defaultMinRoughnessFormat              = 3
  integer  , parameter :: defaultEffectiveWeightFormat           = 0
  integer  , parameter :: defaultBoundKernelSizeFormat           = 0
  real(fp) , parameter :: defaultIsotropicThreshold              = 0.85_fp
  logical  , parameter :: defaultUseGlobalSmoothing              = .false.
  real(fp) , parameter :: defaultMinSizeFactor                   = 1.2_fp
  real(fp) , parameter :: defaultMaxSizeFactor                   = 0.5_fp 
  real(fp) , parameter :: defaultBorderFraction                  = 0.05_fp
  real(fp) , parameter :: defaultMaxSigmaGrowth                  = 1.5_fp
  character(len=*), parameter :: defaultOutputFileName           = 'gpkde.out'

  ! Module variables defined after initialization
  integer               :: nDim
  integer, dimension(3) :: dimensionMask = (/1,1,1/)
  real(fp)              :: fNDim
  real(fp)              :: oneOverNDimPlusFour
  real(fp)              :: minusOneOverNDimPlusSix
  real(fp)              :: onePlusNDimQuarter

  ! Set default access to private
  private

  ! Arrays, related to active bins
  real(fp) , dimension(:,:,:), allocatable, target   :: densityGrid
  real(fp) , dimension(:,:)  , allocatable           :: kernelSmoothing
  real(fp) , dimension(:)    , allocatable           :: kernelSmoothingScale
  real(fp) , dimension(:,:)  , allocatable           :: kernelSmoothingShape
  real(fp) , dimension(:)    , allocatable           :: kernelSigmaSupportScale
  real(fp) , dimension(:,:)  , allocatable           :: curvatureBandwidth
  real(fp) , dimension(:)    , allocatable           :: densityEstimateArray 
  real(fp) , dimension(:)    , allocatable           :: nEstimateArray
  real(fp) , dimension(:)    , allocatable, target   :: roughnessXXArray
  real(fp) , dimension(:)    , allocatable, target   :: roughnessYYArray
  real(fp) , dimension(:)    , allocatable, target   :: roughnessZZArray
  real(fp) , dimension(:)    , allocatable           :: netRoughnessArray
  type(GridCellType), dimension(:), allocatable, target :: activeGridCellsMod
  
  ! Main object
  type, public :: GridProjectedKDEType
  
    ! Histogram and kernel databases 
    type( HistogramType )                                              :: histogram
    type( KernelMultiGaussianType )    , dimension(:,:,:), allocatable :: kernelDatabase
    type( KernelMultiGaussianType )    , dimension(:,:)  , allocatable :: kernelDatabaseFlat
    type( KernelSecondDerivativeXType ), dimension(:)    , allocatable :: kernelSDXDatabase
    type( KernelSecondDerivativeYType ), dimension(:)    , allocatable :: kernelSDYDatabase
    type( KernelSecondDerivativeZType ), dimension(:)    , allocatable :: kernelSDZDatabase

    ! Dimensionality
    integer, dimension(3)                      :: dimensionMask
    integer, dimension(:), allocatable         :: dimensions
    ! For 1d-2d mapping
    integer                                    :: idDim1, idDim2
    class( KernelType ), dimension(:), pointer :: kernelSDDatabase1
    class( KernelType ), dimension(:), pointer :: kernelSDDatabase2
    real(fp), dimension(:), pointer            :: roughness11Array
    real(fp), dimension(:), pointer            :: roughness22Array

    ! Grid properties 
    real(fp), dimension(3) :: binSize
    real(fp), dimension(3) :: domainSize
    real(fp), dimension(3) :: domainOrigin
    integer , dimension(3) :: domainGridSize
    integer , dimension(3) :: deltaBinsOrigin
    integer , dimension(3) :: nBins
    logical                :: adaptGridToCoords
    real(fp)               :: borderFraction
    logical                :: slicedReconstruction 
    integer                :: slicedDimension

    ! Variables
    real(fp), dimension(:,:,:), pointer :: densityEstimateGrid => null()
    real(fp), dimension(:,:,:), pointer :: histogramDensity    => null()
    
    ! Kernel database params 
    real(fp), dimension(3) :: deltaHOverDelta
    real(fp), dimension(3) :: minHOverDelta
    real(fp), dimension(3) :: maxHOverDelta
    integer , dimension(3) :: nDeltaHOverDelta ! Computed at kernel databases
    logical                :: logKernelDatabase 
    logical                :: databaseOptimization 
    logical                :: flatKernelDatabase ! Deprecate ?
    
    ! Optimization
    integer                :: nOptimizationLoops
    real(fp)               :: densityRelativeConvergence
    real(fp)               :: minLimitRoughness = fZERO
    real(fp)               :: minRelativeRoughness
    real(fp)               :: minRoughnessLengthScale
    logical                :: minRoughnessLengthScaleAsSigma
    real(fp), dimension(3) :: initialSmoothing
    real(fp)               :: isotropicThreshold
    integer                :: boundKernelSizeFormat
    logical                :: boundKernels       = .true.
    logical                :: useGlobalSmoothing = .false.
    logical                :: isotropic          = .false. 
    real(fp)               :: initialSmoothingFactor
    real(fp), dimension(3) :: initialSmoothingArray

    ! Eventually could be used for calculating smoothing 
    ! at subsequent reconstructions
    logical :: firstRun  = .true. 
    real(fp), dimension(3) :: averageKernelSmoothing = fZERO

    ! Protocol for selection of initial smoothing
    integer :: initialSmoothingSelection
    ! Min roughness format
    integer :: minRoughnessFormat

    ! Distribution statistics
    real(fp), dimension(3) :: meanCoords
    real(fp), dimension(3) :: stdCoords
    real(fp)               :: stdSigmaScale
    real(fp)               :: hSigmaScale

    ! Limit max kernel size to fit consistently inside 
    ! the reconstruction grid.
    real(fp), dimension(3) :: maxKernelSize   
    real(fp), dimension(3) :: maxKernelSDSize
    integer                :: maxSizeDimId
    ! Limit min kernel size to at least have 2 cells 
    ! of positive shape.
    real(fp), dimension(3) :: minKernelSize   
    real(fp), dimension(3) :: minKernelSDSize
    integer                :: minSizeDimId
    ! Limit the relative growth of the support kernel
    real(fp)               :: maxSigmaGrowth

    ! Report to outUnit
    logical            :: reportToOutUnit = .false.
    integer            :: outFileUnit
    character(len=200) :: outFileName

    ! Constants defined after initialization of module dimensions
    ! Move out ?
    real(fp) :: supportDimensionConstant
    real(fp) :: alphaDimensionConstant
    real(fp) :: betaDimensionConstant
    
    ! Bins to compute
    integer, dimension(:,:), pointer    :: computeBinIds
    integer                             :: nComputeBins = 0
    character( len=300 )                :: outputFileName 
    real(fp), dimension(:,:,:), pointer :: histogramCounts  => null()
    real(fp), dimension(:,:,:), pointer :: histogramWCounts => null()

    ! Bin vector coordinates
    real(fp), dimension(:), allocatable :: coordinatesX
    real(fp), dimension(:), allocatable :: coordinatesY
    real(fp), dimension(:), allocatable :: coordinatesZ


    ! Interfaces
    procedure( SetKernelInterface )  , pass, pointer :: SetKernel      => null()
    procedure( SetKernelInterface )  , pass, pointer :: SetKernelSigma => null()
    procedure( SetKernelInterface )  , pass, pointer :: SetKernelSD1D  => null()
    procedure( SetKernelInterface2D ), pass, pointer :: SetKernelSD2D  => null()
    procedure( SetKernelInterface3D ), pass, pointer :: SetKernelSD3D  => null()
    procedure( SetKernelSDInterface ), pass, pointer :: SetKernelSD    => null()
    procedure( ComputeNetRoughness ) , pass, pointer :: ComputeNetRoughnessEstimate      => null()
    procedure( ComputeIndexes )      , pass, pointer :: ComputeKernelDatabaseIndexes     => null()
    procedure( ComputeFlatIndexes )  , pass, pointer :: ComputeKernelDatabaseFlatIndexes => null()
     
  ! GridProjectedKDEType contains
  contains
  
    ! Procedures
    procedure :: Initialize                      => prInitialize 
    procedure :: Reset                           => prReset 
    procedure :: UpdateBinSize                   => prUpdateBinSize
    procedure :: InitializeModuleDimensions      => prInitializeModuleDimensions
    procedure :: InitializeModuleConstants       => prInitializeModuleConstants
    procedure :: InitializeNetRoughnessFunction  => prInitializeNetRoughnessFunction
    procedure :: InitializeKernelDatabaseFlat    => prInitializeKernelDatabaseFlat
    procedure :: DropKernelDatabase              => prDropKernelDatabase
    procedure :: ComputeDensity                  => prComputeDensity
    procedure :: ComputeDensityOptimization      => prComputeDensityOptimization
    procedure :: ComputeCurvatureKernelBandwidth => prComputeCurvatureBandwidth 
    procedure :: ComputeOptimalSmoothingAndShape => prComputeOptimalSmoothingAndShape
    procedure :: GenerateVectorCoordinates       => prGenerateVectorCoordinates
    procedure :: ExportDensity                   => prExportDensity
    procedure :: ExportDensityUnit               => prExportDensityUnit
    procedure :: ExportDensityBinary             => prExportDensityBinary
    procedure :: ExportDensityUnitBinary         => prExportDensityUnitBinary
  
  end type GridProjectedKDEType
  

  ! Interfaces
  abstract interface
  
    ! ComputeIndexes
    function ComputeIndexes( this, smoothing ) result(indexes)
      import GridProjectedKDEType
      import fp
      implicit none
      class( GridProjectedKDEType )      :: this
      real(fp), dimension(3), intent(in) :: smoothing
      integer, dimension(3)              :: indexes 
      integer :: nd 
    end function ComputeIndexes
  
  
    ! ComputeFlatIndexes
    subroutine ComputeFlatIndexes( this, smoothing, flatDBIndexes, transposeKernel )
      import GridProjectedKDEType
      import fp
      implicit none
      class( GridProjectedKDEType )        :: this
      real(fp), dimension(3), intent(in)   :: smoothing
      integer, dimension(2), intent(inout) :: flatDBIndexes
      logical, intent(inout)               :: transposeKernel
      integer, dimension(3) :: indexes 
      integer :: nd
    end subroutine ComputeFlatIndexes
  
  
    ! NetRoughness
    subroutine ComputeNetRoughness( this, activeGridCells, curvatureBandwidth, &
                         roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                   netRoughnessArray, kernelSigmaSupportScale, &
                                               kernelSDX, kernelSDY, kernelSDZ ) 
      import GridProjectedKDEType
      import GridCellType
      import KernelSecondDerivativeXType
      import KernelSecondDerivativeYType
      import KernelSecondDerivativeZType
      import fp
      implicit none 
      class( GridProjectedKDEType ), target :: this
      type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
      real(fp), dimension(:,:), intent(in)                   :: curvatureBandwidth
      real(fp), dimension(:), intent(inout), target          :: roughnessXXArray
      real(fp), dimension(:), intent(inout), target          :: roughnessYYArray
      real(fp), dimension(:), intent(inout), target          :: roughnessZZArray
      real(fp), dimension(:), intent(inout)                  :: netRoughnessArray
      real(fp), dimension(:), intent(in)                     :: kernelSigmaSupportScale
      type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
      type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
      type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
    end subroutine ComputeNetRoughness
  

    ! SetKernelInterface
    subroutine SetKernelInterface( this, gridCell, kernel, smoothing )
      import GridProjectedKDEType
      import GridCellType
      import KernelType
      import fp
      implicit none 
      class( GridProjectedKDEType ), target      :: this
      type( GridCellType ), intent(inout)        :: gridCell
      class( KernelType ), target, intent(inout) :: kernel
      real(fp), dimension(3), intent(in)         :: smoothing
    end subroutine SetKernelInterface


    ! SetKernelInterface2D
    subroutine SetKernelInterface2D( this, gridCell, kernel1, kernel2, smoothing )
      import GridProjectedKDEType
      import GridCellType
      import KernelType
      import fp
      implicit none 
      class( GridProjectedKDEType ), target      :: this
      type( GridCellType ), intent(inout)        :: gridCell
      class( KernelType ), target, intent(inout) :: kernel1
      class( KernelType ), target, intent(inout) :: kernel2
      real(fp), dimension(3), intent(in)         :: smoothing
    end subroutine SetKernelInterface2D
  

    ! SetKernelInterface3D
    subroutine SetKernelInterface3D( this, gridCell, kernel1, kernel2, kernel3, smoothing )
      import GridProjectedKDEType
      import GridCellType
      import KernelType
      import fp
      implicit none 
      class( GridProjectedKDEType ), target      :: this
      type( GridCellType ), intent(inout)        :: gridCell
      class( KernelType ), target, intent(inout) :: kernel1
      class( KernelType ), target, intent(inout) :: kernel2
      class( KernelType ), target, intent(inout) :: kernel3
      real(fp), dimension(3), intent(in)         :: smoothing
    end subroutine SetKernelInterface3D


    ! SetKernelSDInterface
    subroutine SetKernelSDInterface(& 
      this, gridCell, kernel, smoothing, kernelDatabase, dimId )
      import GridProjectedKDEType
      import GridCellType
      import KernelType
      import fp
      implicit none 
      class( GridProjectedKDEType ), target             :: this
      type( GridCellType ),               intent(inout) :: gridCell
      class( KernelType ) , target,       intent(inout) :: kernel
      real(fp)            ,               intent(in)    :: smoothing
      class( KernelType ) , dimension(:), intent(in)    :: kernelDatabase
      integer             ,               intent(in)    :: dimId
    end subroutine SetKernelSDInterface

  end interface


! GridProjectedKDEModule contains
contains
  ! Subroutines !

  ! Some arguments candidates to be deprecated
  ! - logKernelDatabase
  subroutine prInitialize( this,& 
                domainSize, binSize, domainOrigin, &
                adaptGridToCoords, borderFraction, & 
            slicedReconstruction, slicedDimension, &
                        initialSmoothingSelection, & 
         initialSmoothing, initialSmoothingFactor, &
         nOptimizationLoops, databaseOptimization, &
    minHOverDelta, maxHOverDelta, deltaHOverDelta, &
                                logKernelDatabase, &
                          interpretAdvancedParams, &
                 minRoughnessFormat, minRoughness, & 
    minRelativeRoughness, minRoughnessLengthScale, &
                            effectiveWeightFormat, & 
                            boundKernelSizeFormat, & 
                               isotropicThreshold, & 
                                   maxSigmaGrowth, & 
                                      outFileName  )
    !---------------------------------------------------------------------------
    ! Initialize the module, assign default parameters,
    ! configures the reconstruction grid, module dimensions and others.
    !---------------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------------
    use PrecisionModule, only : fp 
    implicit none
    ! input
    class( GridProjectedKDEType ) :: this
    ! Reconstruction grid parameters
    real(fp), dimension(3), intent(in), optional :: domainSize
    real(fp), dimension(3), intent(in), optional :: binSize
    real(fp), dimension(3), intent(in), optional :: domainOrigin
    logical               , intent(in), optional :: adaptGridToCoords
    real(fp)              , intent(in), optional :: borderFraction
    ! Sliced reconstruction
    logical               , intent(in), optional :: slicedReconstruction
    integer               , intent(in), optional :: slicedDimension
    ! Initial smoothing
    real(fp), dimension(3), intent(in), optional :: initialSmoothing
    real(fp)              , intent(in), optional :: initialSmoothingFactor
    integer               , intent(in), optional :: initialSmoothingSelection
    ! Number of optimization loops
    integer               , intent(in), optional :: nOptimizationLoops 
    ! Kernel database parameters
    logical               , intent(in), optional :: databaseOptimization
    real(fp)              , intent(in), optional :: minHOverDelta
    real(fp)              , intent(in), optional :: maxHOverDelta
    real(fp)              , intent(in), optional :: deltaHOverDelta
    logical               , intent(in), optional :: logKernelDatabase    ! Deprecate ? 
    ! Advanced parameters
    logical , intent(in), optional :: interpretAdvancedParams
    integer , intent(in), optional :: minRoughnessFormat
    real(fp), intent(in), optional :: minRoughness
    real(fp), intent(in), optional :: minRelativeRoughness
    real(fp), intent(in), optional :: minRoughnessLengthScale
    integer , intent(in), optional :: effectiveWeightFormat
    integer , intent(in), optional :: boundKernelSizeFormat
    real(fp), intent(in), optional :: isotropicThreshold
    real(fp), intent(in), optional :: maxSigmaGrowth
    ! General use, indexes
    integer :: nd
    ! The analog to a listUnit, reports
    character(len=200), intent(in), optional :: outFileName
    ! local
    integer :: isThisFileOpen
    logical :: advancedOptions
    character(len=30) :: outfmt
    !---------------------------------------------------------------------------

    ! Enable reporting to outUnit if given 
    if( present( outFileName ) ) then
      isThisFileOpen = -1
      inquire( file=outFileName, number=isThisFileOpen )
      if ( isThisFileOpen .gt. 0 ) then 
        this%reportToOutUnit = .true.
        this%outFileUnit = isThisFileOpen
        this%outFileName = outFileName
        write( this%outFileUnit, * )
        write( this%outFileUnit, '(A)' ) '-----------------------'
        write( this%outFileUnit, '(A)' ) ' GPKDE is initializing '
      end if
    else if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, * )
      write( this%outFileUnit, '(A)' ) '-----------------------'
      write( this%outFileUnit, '(A)' ) ' GPKDE is initializing '
    end if 

    ! Reconstruction grid parameters !

    ! adaptGridToCoords 
    if ( present(adaptGridToCoords) ) then
      this%adaptGridToCoords = adaptGridToCoords
    else
      this%adaptGridToCoords = defaultAdaptGridToCoords
    end if

    ! borderFraction
    if ( present(borderFraction) ) then
      this%borderFraction = borderFraction
    else
      this%borderFraction = defaultBorderFraction
    end if

    ! domainSize 
    if ( present( domainSize ) ) then
      ! Some validation
      ! Forgive this error if adaptGridToCoords and infer later
      if ( .not. this%adaptGridToCoords ) then 
       if ( all( domainSize.eq.fZERO ) ) then 
        write(*,*) 'Error: all the given values for domain size are zero and grid is not adapted to coordinates.'
        stop
       end if 
      end if 
      if ( any( domainSize.lt.fZERO ) ) then 
        write(*,*) 'Error: some values of domain size are negative. They should be positive.' 
        stop
      end if 
      this%domainSize = domainSize
    else
      if ( this%adaptGridToCoords ) then 
        ! Forgive as it will be infered later
        this%domainSize = fZERO
      else
        write(*,*) 'Error: if adapt grid to coords is false, a domain size is needed while initializing GPKDE.'
        stop
      end if  
    end if  

    ! binSize 
    if ( present( binSize ).and.(.not.all(this%domainSize.eq.fZERO)) ) then 

      ! Stop if all bin sizes are zero
      if ( all( binSize .lt. fZERO ) ) then 
        write(*,*) 'Error: while initializing GPKDE, all binSizes are .lt. 0. Stop.'
        stop 
      end if 
      ! Initialize reconstruction grid parameters 
      where( binSize .ne. fZERO ) 
        this%domainGridSize = int( this%domainSize/binSize + 0.5 )
      elsewhere
        this%domainGridSize = 1
      end where
      ! Stop if any the domainGridSize .lt. 1
      if ( any( this%domainGridSize .lt. 1 ) ) then 
        write(*,*) 'Error: while initializing GPKDE, some domainGridSize  .lt. 1. Stop.'
        stop 
      end if
      this%binSize = binSize
    else if ( present( binSize ) ) then 
      ! Assign and relay init to UpdateBinSize
      ! Stop if all bin sizes are zero
      if ( all( binSize .lt. fZERO ) ) then 
        write(*,*) 'Error: while initializing GPKDE, all binSizes are .lt. 0. Stop.'
        stop 
      end if 
      this%binSize = binSize
    end if 

    ! domainOrigin
    if ( present( domainOrigin ) ) then 
      this%domainOrigin = domainOrigin
    else 
      this%domainOrigin = (/0.0_fp,0.0_fp,0.0_fp/)
    end if
    
    ! Save slicing parameters
    if ( present( slicedReconstruction ) ) then 
      this%slicedReconstruction = slicedReconstruction
    else 
      this%slicedReconstruction = .false.
    end if 
    if ( this%slicedReconstruction ) then 
      if ( present( slicedDimension ) ) then 
        this%slicedDimension = slicedDimension
      else 
        this%slicedDimension = 0
      end if
      if ( (this%slicedDimension.lt.0).or.(this%slicedDimension.gt.3) ) then 
        write(*,*) 'Error: invalid value for sliced dimension. It should be a valid dimension index. Stop.'
        stop
      end if 
      ! disable sliced reconstruction
      if ( this%slicedDimension.eq.0 ) then 
        this%slicedReconstruction = .false.
      end if  
    else
      this%slicedDimension = 0 
    end if 

    ! initialSmoothingFactor
    if ( present( initialSmoothingFactor ) ) then 
      this%initialSmoothingFactor = initialSmoothingFactor
    else
      this%initialSmoothingFactor = defaultInitialSmoothingFactor
    end if 

    ! initialSmoothingArray
    this%initialSmoothingArray = fZERO
    if ( present( initialSmoothing ) ) then
      this%initialSmoothingArray = initialSmoothing
    end if 

    ! Bin size related 
    if ( present( binSize ).and.(.not.all(this%domainSize.eq.fZERO)) ) then 

      ! Depending on domainGridSize, is the number of dimensions of the GPDKE
      ! reconstruction process. If any nBins is 1, then that dimension
      ! is compressed. e.g. nBins = (10,1,20), then it is a 2D reconstruction
      ! process where dimensions 'x' and 'z' define the 2D plane. This is not
      ! necessarily the same for the computation of histograms, where determination 
      ! of a particle inside the grid is related to the binSize. If a given binSize
      ! is zero, then histogram computation does not consider this dimension.
      ! If nBins .eq. 1 and binSize .gt. 0 then dimension is considered as valid,
      ! and compared against the origin.

      ! Initialize module dimensions
      ! Will modify the number of dimensions based on slicing parameters and report.
      call this%InitializeModuleDimensions( nDim, dimensionMask ) 

      ! Initialize module constants, uses nDim
      call this%InitializeModuleConstants()

      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(2X,A)' ) 'Initializing Histogram'
      end if

      ! Initialize histogram !
      if ( this%adaptGridToCoords ) then
        ! Skip histogram grid allocation in order 
        ! to adapt to the given particle coordinates 
        call this%histogram%Initialize(     &
         this%domainGridSize, this%binSize, &
            domainOrigin=this%domainOrigin, &
                   adaptGridToCoords=.true. )
        if ( this%reportToOutUnit ) then 
          write( this%outFileUnit, '(3X,A)' ) 'Histogram grid will not follow domain limits, will adapt to data points.'
        end if
      else
        ! Allocate grid according to nBins
        call this%histogram%Initialize(     &
         this%domainGridSize, this%binSize, &
             domainOrigin=this%domainOrigin )
        ! nBins as domainGridSize
        this%nBins = this%domainGridSize
        this%deltaBinsOrigin = 0
        ! Allocate matrix for density 
        if ( allocated( densityGrid ) ) deallocate( densityGrid )
        allocate( densityGrid(this%nBins(1), this%nBins(2), this%nBins(3)) )
        if ( this%reportToOutUnit ) then 
          write( this%outFileUnit, '(3X,A)' ) 'Histogram grid will follow domain grid size.'
        end if
      end if
      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(3X,A)' ) 'Histogram determines dimensions to be analyzed based on bin sizes.'
        write( this%outFileUnit, '(3X,A,I1,A)')&
                'Will compute Histogram considering ', this%histogram%nDim, ' dimensions.'
      end if  

      ! Process further arguments !

      ! initialSmoothing
      if ( present( initialSmoothingSelection ) ) then 
        this%initialSmoothingSelection = initialSmoothingSelection
      else
        this%initialSmoothingSelection = defaultInitialSmoothingSelection
      end if 
      this%initialSmoothing(:) = fZERO
      select case(this%initialSmoothingSelection) 
      case(0)
        ! Choose from global estimate of Silverman (1986)
        continue
      case(1)
        if ( present( initialSmoothingFactor ) ) then 
          this%initialSmoothing = initialSmoothingFactor*this%histogram%binDistance
        else
          this%initialSmoothing = defaultInitialSmoothingFactor*this%histogram%binDistance
        end if 
      case(2)
        if ( present( initialSmoothing ) ) then
          this%initialSmoothing = initialSmoothing
        else
          this%initialSmoothing = defaultInitialSmoothingFactor*this%histogram%binDistance
        end if 
      case default
        write(*,*) 'Error: Initial smoothing selection method not implemented. Stop.'
        stop
      end select
      do nd=1,3
        if ( dimensionMask(nd) .eq. 0 ) then 
          this%initialSmoothing(nd) = fZERO
        end if 
      end do

    end if 

    ! nOptimizationLoops
    if ( present( nOptimizationLoops ) ) then 
      this%nOptimizationLoops = nOptimizationLoops
    else 
      this%nOptimizationLoops = defaultNOptLoops
    end if

    ! Kernel database 
    if ( present( databaseOptimization ) ) then 
      this%databaseOptimization = databaseOptimization
    else 
      this%databaseOptimization = defaultDatabaseOptimization
    end if

    ! Bound kernel size format 
    if ( present( boundKernelSizeFormat ) ) then 
      this%boundKernelSizeFormat = boundKernelSizeFormat 
    else
      this%boundKernelSizeFormat = defaultBoundKernelSizeFormat
    end if 

    ! If kernel database or bound kernels by database params
    if (&
      (this%databaseOptimization).or.   & 
      (this%boundKernelSizeFormat.eq.1) ) then 

     ! Process database discretization parameters 
     if ( present( minHOverDelta ) ) then
       if ( minHOverDelta.gt.fZERO )then 
         this%minHOverDelta = minHOverDelta
       else
         write(*,*) 'Error: Invalid value for minHOverDelta: should be greater than zero. Stop.'
         stop
       end if
     else 
       this%minHOverDelta = defaultMinHOverDelta
     end if

     if ( present( maxHOverDelta ) ) then
      if ( maxHOverDelta.gt.fZERO ) then 
       if ( maxHOverDelta.le.this%minHOverDelta(1) ) then ! minhoverdelta is an array in memory
        if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, *) 'Error: Invalid value for maxHOverDelta: should be greater than minHOverDelta.'
        write( this%outFileUnit, *) '  Value of minHOverDelta: ', this%minHOverDelta(1)
        write( this%outFileUnit, *) '  Value of maxHOverDelta: ', maxHOverDelta
        end if  
        write(*,*) 'Error: Invalid value for maxHOverDelta: should be greater than minHOverDelta. Stop.'
        stop
       end if
       this%maxHOverDelta = maxHOverDelta
      else
       write(*,*) 'Error: Invalid value for maxHOverDelta: should be greater than zero. Stop.'
       stop
      end if
     else 
      this%maxHOverDelta = defaultMaxHOverDelta
      if ( this%maxHOverDelta(1).le.this%minHOverDelta(1) ) then 
       if ( this%reportToOutUnit ) then 
       write( this%outFileUnit, *) 'Error: Invalid value for maxHOverDelta: should be greater than minHOverDelta.'
       write( this%outFileUnit, *) '  Value of minHOverDelta: ', this%minHOverDelta(1) 
       write( this%outFileUnit, *) '  Value of maxHOverDelta: ', this%maxHOverDelta(1)
       end if  
       write(*,*) 'Error: Invalid value for maxHOverDelta: should be greater than minHOverDelta. Stop.'
       stop
      end if
     end if
     
     if ( this%databaseOptimization ) then 
      if ( present( deltaHOverDelta ) ) then 
       if ( deltaHOverDelta.gt.fZERO ) then 
        if ( deltaHOverDelta.ge.this%maxHOverDelta(1) ) then ! maxhoverdelta is an array in memory
         if ( this%reportToOutUnit ) then 
         write( this%outFileUnit, *) 'Error: Invalid value for deltaHOverDelta: should be less than maxHOverDelta.'
         write( this%outFileUnit, *) '  Value of maxHOverDelta  : ', this%maxHOverDelta(1)
         write( this%outFileUnit, *) '  Value of deltaHOverDelta: ', deltaHOverDelta
         end if  
         write(*,*) 'Error: Invalid value for deltaHOverDelta: should be less than maxHOverDelta. Stop.'
         stop
        end if
        this%deltaHOverDelta = deltaHOverDelta
       else
        write(*,*) 'Error: Invalid value for deltaHOverDelta: should be greater than zero. Stop.'
        stop
       end if
      else 
        this%deltaHOverDelta = defaultDeltaHOverDelta
        if ( this%deltaHOverDelta(1).ge.this%maxHOverDelta(1) ) then 
         if ( this%reportToOutUnit ) then 
         write( this%outFileUnit, *) 'Error: Invalid value for deltaHOverDelta: should be less than maxHOverDelta.'
         write( this%outFileUnit, *) '  Value of maxHOverDelta  : ', this%maxHOverDelta(1)
         write( this%outFileUnit, *) '  Value of deltaHOverDelta: ', this%deltaHOverDelta(1)
         end if  
         write(*,*) 'Error: Invalid value for deltaHOverDelta: should be less than maxHOverDelta. Stop.'
         stop
        end if
      end if
     end if 

    else
     ! initialize with defaults
     this%minHOverDelta   = defaultMinHOverDelta
     this%maxHOverDelta   = defaultMaxHOverDelta
     this%deltaHOverDelta = defaultDeltaHOverDelta
    end if

    if ( present( logKernelDatabase ) ) then ! Deprecate ? 
      this%logKernelDatabase = logKernelDatabase
    else 
      this%logKernelDatabase = defaultLogKernelDatabase
    end if

    ! Effective weight format 
    ! Effective weight format is defined as zero by default at histogram  
    if ( present(effectiveWeightFormat) ) then 
      this%histogram%effectiveWeightFormat = effectiveWeightFormat   
    else
      this%histogram%effectiveWeightFormat = defaultEffectiveWeightFormat   
    end if 

    ! Bin size related, and domain size 
    if ( present( binSize ).and.(.not.all(this%domainSize.eq.fZERO)) ) then 

      ! Determine kernel bounding  
      select case(this%boundKernelSizeFormat)
      ! 1: Bounding values given by user
      case(1)
       ! Assign max kernel sizes based on provided values of maxHOverDelta
       this%maxKernelSize(:) = fZERO
       do nd=1,3
        if ( this%dimensionMask(nd).eq.0 ) cycle
        this%maxKernelSize(nd) = this%binSize(nd)*maxHOverDelta
       end do
       ! As the sigma kernel is isotropic, maxSizeDimId is given by the more restrictive dimension. 
       this%maxSizeDimId = minloc( this%maxKernelSize, dim=1, mask=(this%maxKernelSize.gt.fZERO) )
       this%maxKernelSDSize(:) = fZERO
       do nd=1,3
        if ( this%dimensionMask(nd).eq.0 ) cycle
        this%maxKernelSDSize(nd) = this%binSize(nd)*maxHOverDelta
       end do
       ! Assign min kernel sizes based on provided values of minHOverDelta
       this%minKernelSize(:) = fZERO
       do nd=1,3
        if ( this%dimensionMask(nd).eq.0 ) cycle
        this%minKernelSize(nd) = this%binSize(nd)*minHOverDelta
       end do
       ! As the sigma kernel is isotropic, maxSizeDimId is given by the more restrictive dimension. 
       this%minSizeDimId = maxloc( this%minKernelSize, dim=1, mask=(this%minKernelSize.gt.fZERO) )
       this%minKernelSDSize(:) = fZERO
       do nd=1,3
        if ( this%dimensionMask(nd).eq.0 ) cycle
        this%minKernelSDSize(nd) = this%binSize(nd)*minHOverDelta
       end do
      ! 2: Unbounded 
      case(2)
        this%boundKernels = .false.
      ! 0: domain constraints
      case default
       ! Assign max kernel sizes, consistent with domain dimensions
       ! kernel ranges and bin sizes.
       this%maxKernelSize(:) = fZERO
       do nd=1,3
        if ( this%dimensionMask(nd).eq.0 ) cycle
        this%maxKernelSize(nd) = &
          this%binSize(nd)*(defaultMaxSizeFactor*this%domainGridSize(nd) - 1)/real(defaultKernelRange,fp)
       end do
       ! As the sigma kernel is isotropic, the maxSizeDimId 
       ! is given by the more restrictive dimension. 
       this%maxSizeDimId = minloc( this%maxKernelSize, dim=1, mask=(this%maxKernelSize.gt.fZERO) )
       this%maxKernelSDSize(:) = fZERO
       do nd=1,3
        if ( this%dimensionMask(nd).eq.0 ) cycle
        this%maxKernelSDSize(nd) = & 
          this%binSize(nd)*(defaultMaxSizeFactor*this%domainGridSize(nd) - 1)/real(defaultKernelSDRange,fp)
       end do
       ! Assign min kernel sizes, ensuring at least 2 positive shape cells, 
       ! Positive shape is obtained as ceiling 
       this%minKernelSize(:) = fZERO
       do nd=1,3
        if ( this%dimensionMask(nd).eq.0 ) cycle
        this%minKernelSize(nd) = defaultMinSizeFactor*this%binSize(nd)/real(defaultKernelRange,fp)
       end do
       ! As the sigma kernel is isotropic, the minSizeDimId 
       ! is given by the more restrictive dimension. 
       this%minSizeDimId = maxloc( this%minKernelSize, dim=1, mask=(this%minKernelSize.gt.fZERO) )
       this%minKernelSDSize(:) = fZERO
       do nd=1,3
        if ( this%dimensionMask(nd).eq.0 ) cycle
        this%minKernelSDSize(nd) = defaultMinSizeFactor*this%binSize(nd)/real(defaultKernelSDRange,fp)
       end do
      end select

    end if 

    ! Process advanced parameters !
     
    advancedOptions = .false.
    if ( present(interpretAdvancedParams) ) then
      advancedOptions = interpretAdvancedParams
    end if 
    if ( advancedOptions ) then
      ! Min roughness format 
      if ( present( minRoughnessFormat ) ) then 
        this%minRoughnessFormat = minRoughnessFormat
      else
        this%minRoughnessFormat = defaultMinRoughnessFormat
      end if 
      ! Isotropic threshold
      if ( present(isotropicThreshold) ) then 
        this%isotropicThreshold = isotropicThreshold
      else
        this%isotropicThreshold = defaultIsotropicThreshold
      end if
      ! Max sigma growth
      if ( present(maxSigmaGrowth) ) then 
        this%maxSigmaGrowth = maxSigmaGrowth
      else
        this%maxSigmaGrowth = defaultMaxSigmaGrowth
      end if
    else
      ! Should assign eveything to default values
      this%minRoughnessFormat    = defaultMinRoughnessFormat
      this%isotropicThreshold    = defaultIsotropicThreshold
      this%maxSigmaGrowth        = defaultMaxSigmaGrowth
    end if 

    ! Interpret roughness parameters according to format
    select case(this%minRoughnessFormat)
    ! 0: Gaussian: do nothing
    ! 1: Requires relative roughness and length scale
    case(1)
      if ( present(minRelativeRoughness) ) then 
        this%minRelativeRoughness = minRelativeRoughness
      else
        this%minRelativeRoughness = defaultMinRelativeRoughness
      end if 
      if ( present(minRoughnessLengthScale) ) then 
        this%minRoughnessLengthScale = minRoughnessLengthScale
        this%minRoughnessLengthScaleAsSigma = .false.
      else
        this%minRoughnessLengthScaleAsSigma = .true.
      end if 
    ! 2: as minRoughness
    case(2)
      if ( present(minRoughness) ) then 
        this%minLimitRoughness = minRoughness
      else
        this%minLimitRoughness = defaultMinLimitRoughness
      end if 
    ! 3: Do nothing
    end select

    ! Need more reports for roughnesses and eventually min/max kernel sizes

    ! Related to bin and domain size 
    if ( present( binSize ).and.(.not.all(this%domainSize.eq.fZERO)) ) then 

      ! Logging
      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(3X,A)') 'Grid parameters'
        write( this%outFileUnit, '(3X,A)') '---------------'
        outfmt = '(3X,A,3(1X,es18.9e3))'
        write( this%outFileUnit, outfmt) '- binSize            :', this%binSize
        write( this%outFileUnit, outfmt) '- domainSize         :', this%domainSize
        write( this%outFileUnit, outfmt) '- domainOrigin       :', this%domainOrigin
        outfmt = '(3X,A,3(1X,I9))'
        write( this%outFileUnit, outfmt) '- domainGridSize     :', this%domainGridSize
        write( this%outFileUnit, '(3X,A)') '---------------'
        write( this%outFileUnit, '(3X,A)')      'Dimensionality for reconstruction is determined from domain grid size.'
        write( this%outFileUnit, '(3X,A,I2,A)') 'Will perform reconstruction in ', nDim, ' dimensions.'
        if ( this%initialSmoothingSelection.ge.1 ) then 
        outfmt = '(3X,A,3(1X,es18.9e3))'
        write( this%outFileUnit, outfmt) '- initialSmoothing   :', this%initialSmoothing
        end if 
      end if  

      ! Initialize kernel database
      if ( this%databaseOptimization ) then
        call this%InitializeKernelDatabaseFlat( this%minHOverDelta(1), &
                                                this%maxHOverDelta(1), &
                                              this%deltaHOverDelta(1), &
                                                this%logKernelDatabase  )
        ! Pointers for SetKernel
        this%SetKernel => prSetKernelFromDatabase
        this%SetKernelSigma => prSetKernelSigmaFromDatabase
      else
        ! Pointers for SetKernel
        this%SetKernel => prSetKernelBrute
        this%SetKernelSigma => prSetKernelSigmaBrute
      end if 

      ! Initialize net roughness function
      call this%InitializeNetRoughnessFunction( nDim )

      ! Report intialization
      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(A)' ) ' GPKDE is initialized  '
        write( this%outFileUnit, '(A)' ) '-----------------------'
        write( this%outFileUnit,  *    )
        flush( this%outFileUnit ) 
      end if
    
    else

      ! Report intialization
      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(A)' ) ' GPKDE is initialized without a predefined grid, defined later. '
        write( this%outFileUnit, '(A)' ) '----------------------------------------------------------------'
        write( this%outFileUnit,  *    )
        flush( this%outFileUnit ) 
      end if

    end if 

    ! Done

  end subroutine prInitialize


  subroutine prReset( this )
    !------------------------------------------------------------------------------
    ! It could be improved, review
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ) :: this
    !------------------------------------------------------------------------------

    call this%histogram%Reset()
    dimensionMask = (/1,1,1/)
    nDim          = 3
    fNDim         = 3.0_fp
    this%dimensionMask = dimensionMask
    this%idDim1   = 0
    this%idDim2   = 0
    this%binSize  = fZERO
    this%domainSize = fZERO
    this%domainOrigin = fZERO
    this%domainGridSize = 0
    this%deltaBinsOrigin = 0
    this%nBins = 0
    this%adaptGridToCoords = .false.
    this%borderFraction = defaultBorderFraction
    this%slicedReconstruction = .false.
    this%slicedDimension = 0


    if ( allocated( this%dimensions ) ) deallocate( this%dimensions  ) 

    if ( allocated(  densityGrid             )) deallocate(  densityGrid             )
    if ( allocated(  kernelSmoothing         )) deallocate(  kernelSmoothing         )
    if ( allocated(  kernelSmoothingScale    )) deallocate(  kernelSmoothingScale    )
    if ( allocated(  kernelSmoothingShape    )) deallocate(  kernelSmoothingShape    )
    if ( allocated(  kernelSigmaSupportScale )) deallocate(  kernelSigmaSupportScale )
    if ( allocated(  curvatureBandwidth      )) deallocate(  curvatureBandwidth      )
    if ( allocated(  densityEstimateArray    )) deallocate(  densityEstimateArray    )
    if ( allocated(  nEstimateArray          )) deallocate(  nEstimateArray          )
    if ( allocated(  roughnessXXArray        )) deallocate(  roughnessXXArray        )
    if ( allocated(  roughnessYYArray        )) deallocate(  roughnessYYArray        )
    if ( allocated(  roughnessZZArray        )) deallocate(  roughnessZZArray        )
    if ( allocated(  netRoughnessArray       )) deallocate(  netRoughnessArray       )
    if ( allocated(  activeGridCellsMod      )) deallocate(  activeGridCellsMod      )

    this%kernelSDDatabase1   => null()
    this%kernelSDDatabase2   => null()
    this%roughness11Array    => null()
    this%roughness22Array    => null()
    if ( allocated( this%kernelDatabase     ) )deallocate( this%kernelDatabase     )
    if ( allocated( this%kernelDatabaseFlat ) )deallocate( this%kernelDatabaseFlat )
    if ( allocated( this%kernelSDXDatabase  ) )deallocate( this%kernelSDXDatabase  ) 
    if ( allocated( this%kernelSDYDatabase  ) )deallocate( this%kernelSDYDatabase  )
    if ( allocated( this%kernelSDZDatabase  ) )deallocate( this%kernelSDZDatabase  ) 

    this%densityEstimateGrid  => null()
    this%histogramDensity     => null()   
    this%computeBinIds        => null()
  
    if ( allocated( this%coordinatesX ) ) deallocate( this%coordinatesX )
    if ( allocated( this%coordinatesY ) ) deallocate( this%coordinatesY )
    if ( allocated( this%coordinatesZ ) ) deallocate( this%coordinatesZ )


    this%nComputeBins         = 0
    this%outputFileName       = ""
    this%deltaHOverDelta      = defaultDeltaHOverDelta
    this%minHOverDelta        = defaultMinHOverDelta
    this%maxHOverDelta        = defaultMaxHOverDelta
    this%databaseOptimization = .false.

    this%nOptimizationLoops             = defaultNOptLoops
    this%densityRelativeConvergence     = defaultRelativeErrorConvergence
    this%minLimitRoughness              = fZERO
    this%minRelativeRoughness           = defaultMinRelativeRoughness
    this%minRoughnessLengthScale        = fZERO
    this%minRoughnessLengthScaleAsSigma = .true.
    this%initialSmoothing               = fZERO
    this%isotropicThreshold             = defaultIsotropicThreshold
    this%boundKernelSizeFormat          = defaultBoundKernelSizeFormat
    this%boundKernels                   = .true.
    this%useGlobalSmoothing             = .false.
    this%isotropic                      = .false. 
    this%initialSmoothingFactor         = defaultInitialSmoothingFactor
    this%initialSmoothingArray          = fZERO

    this%firstRun = .true.
    this%averageKernelSmoothing = fZERO
    this%initialSmoothingSelection = defaultInitialSmoothingSelection
    this%minRoughnessFormat = defaultMinRoughnessFormat
    this%meanCoords = fZERO 
    this%stdCoords  = fZERO
    this%stdSigmaScale = fZERO
    this%hSigmaScale   = fZERO
        
    this%maxKernelSize   = fZERO 
    this%maxKernelSDSize = fZERO
    this%maxSizeDimId    = 0
    this%minKernelSize   = fZERO
    this%minKernelSDSize = fZERO
    this%minSizeDimId    = 0
    this%maxSigmaGrowth  = defaultMaxSigmaGrowth
    
    this%reportToOutUnit = .false.
    this%outFileUnit     = 0
    this%outFileName     = ""
    this%supportDimensionConstant = fZERO
    this%alphaDimensionConstant   = fZERO
    this%betaDimensionConstant    = fZERO
     

    this%SetKernel      => null()
    this%SetKernelSigma => null()
    this%SetKernelSD1D  => null()
    this%SetKernelSD2D  => null()
    this%SetKernelSD3D  => null()
    this%SetKernelSD    => null()
    this%ComputeNetRoughnessEstimate       => null()
    this%ComputeKernelDatabaseIndexes      => null()
    this%ComputeKernelDatabaseFlatIndexes  => null()

    if ( associated(this%histogramCounts)  ) this%histogramCounts  => null()
    if ( associated(this%histogramWCounts) ) this%histogramWCounts => null()


  end subroutine prReset

  
  subroutine prUpdateBinSize( this, binSize )
    !---------------------------------------------------------------------------
    ! Performs the part of initialization related to the bin size selection. 
    !---------------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------------
    use PrecisionModule, only : fp 
    implicit none
    ! input
    class( GridProjectedKDEType ) :: this
    ! Reconstruction grid parameters
    real(fp), dimension(3), intent(in) :: binSize
    ! General use, indexes
    integer :: nd
    ! local
    character(len=30) :: outfmt
    !---------------------------------------------------------------------------

    ! Validate bin size 
    if ( all( binSize .le. fZERO ) ) then 
      write(*,*) 'Error: while updating bin size at GPKDE, all values are .lt. 0. Stop.'
      stop 
    end if
    if ( any( binSize .lt. fZERO ) ) then 
      write(*,*) 'Error: while updating bin size at GPKDE, some values are negative. Stop.'
      stop 
    end if
    ! Validate domain size 
    if ( all( this%domainSize .eq. fZERO ) ) then 
      write(*,*) 'Error: while updating bin size at GPKDE. All dimensions are zero. Stop.'
      stop 
    end if

    ! Initialize reconstruction grid parameters !

    where( binSize .ne. fZERO )
      ! domainSize expected to be already defined into the object 
      this%domainGridSize = int( this%domainSize/binSize + 0.5 )
    elsewhere
      this%domainGridSize = 1
    end where

    ! Stop if any the domainGridSize .lt. 1
    if ( any( this%domainGridSize .lt. 1 ) ) then 
      write(*,*) 'Error: while initializing GPKDE, some domainGridSize  .lt. 1. Stop.'
      stop 
    end if
    this%binSize = binSize

    ! Depending on domainGridSize, is the number of dimensions of the GPDKE
    ! reconstruction process. If any nBins is 1, then that dimension
    ! is compressed. e.g. nBins = (10,1,20), then it is a 2D reconstruction
    ! process where dimensions 'x' and 'z' define the 2D plane. This is not
    ! necessarily the same for the computation of histograms, where determination 
    ! of a particle inside the grid is related to the binSize. If a given binSize
    ! is zero, then histogram computation does not consider this dimension.
    ! If nBins .eq. 1 and binSize .gt. 0 then dimension is considered as valid,
    ! and compared against the origin.

    ! Initialize module dimensions
    ! Will modify the number of dimensions based on slicing parameters and report.
    ! dimensionMask and nDim are defined at the module level
    call this%InitializeModuleDimensions( nDim, dimensionMask ) 

    ! Initialize module constants, uses nDim
    call this%InitializeModuleConstants()

    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, '(2X,A)' ) 'Initializing Histogram'
    end if

    ! Initialize histogram !
    if ( this%adaptGridToCoords ) then
      ! Skip histogram grid allocation in order 
      ! to adapt to the given particle coordinates 
      call this%histogram%Initialize(     &
       this%domainGridSize, this%binSize, &
          domainOrigin=this%domainOrigin, &
                 adaptGridToCoords=.true. )
      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(3X,A)' ) 'Histogram grid will not follow domain limits, will adapt to data points.'
      end if
    else
      ! Allocate grid according to nBins
      call this%histogram%Initialize(     &
       this%domainGridSize, this%binSize, &
           domainOrigin=this%domainOrigin )
      ! nBins as domainGridSize
      this%nBins = this%domainGridSize
      this%deltaBinsOrigin = 0
      ! Allocate matrix for density 
      if ( allocated( densityGrid ) ) deallocate( densityGrid )
      allocate( densityGrid(this%nBins(1), this%nBins(2), this%nBins(3)) )
      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(3X,A)' ) 'Histogram grid will follow domain grid size.'
      end if
    end if
    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, '(3X,A)' ) 'Histogram determines dimensions to be analyzed based on bin sizes.'
      write( this%outFileUnit, '(3X,A,I1,A)')&
              'Will compute Histogram considering ', this%histogram%nDim, ' dimensions.'
    end if  

    ! Further processes !

    ! initialSmoothing
    select case(this%initialSmoothingSelection) 
    case(0)
      ! Choose from global estimate of Silverman (1986)
      continue
    case(1)
      this%initialSmoothing = this%initialSmoothingFactor*this%histogram%binDistance
    case(2)
      if ( all( this%initialSmoothingArray .le. fZERO ) ) then 
        this%initialSmoothing = defaultInitialSmoothingFactor*this%histogram%binDistance
      else
        this%initialSmoothing = this%initialSmoothingArray
      end if 
    case default
      write(*,*) 'Error: Initial smoothing selection method not implemented. Stop.'
      stop
    end select
    do nd=1,3
      if ( dimensionMask(nd) .eq. 0 ) then 
        this%initialSmoothing(nd) = fZERO
      end if 
    end do

    ! Determine kernel bounding  
    select case(this%boundKernelSizeFormat)
    ! 1: Bounding values given by user
    case(1)
     ! Assign max kernel sizes based on provided values of maxHOverDelta
     this%maxKernelSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%maxKernelSize(nd) = this%binSize(nd)*this%maxHOverDelta(1)
     end do
     ! As the sigma kernel is isotropic, maxSizeDimId is given by the more restrictive dimension. 
     this%maxSizeDimId = minloc( this%maxKernelSize, dim=1, mask=(this%maxKernelSize.gt.fZERO) )
     this%maxKernelSDSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%maxKernelSDSize(nd) = this%binSize(nd)*this%maxHOverDelta(1)
     end do
     ! Assign min kernel sizes based on provided values of minHOverDelta
     this%minKernelSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%minKernelSize(nd) = this%binSize(nd)*this%minHOverDelta(1)
     end do
     ! As the sigma kernel is isotropic, maxSizeDimId is given by the more restrictive dimension. 
     this%minSizeDimId = maxloc( this%minKernelSize, dim=1, mask=(this%minKernelSize.gt.fZERO) )
     this%minKernelSDSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%minKernelSDSize(nd) = this%binSize(nd)*this%minHOverDelta(1)
     end do
    ! 2: Unbounded 
    case(2)
      this%boundKernels = .false.
    ! 0: domain constraints
    case default
     ! Assign max kernel sizes, consistent with domain dimensions
     ! kernel ranges and bin sizes.
     this%maxKernelSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%maxKernelSize(nd) = &
        this%binSize(nd)*(defaultMaxSizeFactor*this%domainGridSize(nd) - 1)/real(defaultKernelRange,fp)
     end do
     ! As the sigma kernel is isotropic, the maxSizeDimId 
     ! is given by the more restrictive dimension. 
     this%maxSizeDimId = minloc( this%maxKernelSize, dim=1, mask=(this%maxKernelSize.gt.fZERO) )
     this%maxKernelSDSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%maxKernelSDSize(nd) = & 
        this%binSize(nd)*(defaultMaxSizeFactor*this%domainGridSize(nd) - 1)/real(defaultKernelSDRange,fp)
     end do
     ! Assign min kernel sizes, ensuring at least 2 positive shape cells, 
     ! Positive shape is obtained as ceiling 
     this%minKernelSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%minKernelSize(nd) = defaultMinSizeFactor*this%binSize(nd)/real(defaultKernelRange,fp)
     end do
     ! As the sigma kernel is isotropic, the minSizeDimId 
     ! is given by the more restrictive dimension. 
     this%minSizeDimId = maxloc( this%minKernelSize, dim=1, mask=(this%minKernelSize.gt.fZERO) )
     this%minKernelSDSize(:) = fZERO
     do nd=1,3
      if ( this%dimensionMask(nd).eq.0 ) cycle
      this%minKernelSDSize(nd) = defaultMinSizeFactor*this%binSize(nd)/real(defaultKernelSDRange,fp)
     end do
    end select


    ! Logging
    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, '(3X,A)') 'Grid parameters'
      write( this%outFileUnit, '(3X,A)') '---------------'
      outfmt = '(3X,A,3(1X,es18.9e3))'
      write( this%outFileUnit, outfmt) '- binSize            :', this%binSize
      write( this%outFileUnit, outfmt) '- domainSize         :', this%domainSize
      write( this%outFileUnit, outfmt) '- domainOrigin       :', this%domainOrigin
      outfmt = '(3X,A,3(1X,I9))'
      write( this%outFileUnit, outfmt) '- domainGridSize     :', this%domainGridSize
      write( this%outFileUnit, '(3X,A)') '---------------'
      write( this%outFileUnit, '(3X,A)')      'Dimensionality for reconstruction is determined from domain grid size.'
      write( this%outFileUnit, '(3X,A,I2,A)') 'Will perform reconstruction in ', nDim, ' dimensions.'
      if ( this%initialSmoothingSelection.ge.1 ) then 
      outfmt = '(3X,A,3(1X,es18.9e3))'
      write( this%outFileUnit, outfmt) '- initialSmoothing   :', this%initialSmoothing
      end if 
    end if  

    ! Initialize kernel database
    if ( this%databaseOptimization ) then
      ! Related to nDim, which is obtained after binsize is given
      ! it should avoid multiple initialization
      if ( .not. allocated( this%kernelDatabaseFlat ) ) then 
        call this%InitializeKernelDatabaseFlat( this%minHOverDelta(1), &
                                                this%maxHOverDelta(1), &
                                              this%deltaHOverDelta(1), &
                                                this%logKernelDatabase  )
        ! Pointers for SetKernel
        this%SetKernel => prSetKernelFromDatabase
        this%SetKernelSigma => prSetKernelSigmaFromDatabase
      end if 
    else
      ! Pointers for SetKernel
      this%SetKernel => prSetKernelBrute
      this%SetKernelSigma => prSetKernelSigmaBrute
    end if 

    ! Initialize net roughness function
    call this%InitializeNetRoughnessFunction( nDim )


    !! Report intialization
    !if ( this%reportToOutUnit ) then 
    !  write( this%outFileUnit, '(A)' ) ' GPKDE is initialized  '
    !  write( this%outFileUnit, '(A)' ) '-----------------------'
    !  write( this%outFileUnit,  *    )
    !  flush( this%outFileUnit ) 
    !end if

    ! Done


  end subroutine prUpdateBinSize 



  subroutine prAllocateArrays( this, nComputeBins,      &
                              inkernelSmoothing,        &
                              inkernelSmoothingScale,   &
                              inkernelSmoothingShape,   &
                              inkernelSigmaSupportScale,&
                              incurvatureBandwidth,     &
                              indensityEstimateArray,   &
                              innEstimateArray,         &
                              inroughnessXXArray,       &
                              inroughnessYYArray,       &
                              inroughnessZZArray,       &
                              innetRoughnessArray,      &
                              activeGridCellsIn         )
    !------------------------------------------------------------------------------
    ! It could be improved, review
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ) :: this
    integer, intent(in) :: nComputeBins
    real(fp), dimension(:,:), allocatable, intent(out) :: inkernelSmoothing
    real(fp), dimension(:)  , allocatable, intent(out) :: inkernelSmoothingScale
    real(fp), dimension(:,:), allocatable, intent(out) :: inkernelSmoothingShape
    real(fp), dimension(:)  , allocatable, intent(out) :: inkernelSigmaSupportScale
    real(fp), dimension(:,:), allocatable, intent(out) :: incurvatureBandwidth
    real(fp), dimension(:)  , allocatable, intent(out) :: indensityEstimateArray 
    real(fp), dimension(:)  , allocatable, intent(out) :: innEstimateArray
    real(fp), dimension(:)  , allocatable, intent(out) :: inroughnessXXArray
    real(fp), dimension(:)  , allocatable, intent(out) :: inroughnessYYArray
    real(fp), dimension(:)  , allocatable, intent(out) :: inroughnessZZArray
    real(fp), dimension(:)  , allocatable, intent(out) :: innetRoughnessArray
    real(fp), dimension(:,:), allocatable :: lockernelSmoothing
    real(fp), dimension(:)  , allocatable :: lockernelSmoothingScale
    real(fp), dimension(:,:), allocatable :: lockernelSmoothingShape
    real(fp), dimension(:)  , allocatable :: lockernelSigmaSupportScale
    real(fp), dimension(:,:), allocatable :: loccurvatureBandwidth
    real(fp), dimension(:)  , allocatable :: locdensityEstimateArray 
    real(fp), dimension(:)  , allocatable :: locnEstimateArray
    real(fp), dimension(:)  , allocatable :: locroughnessXXArray
    real(fp), dimension(:)  , allocatable :: locroughnessYYArray
    real(fp), dimension(:)  , allocatable :: locroughnessZZArray
    real(fp), dimension(:)  , allocatable :: locnetRoughnessArray
    type( GridCellType ), dimension(:), allocatable, intent(out) :: activeGridCellsIn
    type( GridCellType ), dimension(:), allocatable :: activeGridCellsLocal
    !------------------------------------------------------------------------------

    ! Allocate arrays
    allocate(      lockernelSmoothing( 3, nComputeBins ) )
    allocate( lockernelSmoothingShape( 3, nComputeBins ) )
    allocate(   loccurvatureBandwidth( 3, nComputeBins ) )
    allocate(    lockernelSmoothingScale( nComputeBins ) )
    allocate( lockernelSigmaSupportScale( nComputeBins ) )
    allocate(    locdensityEstimateArray( nComputeBins ) )
    allocate(          locnEstimateArray( nComputeBins ) )
    allocate(       locnetRoughnessArray( nComputeBins ) )
    allocate(       activeGridCellsLocal( nComputeBins ) )

    call move_alloc(          lockernelSmoothing,        inkernelSmoothing )
    call move_alloc(     lockernelSmoothingShape,   inkernelSmoothingShape )
    call move_alloc(       loccurvatureBandwidth,     incurvatureBandwidth )
    call move_alloc(     lockernelSmoothingScale,   inkernelSmoothingScale )
    call move_alloc(  lockernelSigmaSupportScale,inkernelSigmaSupportScale )
    call move_alloc(     locdensityEstimateArray,   indensityEstimateArray )
    call move_alloc(           locnEstimateArray,         innEstimateArray )
    call move_alloc(        locnetRoughnessArray,      innetRoughnessArray )
    call move_alloc(        activeGridCellsLocal,       activeGridCellsIn  )

    ! roughness
    if (this%dimensionMask(1).eq.1) then 
      allocate(locroughnessXXArray(nComputeBins))
      call move_alloc(locroughnessXXArray,inroughnessXXArray)
    end if 
    if (this%dimensionMask(2).eq.1) then 
      allocate(locroughnessYYArray(nComputeBins))
      call move_alloc(locroughnessYYArray,inroughnessYYArray)
    end if 
    if (this%dimensionMask(3).eq.1) then 
      allocate(locroughnessZZArray(nComputeBins))
      call move_alloc(locroughnessZZArray,inroughnessZZArray)
    end if 


  end subroutine prAllocateArrays


  subroutine prInitializeModuleDimensions( this, nDim, dimensionMask )
  !-----------------------------------------------------------------
  !
  !-----------------------------------------------------------------
  ! Specifications 
  !-----------------------------------------------------------------
  class( GridProjectedKDEType ), target :: this 
  integer, intent(inout)                :: nDim
  integer, dimension(3), intent(inout)  :: dimensionMask
  integer :: n, nd, currentDim, dcount
  !-----------------------------------------------------------------

    ! Determine dimensions based on number of bins
    do n = 1,3
      if (this%domainGridSize(n) .eq. 1) dimensionMask(n) = 0 
    end do 
    nDim  = sum(dimensionMask)
    fNDim = real(nDim,fp)
    this%dimensionMask = dimensionMask
    if ( nDim .le. 0 ) then 
      write(*,*) 'Error: while initializing GPKDE dimensions. nDim .le. 0. Stop.'
      stop
    end if 

    ! Now, contrast the dimensionality information with 
    ! parameters for sliced reconstruction 
    if ( (nDim.eq.1).and.(this%slicedReconstruction) ) then 
      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(1X,A)' ) 'Number of dimensions is 1. Will ignore sliced reconstruction parameters.'
      end if 
      this%slicedReconstruction = .false.
      this%slicedDimension = 0
    end if 

    ! If remained as true, slicedDimension was a valid index.
    ! Verify that the dimension is active
    if ( this%slicedReconstruction ) then
      if ( this%dimensionMask(this%slicedDimension).eq.0 ) then 
        if ( this%reportToOutUnit ) then 
         write( this%outFileUnit, '(1X,A)' ) 'Sliced dimension is inactive for kernels, will ignore sliced reconstruction.'
        end if 
        this%slicedReconstruction = .false.
        this%slicedDimension = 0
      end if 
    end if 

    ! If remained as true, then the slicing dimension is 
    ! active. Substract one to the number of dimensions 
    ! and modify dimension mask to make it "inactive" (compress kernels).
    if ( this%slicedReconstruction ) then
      dimensionMask(this%slicedDimension) = 0
      nDim  = sum(dimensionMask)
      fNDim = real(nDim,fp)
      this%dimensionMask = dimensionMask
      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(1X,A,I2)' ) 'Reconstruction is sliced in dimension ', this%slicedDimension 
      end if
    end if 
    ! Downstream this point, the initialization of kernels, dimensions 
    ! and constants, is consistent for the purposes of sliced reconstruction.

    ! Initialize dimensions in kernel module
    call InitializeKernelDimensions(dimensionMask)

    ! Identify directions, the OLD way
    ! 1D
    if ( nDim .eq. 1 ) then 
      ! Relate x,y,z dimensions to 1 dimensions
      do nd = 1,3
        if ( this%dimensionMask( nd ) .eq. 0 ) cycle
        select case(nd) 
          case (1)
            this%idDim1 = nd
          case (2)
            this%idDim1 = nd
          case (3)
            this%idDim1 = nd
        end select   
        ! Use the first found
        exit
      end do
    end if

    ! 2D
    if ( nDim .eq. 2 ) then
      ! Relate x,y,z dimensions to 1,2 dimensions
      do nd = 1,3
        if ( this%dimensionMask( nd ) .eq. 0 ) cycle
        currentDim = sum( this%dimensionMask(1:nd) )
        select case(nd) 
          case (1)
            this%idDim1 = nd
          case (2)
            if ( currentDim .eq. 1 ) then 
              this%idDim1 = nd
            else if ( currentDim .eq. 2 ) then
              this%idDim2 = nd
            end if
          case (3)
            this%idDim2 = nd
        end select   
      end do
    end if

    ! The NEW way
    if ( allocated( this%dimensions ) ) deallocate( this%dimensions )
    allocate( this%dimensions( nDim  ) )
    dcount= 0
    do nd = 1, 3
      if ( this%dimensionMask(nd) .eq. 0 ) cycle
      dcount = dcount + 1
      this%dimensions(dcount) = nd
    end do


    ! Done
    return


  end subroutine prInitializeModuleDimensions 


  subroutine prInitializeModuleConstants( this )
    !------------------------------------------------------------------------------
    ! Constants:
    !   - supportDimensionConstant for Eq. (23) in Sole-Mari et al. (2019)
    !   - alphaDimensionConstant
    !     and betaDimensionConstant, Eq. (28) in Sole-Mari et al. (2019)
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    class( GridProjectedKDEType ) :: this 
    !------------------------------------------------------------------------------

    ! Compute constants
    this%supportDimensionConstant = ( ( fNDim + fTWO )*( fEIGHT*pi )**( 0.5*fNDim ) )**( 0.25 )

    this%alphaDimensionConstant =& 
        ( ( fONE + fTWO**(0.5*fNDim + fTWO) )/( fTHREE*fTWO**( fFOUR/( fNDim + fFOUR ) ) ) )**(&
        fONE/(fNDim + fSIX) )*( fNDim + fTWO )**( fONE/(fNDim + fFOUR) )/( ( fNDim + fFOUR )**( fONE/(fNDim + fSIX) ) )

    this%betaDimensionConstant  = fTWO/( fNDim + fFOUR)/( fNDim + fSIX ) 

    ! Appear recurrently in expressions
    oneOverNDimPlusFour     = fONE/( fNDim + fFOUR )
    minusOneOverNDimPlusSix = -fONE/( fNDim + fSIX )
    onePlusNDimQuarter      = fONE + 0.25*fNDim

    ! Done
    return

  end subroutine prInitializeModuleConstants 


  ! NET ROUGHNESS
  ! net roughness
  ! initialize
  subroutine prInitializeNetRoughnessFunction( this, nDim )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    class( GridProjectedKDEType ), target :: this 
    integer, intent(in)                   :: nDim
    !------------------------------------------------------------------------------

    ! Assign interface depending on dimensionality
    if ( nDim .eq. 1 ) then 
      this%ComputeNetRoughnessEstimate => prComputeNetRoughness1D
      if ( this%databaseOptimization ) then 
        this%SetKernelSD1D => prSetKernelSD1DFromDatabase
      else
        this%SetKernelSD1D => prSetKernelSD1DBrute
      end if 
    end if
    if ( nDim .eq. 2 ) then
      this%ComputeNetRoughnessEstimate => prComputeNetRoughness2D
      if ( this%databaseOptimization ) then 
        this%SetKernelSD2D => prSetKernelSD2DFromDatabase
      else
        this%SetKernelSD2D => prSetKernelSD2DBrute
      end if 
    end if
    if ( nDim .eq. 3 ) then
#ifdef __INTEL_COMPILER
      this%ComputeNetRoughnessEstimate => prComputeNetRoughness3DIndep
#else
      this%ComputeNetRoughnessEstimate => prComputeNetRoughness3D
#endif
      if ( this%databaseOptimization ) then 
        this%SetKernelSD => prSetKernelSDFromDatabase
      else
        this%SetKernelSD => prSetKernelSDBrute
      end if 
    end if

    ! Done
    return

  end subroutine prInitializeNetRoughnessFunction 


  ! NET ROUGHNESS
  ! net roughness
  ! 1D
  subroutine prComputeNetRoughness1D( this, activeGridCells, curvatureBandwidth, &
                           roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                     netRoughnessArray, kernelSigmaSupportScale, &
                                                 kernelSDX, kernelSDY, kernelSDZ ) 
    !-----------------------------------------------------------------------------
    ! Net roughness in 1D
    ! 
    !  - Eq. 13a in Sole-Mari et al. (2019)
    !
    !-----------------------------------------------------------------------------
    ! Specifications 
    !-----------------------------------------------------------------------------
    !input
    class( GridProjectedKDEType ), target                  :: this
    type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
    real(fp), dimension(:,:), intent(in)                   :: curvatureBandwidth
    real(fp), dimension(:), intent(in)                     :: kernelSigmaSupportScale
    type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
    type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
    type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
    ! out
    real(fp), dimension(:), intent(inout), target   :: roughnessXXArray
    real(fp), dimension(:), intent(inout), target   :: roughnessYYArray
    real(fp), dimension(:), intent(inout), target   :: roughnessZZArray
    real(fp), dimension(:), intent(inout)           :: netRoughnessArray
    ! local 
    type( GridCellType ), pointer                   :: gc => null()
    real(fp), dimension(:), pointer                 :: roughness11Array
    real(fp), dimension(:,:,:), allocatable, target :: curvature1
    type( KernelMultiGaussianType )                 :: kernelSigma
    real(fp), dimension(3)                          :: smoothingCarrier
    integer :: n
    !-----------------------------------------------------------------------------
    allocate( curvature1(  this%nBins(1), this%nBins(2), this%nBins(3) )) 
    !-----------------------------------------------------------------------------
  
    ! Initialize
    curvature1  = fZERO

    ! Assign dimension pointers
    select case( this%idDim1 ) 
      case (1)
        roughness11Array => roughnessXXArray
        !$omp parallel do schedule( dynamic, 1 ) & 
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( activeGridCells )          &
        !$omp shared( curvatureBandwidth )       &
        !$omp reduction( +:curvature1 )          &
        !$omp firstprivate( kernelSDX )          &
        !$omp private( n )                       &
        !$omp private( gc )                       
        do n = 1, this%nComputeBins
  
          ! Assign gc pointer 
          gc => activeGridCells(n)
  
          if ( ( curvatureBandwidth(this%idDim1,n) .lt. fZERO ) .or. & 
               ( curvatureBandwidth(this%idDim1,n) /= curvatureBandwidth(this%idDim1,n) ) ) cycle

          ! Set kernel
          call this%SetKernelSD1D( gc, kernelSDX, curvatureBandwidth(:,n) )

          ! Compute curvature
          curvature1( &
                  gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                  gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                  gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
              ) = curvature1( &
                  gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                  gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                  gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
              ) + this%histogramCounts(                             &
              !) + this%histogram%counts(                             &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                          gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                          gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                          gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                   )
        end do
        !$omp end parallel do

      case (2)
        roughness11Array => roughnessYYArray
        !$omp parallel do schedule( dynamic, 1 ) & 
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( activeGridCells )          &
        !$omp shared( curvatureBandwidth )       &
        !$omp reduction( +:curvature1 )          &
        !$omp firstprivate( kernelSDY )          &
        !$omp private( n )                       &
        !$omp private( gc )                       
        do n = 1, this%nComputeBins
  
          ! Assign gc pointer 
          gc => activeGridCells(n)
  
          if ( ( curvatureBandwidth(this%idDim1,n) .lt. fZERO ) .or. & 
               ( curvatureBandwidth(this%idDim1,n) /= curvatureBandwidth(this%idDim1,n) ) ) cycle

          ! Set kernel
          call this%SetKernelSD1D( gc, kernelSDY, curvatureBandwidth(:,n) )

          ! Compute curvature
          curvature1( &
                  gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                  gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                  gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
              ) = curvature1( &
                  gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                  gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                  gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
              ) + this%histogramCounts(                             &
              !) + this%histogram%counts(                             &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                          gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                          gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                          gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                   )
        end do
        !$omp end parallel do

      case (3)
        roughness11Array => roughnessZZArray
        !$omp parallel do schedule( dynamic, 1 ) & 
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( activeGridCells )          &
        !$omp shared( curvatureBandwidth )       &
        !$omp reduction( +:curvature1 )          &
        !$omp firstprivate( kernelSDZ )          &
        !$omp private( n )                       &
        !$omp private( gc )                       
        do n = 1, this%nComputeBins
  
          ! Assign gc pointer 
          gc => activeGridCells(n)
  
          if ( ( curvatureBandwidth(this%idDim1,n) .lt. fZERO ) .or. & 
               ( curvatureBandwidth(this%idDim1,n) /= curvatureBandwidth(this%idDim1,n) ) ) cycle

          ! Set kernel
          call this%SetKernelSD1D( gc, kernelSDZ, curvatureBandwidth(:,n) )

          ! Compute curvature
          curvature1( &
                  gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                  gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                  gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
              ) = curvature1( &
                  gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                  gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                  gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
              ) + this%histogramCounts(                             &
              !) + this%histogram%counts(                             &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                          gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                          gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                          gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                   )
        end do
        !$omp end parallel do
    end select
    call kernelSDX%ResetMatrix()
    call kernelSDY%ResetMatrix()
    call kernelSDZ%ResetMatrix()

    curvature1 = curvature1/this%histogram%binVolume
    ! Matrix from curvature kernels is delta**2*KernelVMatrix
    curvature1 = curvature1/( this%binSize(this%idDim1)**fTWO )

    ! Product curvature
    curvature1 = curvature1*curvature1

    ! Net roughness
    roughness11Array  = fZERO 
    netRoughnessArray = fZERO
    if ( this%useGlobalSmoothing ) then
      ! If global smoothing, then the integral is over 
      ! the whole domain
      roughness11Array(:) = sum(curvature1)
      netRoughnessArray = roughness11Array

      ! Deallocate
      deallocate( curvature1 )
      
      ! Done
      return
    end if

    ! Continue to local form !

    ! kernelSigma was already computed ? 
    call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )

    smoothingCarrier  = fZERO
    !$omp parallel do schedule( dynamic, 1 ) &
    !$omp default( none )                    &
    !$omp shared( this )                     &
    !$omp shared( activeGridCells )          &
    !$omp shared( curvature1 )               &
    !$omp shared( roughness11Array )         &
    !$omp shared( kernelSigmaSupportScale )  & 
    !$omp firstprivate( kernelSigma )        &
    !$omp firstprivate( smoothingCarrier )   &
    !$omp private( n )                       &
    !$omp private( gc )
    do n = 1, this%nComputeBins

      ! Assign pointer 
      gc => activeGridCells(n)

      smoothingCarrier(:) = kernelSigmaSupportScale(n)
      call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

      ! Compute roughness grid estimates
      roughness11Array( n ) = sum(&
          curvature1(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
    end do
    !$omp end parallel do 
    netRoughnessArray = roughness11Array

    ! Deallocate
    deallocate( curvature1 )

    ! Done
    return

  end subroutine prComputeNetRoughness1D


  ! NET ROUGHNESS
  ! net roughness
  ! 2D
  subroutine prComputeNetRoughness2D( this, activeGridCells, curvatureBandwidth, &
                           roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                     netRoughnessArray, kernelSigmaSupportScale, &
                                                 kernelSDX, kernelSDY, kernelSDZ ) 
    !-----------------------------------------------------------------------------
    ! Net roughness in 2D
    ! 
    !  - Eq. 13b in Sole-Mari et al. (2019)
    ! 
    !-----------------------------------------------------------------------------
    ! Specifications 
    !-----------------------------------------------------------------------------
    ! input
    class( GridProjectedKDEType ), target                  :: this
    type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
    real(fp), dimension(:,:), intent(in)                   :: curvatureBandwidth
    real(fp), dimension(:), intent(in)                     :: kernelSigmaSupportScale
    type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
    type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
    type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
    ! out
    real(fp), dimension(:), intent(inout), target   :: roughnessXXArray
    real(fp), dimension(:), intent(inout), target   :: roughnessYYArray
    real(fp), dimension(:), intent(inout), target   :: roughnessZZArray
    real(fp), dimension(:), intent(inout)           :: netRoughnessArray
    ! local
    type( GridCellType ), pointer                   :: gc => null()
    real(fp), dimension(:), pointer                 :: roughness11Array
    real(fp), dimension(:), pointer                 :: roughness22Array
    real(fp), dimension(:,:,:), allocatable         :: curvature1
    real(fp), dimension(:,:,:), allocatable         :: curvature2
    real(fp), dimension(:,:,:), allocatable         :: curvature12
    type( KernelMultiGaussianType )                 :: kernelSigma
    real(fp), dimension(3)                          :: smoothingCarrier
    integer :: n
    !-----------------------------------------------------------------------------
    allocate( curvature1( this%nBins(1), this%nBins(2), this%nBins(3)  )) 
    allocate( curvature2( this%nBins(1), this%nBins(2), this%nBins(3)  )) 
    allocate( curvature12( this%nBins(1), this%nBins(2), this%nBins(3) )) 
    !-----------------------------------------------------------------------------

    ! Initialize 
    curvature1  = fZERO
    curvature2  = fZERO

    ! Choose curvature computation
    if ( ( this%idDim1 .eq. 1 ) .and. ( this%idDim2 .eq. 2 ) ) then 
      ! XY
      roughness11Array => roughnessXXArray
      roughness22Array => roughnessYYArray
      !$omp parallel do schedule( dynamic, 1 ) & 
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvatureBandwidth )       &
      !$omp reduction( +:curvature1 )          &
      !$omp reduction( +:curvature2 )          &
      !$omp firstprivate( kernelSDX )          &
      !$omp firstprivate( kernelSDY )          &
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins
  
        ! Assign gc pointer 
        gc => activeGridCells(n)
  
        if ( ( any( curvatureBandwidth(:,n) .lt. fZERO ) ) .or. & 
             ( any( curvatureBandwidth(:,n) /= curvatureBandwidth(:,n) ) ) ) cycle

        ! Set kernels        
        call this%SetKernelSD2D( gc, kernelSDX, kernelSDY, curvatureBandwidth( :, n ) )

        ! Compute curvature
        curvature1( &
                gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
            ) = curvature1( &
                gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
            ) + this%histogramCounts(                                &
            !) + this%histogram%counts(                               &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                        gc%kernelSD1XMSpan(1):gc%kernelSD1XMSpan(2), &
                        gc%kernelSD1YMSpan(1):gc%kernelSD1YMSpan(2), & 
                        gc%kernelSD1ZMSpan(1):gc%kernelSD1ZMSpan(2)  & 
                   )
        ! Compute curvature
        curvature2( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) = curvature2( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) + this%histogramCounts(                                &
            !) + this%histogram%counts(                               &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD2Matrix(&
                        gc%kernelSD2XMSpan(1):gc%kernelSD2XMSpan(2), &
                        gc%kernelSD2YMSpan(1):gc%kernelSD2YMSpan(2), & 
                        gc%kernelSD2ZMSpan(1):gc%kernelSD2ZMSpan(2)  & 
                   )
      end do
      !$omp end parallel do
    else if ( ( this%idDim1 .eq. 1 ) .and. ( this%idDim2 .eq. 3 ) ) then
      ! XZ
      roughness11Array => roughnessXXArray
      roughness22Array => roughnessZZArray
      !$omp parallel do schedule( dynamic, 1 ) & 
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvatureBandwidth )       &
      !$omp reduction( +:curvature1 )          &
      !$omp reduction( +:curvature2 )          &
      !$omp firstprivate( kernelSDX )          &
      !$omp firstprivate( kernelSDZ )          &
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins
  
        ! Assign gc pointer 
        gc => activeGridCells(n)
  
        if ( ( any( curvatureBandwidth( :, n ) .lt. fZERO ) ) .or. & 
             ( any( curvatureBandwidth( :, n ) /= curvatureBandwidth( :, n ) ) ) ) cycle

        ! Set kernels        
        call this%SetKernelSD2D( gc, kernelSDX, kernelSDZ, curvatureBandwidth( :, n ) )

        ! Compute curvature
        curvature1( &
                gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
            ) = curvature1( &
                gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
            ) + this%histogramCounts(                                &
            !) + this%histogram%counts(                               &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                        gc%kernelSD1XMSpan(1):gc%kernelSD1XMSpan(2), &
                        gc%kernelSD1YMSpan(1):gc%kernelSD1YMSpan(2), & 
                        gc%kernelSD1ZMSpan(1):gc%kernelSD1ZMSpan(2)  & 
                   )
        ! Compute curvature
        curvature2( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) = curvature2( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) + this%histogramCounts(                                &
            !) + this%histogram%counts(                               &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD2Matrix(&
                        gc%kernelSD2XMSpan(1):gc%kernelSD2XMSpan(2), &
                        gc%kernelSD2YMSpan(1):gc%kernelSD2YMSpan(2), & 
                        gc%kernelSD2ZMSpan(1):gc%kernelSD2ZMSpan(2)  & 
                   )
      end do
      !$omp end parallel do
    else
      ! YZ
      roughness11Array => roughnessYYArray
      roughness22Array => roughnessZZArray
      !$omp parallel do schedule( dynamic, 1 ) & 
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvatureBandwidth )       &
      !$omp reduction( +:curvature1 )          &
      !$omp reduction( +:curvature2 )          &
      !$omp firstprivate( kernelSDY )          &
      !$omp firstprivate( kernelSDZ )          &
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins
  
        ! Assign gc pointer 
        gc => activeGridCells(n)
  
        if ( ( any( curvatureBandwidth( :, n ) .lt. fZERO ) ) .or. & 
             ( any( curvatureBandwidth( :, n ) /= curvatureBandwidth( :, n ) ) ) ) cycle

        ! Set kernels        
        call this%SetKernelSD2D( gc, kernelSDY, kernelSDZ, curvatureBandwidth( :, n ) )

        ! Compute curvature
        curvature1( &
                gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
            ) = curvature1( &
                gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
            ) + this%histogramCounts(                                &
            !) + this%histogram%counts(                               &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                        gc%kernelSD1XMSpan(1):gc%kernelSD1XMSpan(2), &
                        gc%kernelSD1YMSpan(1):gc%kernelSD1YMSpan(2), & 
                        gc%kernelSD1ZMSpan(1):gc%kernelSD1ZMSpan(2)  & 
                   )
        ! Compute curvature
        curvature2( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) = curvature2( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) + this%histogramCounts(                                &
            !) + this%histogram%counts(                               &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD2Matrix(&
                        gc%kernelSD2XMSpan(1):gc%kernelSD2XMSpan(2), &
                        gc%kernelSD2YMSpan(1):gc%kernelSD2YMSpan(2), & 
                        gc%kernelSD2ZMSpan(1):gc%kernelSD2ZMSpan(2)  & 
                   )
      end do
      !$omp end parallel do
    end if
    call kernelSDX%ResetMatrix()
    call kernelSDY%ResetMatrix()
    call kernelSDZ%ResetMatrix()

    curvature1 = curvature1/this%histogram%binVolume
    curvature2 = curvature2/this%histogram%binVolume
    ! Matrix from curvature kernels is delta**2*KernelVMatrix
    curvature1 = curvature1/( this%binSize(this%idDim1)**fTWO )
    curvature2 = curvature2/( this%binSize(this%idDim2)**fTWO )

    ! Product curvatures
    curvature12 = curvature1*curvature2
    curvature1  = curvature1*curvature1
    curvature2  = curvature2*curvature2

    ! Net roughness
    roughness11Array  = fZERO 
    roughness22Array  = fZERO 
    netRoughnessArray = fZERO
    if ( this%useGlobalSmoothing ) then
      ! If global smoothing, then the integral is over 
      ! the whole domain and kernels are considered to 
      ! be isotropic 
      roughness11Array(:) = sum(curvature1)
      roughness22Array(:) = sum(curvature2)
      netRoughnessArray   = roughness11Array + roughness22Array + fTWO*sum(curvature12)

      ! Deallocate
      deallocate( curvature1  ) 
      deallocate( curvature2  ) 
      deallocate( curvature12 ) 

      ! Done
      return
    end if

    ! Continue to local form !

    ! Initialize kernelSigma
    call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )

    smoothingCarrier  = fZERO
    if ( .not. this%isotropic ) then
      ! Anisotropic  
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvature1 )               &
      !$omp shared( curvature2 )               &
      !$omp shared( curvature12 )              &
      !$omp shared( roughness11Array )         &
      !$omp shared( roughness22Array )         &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp shared( netRoughnessArray )        & 
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp private( n )                       &
      !$omp private( gc )
      do n = 1, this%nComputeBins

        gc => activeGridCells(n)
        smoothingCarrier(:) = kernelSigmaSupportScale(n)
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Directional roughness
        roughness11Array( n ) = sum(&
            curvature1(&
                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
            )*gc%kernelSigmaMatrix(&
                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) ))
        roughness22Array( n ) = sum(&
            curvature2(&
                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
            )*gc%kernelSigmaMatrix(&
                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) ))
        ! Net roughness
        netRoughnessArray( n )  = fTWO*sum(&
            curvature12(&
                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2),     &
                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2),     & 
                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)      & 
            )*gc%kernelSigmaMatrix(&
                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2),     &
                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2),     & 
                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) + & 
              fTWO*sqrt(roughness11Array(n)*roughness22Array(n)) 
      end do
      !$omp end parallel do 
    else
      ! Isotropic
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvature1 )               &
      !$omp shared( curvature2 )               &
      !$omp shared( curvature12 )              &
      !$omp shared( roughness11Array )         &
      !$omp shared( roughness22Array )         &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp shared( netRoughnessArray )        & 
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp private( n )                       &
      !$omp private( gc )
      do n = 1, this%nComputeBins

        gc => activeGridCells(n)
        smoothingCarrier(:) = kernelSigmaSupportScale(n)
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Directional roughness
        roughness11Array( n ) = sum(&
            curvature1(&
                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
            )*gc%kernelSigmaMatrix(&
                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughness22Array( n ) = sum(&
            curvature2(&
                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
            )*gc%kernelSigmaMatrix(&
                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) ))
        ! Net roughness
        netRoughnessArray( n )  = fTWO*sum(&
            curvature12(&
                gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2),     &
                gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2),     & 
                gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)      & 
            )*gc%kernelSigmaMatrix(&
                gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2),     &
                gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2),     & 
                gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) + & 
              roughness11Array(n) + roughness22Array(n) 
      end do
      !$omp end parallel do 
    end if 
    call kernelSigma%ResetMatrix() 

    ! Deallocate
    deallocate( curvature1  ) 
    deallocate( curvature2  ) 
    deallocate( curvature12 ) 

    ! Done
    return

  end subroutine prComputeNetRoughness2D


  ! NET ROUGHNESS
  ! net roughness
  ! 3D
  subroutine prComputeNetRoughness3D( this, activeGridCells, curvatureBandwidth, &
                           roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                     netRoughnessArray, kernelSigmaSupportScale, &
                                                 kernelSDX, kernelSDY, kernelSDZ ) 
    !-----------------------------------------------------------------------------
    ! Net roughness in 3D
    !
    !  - Eq. 13c in Sole-Mari et al. (2019)
    !
    !-----------------------------------------------------------------------------
    ! Specifications 
    !-----------------------------------------------------------------------------
    !input
    class( GridProjectedKDEType ), target :: this
    type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
    real(fp), dimension(:,:), intent(in)                   :: curvatureBandwidth
    real(fp), dimension(:), intent(in)                     :: kernelSigmaSupportScale
    type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
    type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
    type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
    ! out
    real(fp), dimension(:), intent(inout), target   :: roughnessXXArray
    real(fp), dimension(:), intent(inout), target   :: roughnessYYArray
    real(fp), dimension(:), intent(inout), target   :: roughnessZZArray
    real(fp), dimension(:), intent(inout)           :: netRoughnessArray
    ! local 
    type( GridCellType ), pointer                   :: gc => null()
    real(fp), dimension(:,:,:), allocatable         ::  curvatureX
    real(fp), dimension(:,:,:), allocatable         ::  curvatureY
    real(fp), dimension(:,:,:), allocatable         ::  curvatureZ
    real(fp), dimension(:,:,:), allocatable         :: curvatureXY
    real(fp), dimension(:,:,:), allocatable         :: curvatureXZ
    real(fp), dimension(:,:,:), allocatable         :: curvatureYZ
    type( KernelMultiGaussianType )                 :: kernelSigma
    real(fp), dimension(3)                          :: smoothingCarrier
    integer :: n 
    real(fp) :: roughnessXY, roughnessXZ, roughnessYZ 
    !------------------------------------------------------------------------------
    allocate(  curvatureX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate(  curvatureY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate(  curvatureZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate( curvatureXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate( curvatureXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate( curvatureYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !------------------------------------------------------------------------------

    ! Initialize 
    curvatureX(:,:,:)  = fZERO
    curvatureY(:,:,:)  = fZERO
    curvatureZ(:,:,:)  = fZERO

    ! Curvatures, kappa
    ! Computed one at time for compatibility 
    ! with intel fortran compiler. Somehow reduction 
    ! of multiple matrices in one loop is not so stable
    ! and may lead to stack memory errors.
    !$omp parallel                           & 
    !$omp default( none )                    &
    !$omp shared( this )                     &
    !$omp shared( activeGridCells )          &
    !$omp shared( curvatureBandwidth )       & 
    !$omp firstprivate( kernelSDX )          &
    !$omp firstprivate( kernelSDY )          &
    !$omp firstprivate( kernelSDZ )          &
    !$omp reduction( +:curvatureX )          &
    !$omp reduction( +:curvatureY )          &
    !$omp reduction( +:curvatureZ )          
    ! CX
    !$omp do schedule( dynamic, 1 )  &
    !$omp private( n )               &
    !$omp private( gc )                       
    do n = 1, this%nComputeBins
  
      ! Assign gc pointer 
      gc => activeGridCells(n)
  
      ! Set kernel 
      call this%SetKernelSD(gc, kernelSDX, curvatureBandwidth( 1, n ), &
                                             this%kernelSDXDatabase, 1 )
      ! Compute curvatures
      curvatureX( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) = curvatureX( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) + this%histogram%counts(                             &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDMatrix(  &
                      gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                      gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                      gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
              )
    end do
    !$omp end do nowait
    ! CY
    !$omp do schedule( dynamic, 1 )  & 
    !$omp private( n )               &
    !$omp private( gc )               
    do n = 1, this%nComputeBins
  
      ! Assign gc pointer 
      gc => activeGridCells(n)
  
      ! Set kernel 
      call this%SetKernelSD(gc, kernelSDY, curvatureBandwidth( 2, n ), &
                                             this%kernelSDYDatabase, 2 )
      ! Compute curvatures
      curvatureY( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) = curvatureY( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) + this%histogram%counts(                             &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDMatrix(  &
                      gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                      gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                      gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
              )
    end do
    !$omp end do nowait
    ! CZ
    !$omp do schedule( dynamic, 1 )  & 
    !$omp private( n )               &
    !$omp private( gc )               
    do n = 1, this%nComputeBins
  
      ! Assign gc pointer 
      gc => activeGridCells(n)
  
      ! Set kernel 
      call this%SetKernelSD(gc, kernelSDZ, curvatureBandwidth( 3, n ), &
                                             this%kernelSDZDatabase, 3 )
      ! Compute curvatures
      curvatureZ( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) = curvatureZ( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) + this%histogram%counts(                             &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDMatrix(  &
                      gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                      gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                      gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
              )
    end do
    !$omp end do

    !$omp end parallel
    call kernelSDX%ResetMatrix()
    call kernelSDY%ResetMatrix()
    call kernelSDZ%ResetMatrix()

    curvatureX = curvatureX/this%histogram%binVolume
    curvatureY = curvatureY/this%histogram%binVolume
    curvatureZ = curvatureZ/this%histogram%binVolume
    ! Matrix from curvature kernels is delta**2*KernelVMatrix
    curvatureX = curvatureX/( this%binSize(1)**fTWO )
    curvatureY = curvatureY/( this%binSize(2)**fTWO )
    curvatureZ = curvatureZ/( this%binSize(3)**fTWO )

    ! Compute curvatures product
    curvatureX  = curvatureX*curvatureX
    curvatureY  = curvatureY*curvatureY
    curvatureZ  = curvatureZ*curvatureZ
    curvatureXY = curvatureX*curvatureY
    curvatureXZ = curvatureX*curvatureZ
    curvatureYZ = curvatureY*curvatureZ

    ! Net roughness
    roughnessXXArray  = fZERO 
    roughnessYYArray  = fZERO 
    roughnessZZArray  = fZERO 
    netRoughnessArray = fZERO
    if ( this%useGlobalSmoothing ) then
      ! If global smoothing, then the integral is over 
      ! the whole domain and kernels are considered to 
      ! be isotropic 
      roughnessXXArray(:) = sum(curvatureX)
      roughnessYYArray(:) = sum(curvatureY)
      roughnessZZArray(:) = sum(curvatureZ)
      netRoughnessArray   = &
        roughnessXXArray + roughnessYYArray + roughnessZZArray + & 
        fTWO*( sum(curvatureXY) + sum(curvatureXZ) + sum(curvatureYZ) )

      ! Deallocate
      deallocate(  curvatureX )
      deallocate(  curvatureY )
      deallocate(  curvatureZ )
      deallocate( curvatureXY )
      deallocate( curvatureXZ )
      deallocate( curvatureYZ )

      ! Done
      return
    end if

    ! Continue to local form !

    ! Initialize kernelSigma
    call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )

    roughnessXY       = fZERO
    roughnessXZ       = fZERO
    roughnessYZ       = fZERO
    smoothingCarrier  = fZERO
    if ( .not. this%isotropic ) then
      ! Anisotropic  
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvatureX )               &
      !$omp shared( curvatureY )               &
      !$omp shared( curvatureZ )               &
      !$omp shared( curvatureXY )              &
      !$omp shared( curvatureYZ )              &
      !$omp shared( curvatureXZ )              &
      !$omp shared( roughnessXXArray )         &
      !$omp shared( roughnessYYArray )         &
      !$omp shared( roughnessZZArray )         &
      !$omp shared( netRoughnessArray )        &
      !$omp firstprivate( roughnessXY )        &
      !$omp firstprivate( roughnessXZ )        &
      !$omp firstprivate( roughnessYZ )        &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins

        gc => activeGridCells(n)
        smoothingCarrier(:) = kernelSigmaSupportScale(n)
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Directional roughness 
        roughnessXXArray( n ) = sum(&
          curvatureX(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYYArray( n ) = sum(&
          curvatureY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessZZArray( n ) = sum(&
          curvatureZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXY = sum(&
          curvatureXY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXZ = sum(&
          curvatureXZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYZ = sum(&
          curvatureYZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        ! Net roughness
        netRoughnessArray( n ) = fTHREE*( roughnessXXArray(n)*roughnessYYArray(n)*roughnessZZArray(n) )**(fONE/fTHREE) 
        if ( roughnessYZ.gt.fZERO ) netRoughnessArray(n) = netRoughnessArray(n) + & 
            fTWO*roughnessYZ*( roughnessXXArray(n)**fTWO/roughnessYYArray(n)/roughnessZZArray(n) )**(fONE/fSIX)
        if ( roughnessXZ.gt.fZERO ) netRoughnessArray(n) = netRoughnessArray(n) + & 
            fTWO*roughnessXZ*( roughnessYYArray(n)**fTWO/roughnessXXArray(n)/roughnessZZArray(n) )**(fONE/fSIX)
        if ( roughnessXY.gt.fZERO ) netRoughnessArray(n) = netRoughnessArray(n) + & 
            fTWO*roughnessXY*( roughnessZZArray(n)**fTWO/roughnessXXArray(n)/roughnessYYArray(n) )**(fONE/fSIX)
      end do
      !$omp end parallel do
    else
      ! Isotropic
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvatureX )               &
      !$omp shared( curvatureY )               &
      !$omp shared( curvatureZ )               &
      !$omp shared( curvatureXY )              &
      !$omp shared( curvatureYZ )              &
      !$omp shared( curvatureXZ )              &
      !$omp shared( roughnessXXArray )         &
      !$omp shared( roughnessYYArray )         &
      !$omp shared( roughnessZZArray )         &
      !$omp shared( netRoughnessArray )        &
      !$omp firstprivate( roughnessXY )        &
      !$omp firstprivate( roughnessXZ )        &
      !$omp firstprivate( roughnessYZ )        &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins

        gc => activeGridCells(n)
        smoothingCarrier(:) = kernelSigmaSupportScale(n)
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Directional roughness
        roughnessXXArray( n ) = sum(&
          curvatureX(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYYArray( n ) = sum(&
          curvatureY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessZZArray( n ) = sum(&
          curvatureZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXY = sum(&
          curvatureXY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXZ = sum(&
          curvatureXZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYZ = sum(&
          curvatureYZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        ! Net roughness
        netRoughnessArray( n ) = roughnessXXArray(n) + fTWO*roughnessXY + fTWO*roughnessXZ + &
                                 roughnessYYArray(n) + fTWO*roughnessYZ + roughnessZZArray(n)
      end do
      !$omp end parallel do
    end if 
    call kernelSigma%ResetMatrix() 

    ! Deallocate
    deallocate(  curvatureX )
    deallocate(  curvatureY )
    deallocate(  curvatureZ )
    deallocate( curvatureXY )
    deallocate( curvatureXZ )
    deallocate( curvatureYZ )


    ! Done
    return


  end subroutine prComputeNetRoughness3D


  ! NET ROUGHNESS WITH INDEPENDENT CURVATURE PARALLEL LOOPS
  ! net roughness
  ! 3D
  subroutine prComputeNetRoughness3DIndep( this, activeGridCells, curvatureBandwidth, &
                                roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                          netRoughnessArray, kernelSigmaSupportScale, &
                                                      kernelSDX, kernelSDY, kernelSDZ ) 
    !-----------------------------------------------------------------------------
    ! Net roughness in 3D
    !
    !  - Eq. 13c in Sole-Mari et al. (2019)
    !
    !-----------------------------------------------------------------------------
    ! Specifications 
    !-----------------------------------------------------------------------------
    !input
    class( GridProjectedKDEType ), target :: this
    type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
    real(fp), dimension(:,:), intent(in)                   :: curvatureBandwidth
    real(fp), dimension(:), intent(in)                     :: kernelSigmaSupportScale
    type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
    type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
    type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
    ! out
    real(fp), dimension(:), intent(inout), target   :: roughnessXXArray
    real(fp), dimension(:), intent(inout), target   :: roughnessYYArray
    real(fp), dimension(:), intent(inout), target   :: roughnessZZArray
    real(fp), dimension(:), intent(inout)           :: netRoughnessArray
    ! local 
    type( GridCellType ), pointer                   :: gc => null()
    real(fp), dimension(:,:,:), allocatable         ::  curvatureX
    real(fp), dimension(:,:,:), allocatable         ::  curvatureY
    real(fp), dimension(:,:,:), allocatable         ::  curvatureZ
    real(fp), dimension(:,:,:), allocatable         :: curvatureXY
    real(fp), dimension(:,:,:), allocatable         :: curvatureXZ
    real(fp), dimension(:,:,:), allocatable         :: curvatureYZ
    type( KernelMultiGaussianType )                 :: kernelSigma
    real(fp), dimension(3)                          :: smoothingCarrier
    integer :: n 
    real(fp) :: roughnessXY, roughnessXZ, roughnessYZ 
    !------------------------------------------------------------------------------
    allocate(  curvatureX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate(  curvatureY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate(  curvatureZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate( curvatureXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate( curvatureXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    allocate( curvatureYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
    !------------------------------------------------------------------------------

    ! Initialize 
    curvatureX(:,:,:)  = fZERO
    curvatureY(:,:,:)  = fZERO
    curvatureZ(:,:,:)  = fZERO

    ! Curvatures, kappa
    ! Computed one at time for compatibility 
    ! with intel fortran compiler. Somehow reduction 
    ! of multiple matrices in one loop is not so stable
    ! and may lead to stack memory errors.
    ! CX
    !$omp parallel do schedule( dynamic, 1 ) &
    !$omp default( none )                    &
    !$omp shared( this )                     &
    !$omp shared( activeGridCells )          &
    !$omp shared( curvatureBandwidth )       & 
    !$omp firstprivate( kernelSDX )          &
    !$omp reduction( +:curvatureX )          &          
    !$omp private( n )                       &
    !$omp private( gc )                       
    do n = 1, this%nComputeBins
  
      ! Assign gc pointer 
      gc => activeGridCells(n)
  
      ! Set kernel 
      call this%SetKernelSD(gc, kernelSDX, curvatureBandwidth( 1, n ), &
                                             this%kernelSDXDatabase, 1 )
      ! Compute curvatures
      curvatureX( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) = curvatureX( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) + this%histogram%counts(                             &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDMatrix(  &
                      gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                      gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                      gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
              )
    end do
    !$omp end parallel do  
    ! CY
    !$omp parallel do schedule( dynamic, 1 ) &
    !$omp default( none )                    &
    !$omp shared( this )                     &
    !$omp shared( activeGridCells )          &
    !$omp shared( curvatureBandwidth )       & 
    !$omp firstprivate( kernelSDY )          &
    !$omp reduction( +:curvatureY )          &
    !$omp private( n )                       &
    !$omp private( gc )                       
    do n = 1, this%nComputeBins
  
      ! Assign gc pointer 
      gc => activeGridCells(n)
  
      ! Set kernel 
      call this%SetKernelSD(gc, kernelSDY, curvatureBandwidth( 2, n ), &
                                             this%kernelSDYDatabase, 2 )
      ! Compute curvatures
      curvatureY( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) = curvatureY( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) + this%histogram%counts(                             &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDMatrix(  &
                      gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                      gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                      gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
              )
    end do
    !$omp end parallel do
    ! CZ
    !$omp parallel do schedule( dynamic, 1 ) &
    !$omp default( none )                    &
    !$omp shared( this )                     &
    !$omp shared( activeGridCells )          &
    !$omp shared( curvatureBandwidth )       & 
    !$omp firstprivate( kernelSDZ )          &
    !$omp reduction( +:curvatureZ )          &
    !$omp private( n )                       &
    !$omp private( gc )               
    do n = 1, this%nComputeBins
  
      ! Assign gc pointer 
      gc => activeGridCells(n)
  
      ! Set kernel 
      call this%SetKernelSD(gc, kernelSDZ, curvatureBandwidth( 3, n ), &
                                             this%kernelSDZDatabase, 3 )
      ! Compute curvatures
      curvatureZ( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) = curvatureZ( &
              gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
              gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
              gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
          ) + this%histogram%counts(                             &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSDMatrix(  &
                      gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                      gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                      gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
              )
    end do
    !$omp end parallel do
    call kernelSDX%ResetMatrix()
    call kernelSDY%ResetMatrix()
    call kernelSDZ%ResetMatrix()

    curvatureX = curvatureX/this%histogram%binVolume
    curvatureY = curvatureY/this%histogram%binVolume
    curvatureZ = curvatureZ/this%histogram%binVolume
    ! Matrix from curvature kernels is delta**2*KernelVMatrix
    curvatureX = curvatureX/( this%binSize(1)**fTWO )
    curvatureY = curvatureY/( this%binSize(2)**fTWO )
    curvatureZ = curvatureZ/( this%binSize(3)**fTWO )

    ! Compute curvatures product
    curvatureX  = curvatureX*curvatureX
    curvatureY  = curvatureY*curvatureY
    curvatureZ  = curvatureZ*curvatureZ
    curvatureXY = curvatureX*curvatureY
    curvatureXZ = curvatureX*curvatureZ
    curvatureYZ = curvatureY*curvatureZ

    ! Net roughness
    roughnessXXArray  = fZERO 
    roughnessYYArray  = fZERO 
    roughnessZZArray  = fZERO 
    netRoughnessArray = fZERO
    if ( this%useGlobalSmoothing ) then
      ! If global smoothing, then the integral is over 
      ! the whole domain and kernels are considered to 
      ! be isotropic 
      roughnessXXArray(:) = sum(curvatureX)
      roughnessYYArray(:) = sum(curvatureY)
      roughnessZZArray(:) = sum(curvatureZ)
      netRoughnessArray   = &
        roughnessXXArray + roughnessYYArray + roughnessZZArray + & 
        fTWO*( sum(curvatureXY) + sum(curvatureXZ) + sum(curvatureYZ) )

      ! Deallocate
      deallocate(  curvatureX )
      deallocate(  curvatureY )
      deallocate(  curvatureZ )
      deallocate( curvatureXY )
      deallocate( curvatureXZ )
      deallocate( curvatureYZ )

      ! Done
      return
    end if

    ! Continue to local form !

    ! Initialize kernelSigma
    call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )

    roughnessXY       = fZERO
    roughnessXZ       = fZERO
    roughnessYZ       = fZERO
    smoothingCarrier  = fZERO
    if ( .not. this%isotropic ) then
      ! Anisotropic  
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvatureX )               &
      !$omp shared( curvatureY )               &
      !$omp shared( curvatureZ )               &
      !$omp shared( curvatureXY )              &
      !$omp shared( curvatureYZ )              &
      !$omp shared( curvatureXZ )              &
      !$omp shared( roughnessXXArray )         &
      !$omp shared( roughnessYYArray )         &
      !$omp shared( roughnessZZArray )         &
      !$omp shared( netRoughnessArray )        &
      !$omp firstprivate( roughnessXY )        &
      !$omp firstprivate( roughnessXZ )        &
      !$omp firstprivate( roughnessYZ )        &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins

        gc => activeGridCells(n)
        smoothingCarrier(:) = kernelSigmaSupportScale(n)
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Directional roughness 
        roughnessXXArray( n ) = sum(&
          curvatureX(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYYArray( n ) = sum(&
          curvatureY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessZZArray( n ) = sum(&
          curvatureZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXY = sum(&
          curvatureXY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXZ = sum(&
          curvatureXZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYZ = sum(&
          curvatureYZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        ! Net roughness
        netRoughnessArray( n ) = fTHREE*( roughnessXXArray(n)*roughnessYYArray(n)*roughnessZZArray(n) )**(fONE/fTHREE) 
        if ( roughnessYZ.gt.fZERO ) netRoughnessArray(n) = netRoughnessArray(n) + & 
            fTWO*roughnessYZ*( roughnessXXArray(n)**fTWO/roughnessYYArray(n)/roughnessZZArray(n) )**(fONE/fSIX)
        if ( roughnessXZ.gt.fZERO ) netRoughnessArray(n) = netRoughnessArray(n) + & 
            fTWO*roughnessXZ*( roughnessYYArray(n)**fTWO/roughnessXXArray(n)/roughnessZZArray(n) )**(fONE/fSIX)
        if ( roughnessXY.gt.fZERO ) netRoughnessArray(n) = netRoughnessArray(n) + & 
            fTWO*roughnessXY*( roughnessZZArray(n)**fTWO/roughnessXXArray(n)/roughnessYYArray(n) )**(fONE/fSIX)
      end do
      !$omp end parallel do
    else
      ! Isotropic
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvatureX )               &
      !$omp shared( curvatureY )               &
      !$omp shared( curvatureZ )               &
      !$omp shared( curvatureXY )              &
      !$omp shared( curvatureYZ )              &
      !$omp shared( curvatureXZ )              &
      !$omp shared( roughnessXXArray )         &
      !$omp shared( roughnessYYArray )         &
      !$omp shared( roughnessZZArray )         &
      !$omp shared( netRoughnessArray )        &
      !$omp firstprivate( roughnessXY )        &
      !$omp firstprivate( roughnessXZ )        &
      !$omp firstprivate( roughnessYZ )        &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp private( n )                       &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins

        gc => activeGridCells(n)
        smoothingCarrier(:) = kernelSigmaSupportScale(n)
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Directional roughness
        roughnessXXArray( n ) = sum(&
          curvatureX(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYYArray( n ) = sum(&
          curvatureY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessZZArray( n ) = sum(&
          curvatureZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXY = sum(&
          curvatureXY(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessXZ = sum(&
          curvatureXZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        roughnessYZ = sum(&
          curvatureYZ(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 
        ! Net roughness
        netRoughnessArray( n ) = roughnessXXArray(n) + fTWO*roughnessXY + fTWO*roughnessXZ + &
                                 roughnessYYArray(n) + fTWO*roughnessYZ + roughnessZZArray(n)
      end do
      !$omp end parallel do
    end if 
    call kernelSigma%ResetMatrix() 

    ! Deallocate
    deallocate(  curvatureX )
    deallocate(  curvatureY )
    deallocate(  curvatureZ )
    deallocate( curvatureXY )
    deallocate( curvatureXZ )
    deallocate( curvatureYZ )


    ! Done
    return


  end subroutine prComputeNetRoughness3DIndep



  subroutine prInitializeKernelDatabaseFlat( this,  &
                      minHOverDelta, maxHOverDelta, &
                deltaHOverDelta, logKernelDatabase, &
                         kernelRange, kernelSDRange )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target :: this
    ! input
    real(fp), intent(in) :: minHOverDelta
    real(fp), intent(in) :: maxHOverDelta
    real(fp), intent(in) :: deltaHOverDelta
    logical , intent(in), optional :: logKernelDatabase
    integer , intent(in), optional :: kernelRange
    integer , intent(in), optional :: kernelSDRange
    ! local
    real(fp), dimension(3) :: inputSmoothing
    real(fp), dimension(:), allocatable :: hOverDelta
    integer :: nDelta
    integer :: i, n, m, o, dbi
    logical :: localLogDatabase
    integer :: localKernelRange
    integer :: localKernelSDRange
    ! Mem control
    real(fp) :: kernelMatrixMemory
    real(fp) :: kernelDBMemory    
    real(fp) :: kernelSDDBMemory  
    ! Time control 
    integer  :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
    real(fp) :: elapsedTime
    character(len=30) :: outfmt
    !------------------------------------------------------------------------------

    ! Needs sanity check for input parameters !
    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, '(2X,A)' ) 'Initializing Kernels Database'
    end if 

    ! Default parameters
    if ( present( logKernelDatabase ) ) then 
      localLogDatabase = logKernelDatabase
    else
      localLogDatabase = defaultLogKernelDatabase
    end if 
    ! Kernel ranges 
    if ( present( kernelRange ) )  then 
      localKernelRange = kernelRange
    else 
      localKernelRange = defaultKernelRange
    end if
    if ( present( kernelSDRange ) ) then 
      localKernelSDRange = kernelSDRange
    else 
      localKernelSDRange = defaultKernelSDRange
    end if 

    ! In the meantime a single nDelta, 
    ! it could be any discretization
    if ( localLogDatabase ) then
      ! LOG FORM
      nDelta      = ceiling( log10( maxHOverDelta/minHOverDelta )/log10( 1 + deltaHOverDelta ) ) + 1
      allocate( hOverDelta( nDelta ) )
      hOverDelta = prGenerateLogSpaceData( minHOverDelta, maxHOverDelta, nDelta )
      ! Assign indexes interfaces
      this%ComputeKernelDatabaseFlatIndexes => prComputeKernelDatabaseFlatIndexesLog
      this%ComputeKernelDatabaseIndexes     => prComputeKernelDatabaseIndexesLog ! meanwhile for SD's
      this%deltaHOverDelta = log( hOverDelta(2)/hOverDelta(1) ) ! Fix this inconsistency, is overwritten
    else 
      ! LINEAR FORM
      nDelta      = floor( ( maxHOverDelta - minHOverDelta )/deltaHOverDelta )
      allocate( hOverDelta( nDelta ) )
      hOverDelta = [ (minHOverDelta + i*deltaHOverDelta, i=0, nDelta ) ]
      ! Assign indexes interface
      this%ComputeKernelDatabaseFlatIndexes => prComputeKernelDatabaseFlatIndexesLinear
      this%ComputeKernelDatabaseIndexes     => prComputeKernelDatabaseIndexesLinear ! meanwhile for SD's
    end if

    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, '(3X,A)') 'Database parameters'
      write( this%outFileUnit, '(3X,A)') '-------------------'
      outfmt = '(3X,A,(1X,es18.9e3))'
      write( this%outFileUnit, outfmt ) '- minHOverDelta         :', minHOverDelta
      write( this%outFileUnit, outfmt ) '- maxHOverDelta         :', maxHOverDelta
      write( this%outFileUnit, outfmt ) '- deltaHOverDelta       :', deltaHOverDelta
      write( this%outFileUnit, '(3X,A)') '-------------------'
    end if 

    ! Assign to the object
    ! Temporarilly the same value for each axis
    this%nDeltaHOverDelta = nDelta

    ! Depending on the number of dimensions
    ! is the required kernel database.
    select case(nDim)
    !1D
    case(1)
      ! Allocate kernel databases
      allocate( this%kernelDatabaseFlat( nDelta, 1 ) )
      
      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(3X,A)' ) 'Database is for kernels 1D '
      end if 

      ! Assign kernelSD db pointer according 
      ! to determined direction
      select case(this%idDim1) 
        case (1)
          allocate( this%kernelSDXDatabase( nDelta ) )
          this%kernelSDDatabase1 => this%kernelSDXDatabase
        case (2)
          allocate( this%kernelSDYDatabase( nDelta ) )
          this%kernelSDDatabase1 => this%kernelSDYDatabase
        case (3)
          allocate( this%kernelSDZDatabase( nDelta ) )
          this%kernelSDDatabase1 => this%kernelSDZDatabase
      end select   

      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Kernel database
      kernelMatrixMemory = fZERO
      kernelDBMemory = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverDelta )               &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelRange )         &
      !$omp reduction( +:kernelDBMemory )      &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( n )                       &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        inputSmoothing(:) = 0
        inputSmoothing( this%idDim1 ) = hOverDelta(n)
        call this%kernelDatabaseFlat( n, 1 )%Initialize( &
          this%binSize, matrixRange=localKernelRange )
        call this%kernelDatabaseFlat( n, 1 )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelDatabaseFlat( n, 1 )%matrix )/1.0e6_fp
        kernelDBMemory     = kernelDBMemory + kernelMatrixMemory
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Second derivatives
      kernelMatrixMemory = fZERO
      kernelSDDBMemory = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverDelta )               &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelSDRange )       &
      !$omp reduction( +:kernelSDDBMemory )    &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( n )                       &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        inputSmoothing(:) = 0
        inputSmoothing( this%idDim1 ) = hOverDelta(n)
        ! 1 
        call this%kernelSDDatabase1( n )%Initialize(& 
          this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDDatabase1( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDDatabase1( n )%matrix )/1.0e6_fp
        kernelSDDBMemory   = kernelSDDBMemory + kernelMatrixMemory
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

    ! 2D
    case(2)
      allocate( this%kernelDatabaseFlat( nDelta*( nDelta + 1 )/2, 1 ) )

      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(3X,A)' ) 'Database is for kernels 2D '
      end if

      ! Assign kernelSD db pointers according 
      ! to determined directions
      select case(this%idDim1) 
        case (1)
          allocate( this%kernelSDXDatabase( nDelta ) )
          this%kernelSDDatabase1 => this%kernelSDXDatabase
        case (2)
          allocate( this%kernelSDYDatabase( nDelta ) )
          this%kernelSDDatabase1 => this%kernelSDYDatabase
        case (3)
          allocate( this%kernelSDZDatabase( nDelta ) )
          this%kernelSDDatabase1 => this%kernelSDZDatabase
      end select   
      select case(this%idDim2) 
        case (1)
          ! This should not be the case, x is always the first
          allocate( this%kernelSDXDatabase( nDelta ) )
          this%kernelSDDatabase2 => this%kernelSDXDatabase
        case (2)
          allocate( this%kernelSDYDatabase( nDelta ) )
          this%kernelSDDatabase2 => this%kernelSDYDatabase
        case (3)
          allocate( this%kernelSDZDatabase( nDelta ) )
          this%kernelSDDatabase2 => this%kernelSDZDatabase
      end select   

      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Kernel database
      kernelMatrixMemory = fZERO
      kernelDBMemory = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverDelta )               &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelRange )         &
      !$omp private( m, dbi )                  &
      !$omp reduction( +:kernelDBMemory )      &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( n )                       &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        do m = 1, min( n, nDelta )
          dbi = n*( n - 1 )/2 + m
          inputSmoothing(:) = fZERO
          inputSmoothing( this%idDim1 ) =  hOverDelta(n)
          inputSmoothing( this%idDim2 ) =  hOverDelta(m)
          call this%kernelDatabaseFlat( dbi, 1 )%Initialize( & 
            this%binSize, matrixRange=localKernelRange )
          call this%kernelDatabaseFlat( dbi, 1 )%SetupMatrix( inputSmoothing*this%binSize )
          kernelMatrixMemory = sizeof( this%kernelDatabaseFlat( dbi, 1 )%matrix )/1.0e6_fp
          kernelDBMemory = kernelDBMemory + kernelMatrixMemory
        end do
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Second derivatives
      kernelMatrixMemory = fZERO
      kernelSDDBMemory = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverDelta )               &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelSDRange )       &
      !$omp reduction( +:kernelSDDBMemory )    &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( n )                       &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        inputSmoothing(:) = 0
        inputSmoothing( this%idDim1 ) = hOverDelta(n)
        inputSmoothing( this%idDim2 ) = hOverDelta(n)
        ! 1 
        call this%kernelSDDatabase1( n )%Initialize(& 
            this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDDatabase1( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDDatabase1( n )%matrix )/1.0e6_fp
        kernelSDDBMemory   = kernelSDDBMemory + kernelMatrixMemory
        ! 2 
        call this%kernelSDDatabase2( n )%Initialize(& 
            this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDDatabase2( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDDatabase1( n )%matrix )/1.0e6_fp
        kernelSDDBMemory   = kernelSDDBMemory + kernelMatrixMemory
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

    ! 3D
    case(3)

      ! Allocate kernel databases
      allocate( this%kernelDatabaseFlat( nDelta*( nDelta + 1 )/2, nDelta ) )
      allocate( this%kernelSDXDatabase( nDelta ) )
      allocate( this%kernelSDYDatabase( nDelta ) )
      allocate( this%kernelSDZDatabase( nDelta ) )

      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(3X,A)' ) 'Database is for kernels 3D '
      end if

      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Kernel database
      kernelMatrixMemory = fZERO
      kernelDBMemory = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverDelta )               &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelRange )         &
      !$omp private( n, m, dbi )               &
      !$omp reduction( +:kernelDBMemory )      &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( o )                       &
      !$omp private( inputSmoothing )
      do o = 1, nDelta
        do n = 1, nDelta
          do m = 1, min( n, nDelta )
            dbi = n*( n - 1 )/2 + m
            inputSmoothing = (/ hOverDelta(n), hOverDelta(m), hOverDelta(o) /) 
            call this%kernelDatabaseFlat( dbi, o )%Initialize( & 
              this%binSize, matrixRange=localKernelRange )
            call this%kernelDatabaseFlat( dbi, o )%SetupMatrix( inputSmoothing*this%binSize )
            kernelMatrixMemory = sizeof( this%kernelDatabaseFlat( dbi, o )%matrix )/1.0e6_fp
            kernelDBMemory = kernelDBMemory + kernelMatrixMemory
          end do
        end do
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
      
      ! TIC
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      ! Second derivatives
      kernelMatrixMemory = fZERO
      kernelSDDBMemory = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( hOverDelta )               &
      !$omp shared( nDelta )                   &
      !$omp shared( localKernelSDRange )       &
      !$omp reduction( +:kernelSDDBMemory )    &
      !$omp private( kernelMatrixMemory )      &
      !$omp private( n )                       &
      !$omp private( inputSmoothing )
      do n = 1, nDelta
        inputSmoothing = (/ hOverDelta(n), hOverDelta(n), hOverDelta(n) /)
        ! X 
        call this%kernelSDXDatabase( n )%Initialize(& 
          this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDXDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDXDatabase( n )%matrix )/1.0e6_fp
        kernelSDDBMemory = kernelSDDBMemory + kernelMatrixMemory
        ! Y
        call this%kernelSDYDatabase( n )%Initialize(& 
          this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDYDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDYDatabase( n )%matrix )/1.0e6_fp
        kernelSDDBMemory = kernelSDDBMemory + kernelMatrixMemory
        ! Z
        call this%kernelSDZDatabase( n )%Initialize(& 
          this%binSize, matrixRange=localKernelSDRange )
        call this%kernelSDZDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )
        kernelMatrixMemory = sizeof( this%kernelSDZDatabase( n )%matrix )/1.0e6_fp
        kernelSDDBMemory = kernelSDDBMemory + kernelMatrixMemory
      end do
      !$omp end parallel do
      ! TOC
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

    end select

    ! Report db sizes
    if ( this%reportToOutUnit ) then 
      write(this%outFileUnit, '(3X,A,E15.1,A)')& 
        ' - Allocated memory Kernel DB    : ', kernelDBMemory, ' MB'
      write(this%outFileUnit, '(3X,A,E15.1,A)')& 
        ' - Allocated memory Kernel SD DB : ', kernelSDDBMemory, ' MB'
    end if 

    deallocate( hOverDelta )

    ! Done
    return


  end subroutine prInitializeKernelDatabaseFlat

  
  subroutine prDropKernelDatabase( this )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ) :: this
    !------------------------------------------------------------------------------

    if ( allocated( this%kernelDatabase     ) ) deallocate( this%kernelDatabase )
    if ( allocated( this%kernelDatabaseFlat ) ) deallocate( this%kernelDatabaseFlat )
    if ( allocated( this%kernelSDXDatabase  ) ) deallocate( this%kernelSDXDatabase )
    if ( allocated( this%kernelSDYDatabase  ) ) deallocate( this%kernelSDYDatabase )
    if ( allocated( this%kernelSDZDatabase  ) ) deallocate( this%kernelSDZDatabase )

    ! Done
    return

  end subroutine prDropKernelDatabase


  subroutine prGenerateVectorCoordinates( this ) 
  !----------------------------------------------------------------------------------------
  ! For a given grid and dimensionality, fill vectors with cell center coordinates
  !----------------------------------------------------------------------------------------
  ! Specifications 
  !----------------------------------------------------------------------------------------
  implicit none
  ! input
  class( GridProjectedKDEType ) :: this
  ! local 
  integer :: nd, m, idbin
  !----------------------------------------------------------------------------------------


    do nd=1,3
      if ( this%dimensionMask(nd).eq.1 ) then
        select case(nd)
        case(1)
          if ( allocated( this%coordinatesX ) ) deallocate( this%coordinatesX )
          allocate( this%coordinatesX(this%nBins(nd)) )
          do m = 1, this%nBins(nd)
             idbin = m+this%deltaBinsOrigin(nd)
             this%coordinatesX(m) = (real(idbin,fp) + 0.5_fp)*this%binSize(nd) + this%domainOrigin(nd)
          end do 
        case(2)
          if ( allocated( this%coordinatesY ) ) deallocate( this%coordinatesY )
          allocate( this%coordinatesY(this%nBins(nd)) )
          do m = 1, this%nBins(nd)
             idbin = m+this%deltaBinsOrigin(nd)
             this%coordinatesY(m) = (real(idbin,fp) + 0.5_fp)*this%binSize(nd) + this%domainOrigin(nd)
          end do 
        case(3)
          if ( allocated( this%coordinatesZ ) ) deallocate( this%coordinatesZ )
          allocate( this%coordinatesZ(this%nBins(nd)) )
          do m = 1, this%nBins(nd)
             idbin = m+this%deltaBinsOrigin(nd)
             this%coordinatesZ(m) = (real(idbin,fp) + 0.5_fp)*this%binSize(nd) + this%domainOrigin(nd)
          end do 
        end select 
      end if 
    end do


  end subroutine prGenerateVectorCoordinates


  ! Density computation manager 
  subroutine prComputeDensity( this, dataPoints, nOptimizationLoops, &
                                     outputFileName, outputFileUnit, &
                               outputColumnFormat, outputDataFormat, &
                     outputDataId, outputDataIdVal, particleGroupId, &
              persistentKernelDatabase, exportOptimizationVariables, & 
                                   skipErrorConvergence, unitVolume, &
                              scalingFactor, histogramScalingFactor, &
                                                  computeRawDensity, &
              weightedHistogram, weights, onlyHistogram, exactPoint, &
                                relativeErrorConvergence, isotropic, &
                                                 useGlobalSmoothing, & 
                                                 histogramBinFormat, & 
                                                    binSizeFraction, &
                                           generateVectorCoordinates ) 
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    ! input
    class( GridProjectedKDEType ), target        :: this
    real(fp), dimension(:,:), intent(in)         :: dataPoints
    integer, intent(in), optional                :: nOptimizationLoops
    character(len=*), intent(in), optional       :: outputFileName
    integer, intent(in), optional                :: outputFileUnit
    integer, intent(in), optional                :: outputColumnFormat
    integer, intent(in), optional                :: outputDataFormat
    integer, intent(in), optional                :: outputDataId
    real(fp), intent(in), optional               :: outputDataIdVal
    integer, intent(in), optional                :: particleGroupId
    logical, intent(in), optional                :: persistentKernelDatabase
    logical, intent(in), optional                :: exportOptimizationVariables
    logical, intent(in), optional                :: skipErrorConvergence
    logical, intent(in), optional                :: unitVolume
    real(fp), intent(in), optional               :: scalingFactor
    real(fp), intent(in), optional               :: histogramScalingFactor
    logical, intent(in), optional                :: computeRawDensity
    logical, intent(in), optional                :: weightedHistogram
    logical, intent(in), optional                :: onlyHistogram
    logical, intent(in), optional                :: exactPoint
    real(fp), dimension(:), intent(in), optional :: weights
    real(fp), intent(in), optional               :: relativeErrorConvergence
    logical, intent(in), optional                :: isotropic
    logical, intent(in), optional                :: useGlobalSmoothing
    integer, intent(in), optional                :: histogramBinFormat 
    real(fp), intent(in), optional               :: binSizeFraction
    logical, intent(in), optional                :: generateVectorCoordinates
    ! local 
    logical               :: persistKDB
    logical               :: locExportOptimizationVariables
    logical               :: locSkipErrorConvergence
    logical               :: locUnitVolume
    real(fp)              :: locScalingFactor
    logical               :: locScaleHistogram
    logical               :: locComputeRawDensity
    real(fp)              :: locHistogramScalingFactor
    logical               :: locWeightedHistogram
    logical               :: locOnlyHistogram 
    logical               :: locExactPoint 
    logical               :: locIsotropic
    integer               :: localNOptimizationLoops
    integer               :: locOutputColumnFormat
    integer               :: locOutputDataFormat
    real(fp)              :: locBinSizeFraction
    integer, dimension(2) :: dataPointsShape
    real(fp)              :: locRelativeErrorConvergence
    character(len=16)     :: timeChar
    character(len=16)     :: spcChar
    integer               :: nd
    ! For determination of sub grid
    real(fp), dimension(3) :: minCoords
    real(fp), dimension(3) :: minSubGridCoords
    real(fp), dimension(3) :: maxCoords
    real(fp), dimension(3) :: maxSubGridCoords
    real(fp), dimension(3) :: deltaCoords 
    real(fp), dimension(3) :: subGridSize
    integer , dimension(3) :: subGridNBins
    real(fp), dimension(3) :: subGridOrigin
    integer , dimension(3) :: subGridOriginIndexes
    integer , dimension(3) :: subGridLimitIndexes
    ! sliced reconstruction
    real(fp), dimension(:,:,:), pointer :: slicedDensity
    integer                             :: sliceId
    integer                             :: activeBins
    ! clock
    real(fp)               :: elapsedTime
    integer                :: clockCountStart, clockCountStop
    integer                :: clockCountRate, clockCountMax
    !------------------------------------------------------------------------------

    ! Initialize optional arguments
    persistKDB = .true.
    locExportOptimizationVariables =.false.
    locSkipErrorConvergence =.false.
    locUnitVolume =.false.
    locScalingFactor = fONE
    locScaleHistogram = .false.
    locComputeRawDensity = .false.
    locHistogramScalingFactor = fONE
    locWeightedHistogram = .false.
    locOnlyHistogram = .false.
    locExactPoint = .false.
    localNOptimizationLoops = this%nOptimizationLoops
    this%isotropic = .false.
    this%useGlobalSmoothing = .false.

    ! Process optional arguments !
    if ( present( nOptimizationLoops ) ) then 
      localNOptimizationLoops = nOptimizationLoops
    end if 
    if ( present( outputFileName ) ) then 
      this%outputFileName = outputFileName
    else
      this%outputFileName = defaultOutputFileName
    end if
    if ( present( persistentKernelDatabase ) ) then
      persistKDB = persistentKernelDatabase
    end if
    if ( present( exportOptimizationVariables ) ) then
      locExportOptimizationVariables = exportOptimizationVariables
    end if
    if ( present( skipErrorConvergence ) ) then
      locSkipErrorConvergence = skipErrorConvergence
    end if
    if ( present( unitVolume ) ) then 
      locUnitVolume = unitVolume
    end if 
    if ( present( computeRawDensity ) ) then 
      locComputeRawDensity = computeRawDensity
    end if 
    if ( present( scalingFactor ) ) then 
      locScalingFactor = scalingFactor
    end if
    if ( present( histogramScalingFactor ) ) then
      locScaleHistogram = .true. 
      locHistogramScalingFactor = histogramScalingFactor
    end if
    if ( present( onlyHistogram ) ) then
      locOnlyHistogram = onlyHistogram
    end if
    if ( present( weightedHistogram ) ) then
      locWeightedHistogram = weightedHistogram
    end if
    if ( present( exactPoint ) ) then
      locExactPoint = exactPoint
    end if
    if ( present( relativeErrorConvergence ) ) then
      locRelativeErrorConvergence = relativeErrorConvergence
    else
      locRelativeErrorConvergence = defaultRelativeErrorConvergence
    end if
    if ( present( isotropic ) ) then
      locIsotropic = isotropic 
      this%isotropic = isotropic
    end if
    if ( present( useGlobalSmoothing ) ) then
      this%useGlobalSmoothing = useGlobalSmoothing
      ! If global smoothing, force isotropic 
      if ( this%useGlobalSmoothing ) then 
        this%isotropic = .true.
      end if 
    end if
    if ( (locWeightedHistogram).and.(.not.present(weights)) ) then 
      write(*,*) 'Error: weightedHistogram requires weights and were not given. Stop.'
      stop
    end if
    ! Verify weights size for weighted reconstruction
    dataPointsShape = shape(dataPoints)
    if ( locWeightedHistogram ) then
     if ( size(weights).ne.dataPointsShape(1) ) then 
      write(*,*) 'Error: given weights are not the same length than datapoints. Stop.'
      stop
     end if
    end if
    if ( dataPointsShape(1).lt.1 ) then 
     write(*,*) 'Error: data points is empty. Stop.'
     stop
    end if
    if ( this%reportToOutUnit ) then
     write(this%outFileUnit, *  )
     write(this%outFileUnit, '(1X,A)' ) 'Histogram info '
    end if 
    ! Until here, independent of sliced reconstruction

    ! If adapt grid to coords and all domain sizes 
    ! were zero, then define grid size
    if ( this%adaptGridToCoords .and. all(this%domainSize.eq.fZERO) ) then
      maxCoords        = maxval( dataPoints, dim=1 ) 
      minCoords        = minval( dataPoints, dim=1 )
      deltaCoords      = abs( maxCoords - minCoords )
      where ( deltaCoords .ne. fZERO )
        ! For the minimum coordinates, substract half the border fraction
        minSubGridCoords   = minCoords - fONE*this%borderFraction*deltaCoords
        ! For the maximum coordinates, add half the border fraction
        maxSubGridCoords   = maxCoords + fONE*this%borderFraction*deltaCoords
        ! Assign domain size relative to the origin 
        this%domainSize = maxSubGridCoords - this%domainOrigin 
      end where
    end if

    ! Now it needs to process the bin size. 
    ! Only for 1D reconstruction !
    if ( present( histogramBinFormat ) ) then 
      if ( present( binSizeFraction ) ) then
        if ( (binSizeFraction.le.fZERO).or.(binSizeFraction.gt.fONE) ) then
          write(*,*) 'Error: bin size fraction should be between 0 and 1. Stop.' 
          stop
        end if 
        locBinSizeFraction = binSizeFraction 
      else
        locBinSizeFraction = fONE
      end if 
      select case( histogramBinFormat ) 
      ! compute from data based on Scott's rule
      case (0)
        call this%histogram%EstimateBinSizeScott( dataPoints, this%binSize )
        this%binSize = locBinSizeFraction*this%binSize
        call this%UpdateBinSize( this%binSize ) 
      ! compute from data based on Freedman-Diaconis rule
      case (1)
        call this%histogram%EstimateBinSizeFD( dataPoints, this%binSize )
        this%binSize = locBinSizeFraction*this%binSize
        call this%UpdateBinSize( this%binSize ) 
      case default
        ! Verify that they were defined
        if ( .not. this%adaptGridToCoords ) then 
          if ( all(this%binSize.le.fZERO) ) then  
            write(*,*) 'Error: bin sizes were not defined. Stop.'
            stop
          end if
        else
          call this%UpdateBinSize( this%binSize ) 
        end if
      end select
    else
      ! Something to verify that it was set
      ! and so on... 
      if ( .not. this%adaptGridToCoords ) then 
        if ( all(this%binSize.le.fZERO) ) then  
          write(*,*) 'Error: bin sizes were not defined. Stop.'
          stop
        end if
      else
        call this%UpdateBinSize( this%binSize ) 
      end if
    end if 
    

    ! Compute sub grid parameters if grids
    ! are to be adapted to the given coordinates 
    if ( this%adaptGridToCoords ) then
      maxCoords        = maxval( dataPoints, dim=1 ) 
      minCoords        = minval( dataPoints, dim=1 )
      deltaCoords      = abs( maxCoords - minCoords )
      minSubGridCoords = this%domainOrigin
      maxSubGridCoords = this%domainSize + this%domainOrigin
      where ( this%binSize .ne. fZERO )
        ! For the minimum coordinates, substract half the border fraction
        minSubGridCoords   = minCoords - 0.5*this%borderFraction*deltaCoords
        ! For the maximum coordinates, add half the border fraction
        maxSubGridCoords   = maxCoords + 0.5*this%borderFraction*deltaCoords
      end where

      ! Limit these coordinates by domain specs
      do nd =1,3
        if ( this%dimensionMask(nd).eq.0 ) cycle
        minSubGridCoords(nd) = max(minSubGridCoords(nd), this%domainOrigin(nd))
        maxSubGridCoords(nd) = min(maxSubGridCoords(nd), this%domainSize(nd) + this%domainOrigin(nd) )
      end do

      ! Determine sub grid dimensions
      subGridSize          = this%domainSize
      subGridNBins         = this%domainGridSize
      subGridOriginIndexes = 0 ! relative to domain indexation
      subGridLimitIndexes  = 0 ! same as above 
      where( this%binSize .gt. fZERO ) 
        subGridOriginIndexes = int((minSubGridCoords-this%domainOrigin)/this%binSize)     ! subestimate
        subGridLimitIndexes  = ceiling((maxSubGridCoords-this%domainOrigin)/this%binSize) ! overestimate
        subGridNBins         = subGridLimitIndexes - subGridOriginIndexes ! it shall verify at least 1 ?
        subGridSize          = subGridNBins*this%binSize
      end where
      subGridOrigin = subGridOriginIndexes*this%binSize

      ! Some health control and eventually reporting
      if ( any(subGridNBins.gt.this%domainGridSize) ) then
        write(*,*)'Error: Inconsistent size of subgrid, is larger than domain size.'
        stop
      end if
      
      ! It looks that all these functionalities fit better into the histogram 
      ! but in the meantime...
      ! Assign nBins as the size computed for the sub grid
      this%nBins                = subGridNBins
      this%deltaBinsOrigin      = subGridOriginIndexes
      this%histogram%gridSize   = subGridNBins
      this%histogram%gridOrigin = subGridOrigin
      this%histogram%nBins      => this%histogram%gridSize
      this%histogram%origin     => this%histogram%gridOrigin 

      ! Allocate the counting grid
      if ( allocated( this%histogram%counts ) ) deallocate( this%histogram%counts ) 
      allocate(this%histogram%counts(subGridNBins(1),subGridNBins(2),subGridNBins(3)))
      this%histogram%counts = 0

      ! Allocate matrix for density 
      if ( allocated( densityGrid ) ) deallocate( densityGrid )
      allocate( densityGrid(this%nBins(1), this%nBins(2), this%nBins(3)) )

      if ( this%reportToOutUnit ) then
       write(this%outFileUnit, '(3X,A,3(1X,I9))') 'Allocated size    :', this%nBins
       flush( this%outFileUnit ) 
      end if
    end if ! adaptGridToCoords


    ! Compute histogram
    ! Needs resolution of effective mass when sliced reconstruction
    if ( locWeightedHistogram ) then 
      select case (this%histogram%effectiveWeightFormat)
      case (2)
        ! This format is mostly for analysis.
        ! In this case histogram%counts store the number of points and wcounts the weights
        if ( allocated( this%histogram%wcounts ) ) deallocate( this%histogram%wcounts ) 
        allocate(this%histogram%wcounts, mold=this%histogram%counts)
        call this%histogram%ComputeCountsAndWeights( dataPoints, weights, locExactPoint )
      case (3)
        ! This format is mostly for analysis.
        if ( allocated( this%histogram%wcounts ) ) deallocate( this%histogram%wcounts ) 
        allocate(this%histogram%wcounts, mold=this%histogram%counts)
        call this%histogram%ComputeEffectiveCountsAndWeights( dataPoints, weights, locExactPoint )
      case default
        ! This is the standard/default format
        ! Cummulative histogram-like quantities
        call this%histogram%ComputeCountsWeighted( dataPoints, weights, locExactPoint )
      end select
    else
      ! Histogram quantities
      call this%histogram%ComputeCounts( dataPoints, locExactPoint )
    end if

    ! More info about the histogram data
    if ( this%reportToOutUnit ) then
     write(this%outFileUnit, '(3X,A,es18.9e3)') 'Max count         :', maxval(this%histogram%counts)
     write(this%outFileUnit, '(3X,A,es18.9e3)') 'Max raw density   :', maxval(this%histogram%counts)/this%histogram%binVolume
     write(this%outFileUnit, '(3X,A,es18.9e3)') 'Min count         :', minval(this%histogram%counts)
     write(this%outFileUnit, '(3X,A,es18.9e3)') 'Min raw density   :', minval(this%histogram%counts)/this%histogram%binVolume
    end if

    ! If only histogram, leave
    if ( locOnlyHistogram ) then 
      if( locWeightedHistogram ) then
        ! Restore histogram to mass
        this%histogram%counts = this%histogram%counts*this%histogram%effectiveMass
      end if
      if ( locComputeRawDensity ) then 
        this%histogram%counts = this%histogram%counts/this%histogram%binVolume
        if ( allocated( this%histogram%wcounts) )&
          this%histogram%wcounts = this%histogram%wcounts/this%histogram%binVolume
      end if 
      return
    end if 

    ! Verify the count
    activeBins = count(this%histogram%counts/=fZERO) 
    if ( activeBins .eq. 0 ) then 
     ! No bins to compute 
     if ( this%reportToOutUnit ) then
      write(this%outFileUnit, *  )
      write(this%outFileUnit, '(A)' ) 'Warning: GPKDE module  '
      write(this%outFileUnit, '(A)' ) 'NO bins to compute. Check origin coordinates or particles. Leaving ComputeDensity.'
      write(this%outFileUnit, *  )
     end if
     ! Leaving  
     return
    else
     if ( this%reportToOutUnit ) then
       write(this%outFileUnit, '(3X,A,es18.9e3)' ) 'Mean raw density  :',& 
         sum(this%histogram%counts)/this%histogram%binVolume/activeBins
       write(this%outFileUnit, '(3X,A,I9)'       ) 'Active bins       :', activeBins
       write(this%outFileUnit, '(3X,A,I9)'       ) 'NPoints           :', this%histogram%nPoints
       write(this%outFileUnit, '(3X,A,es18.9e3)' ) 'NEffective        :', this%histogram%nEffective
     end if 
    end if 

    ! Distribution basic statistics
    this%meanCoords = sum(dataPoints,dim=1)/dataPointsShape(1)
    this%stdCoords  = fZERO
    do nd=1,3
      if ( this%dimensionMask(nd) .eq. 0 ) cycle
      this%stdCoords(nd) = sqrt( sum((dataPoints(:,nd)-this%meanCoords(nd))**fTWO)/dataPointsShape(1) )
    end do 
    this%stdSigmaScale = product( this%stdCoords, mask=(this%dimensionMask.eq.1))
    this%stdSigmaScale = this%stdSigmaScale**(fONE/fNDim)
    ! Selects hSigmaScale based on nPoints instead of nEffective
    this%hSigmaScale   = this%stdSigmaScale*( fFOUR/((fNDim + fTWO)*this%histogram%nPoints) )**(fONE/(fNDim+fFOUR))
    if( this%stdSigmaScale .eq. fZERO ) then 
     if ( this%reportToOutUnit ) then
      write(this%outFileUnit, *  )
      write(this%outFileUnit, '(A)' ) 'Warning: GPKDE module  '
      write(this%outFileUnit, '(A)' ) 'Standard deviation is zero. Will continue and lets see what happens.'
      write(this%outFileUnit, *  )
     end if
    else
     if ( this%reportToOutUnit ) then
      write(this%outFileUnit, *  )
      write(this%outFileUnit, '(1X,A)' ) 'Data points info'
      write(this%outFileUnit, '(3X,A,3(1X,es18.9e3))'     ) 'Mean coordinates                 :', this%meanCoords
      write(this%outFileUnit, '(3X,A,3(1X,es18.9e3))'     ) 'Std. dev. coordinates            :', this%stdCoords
      write(this%outFileUnit, '(3X,A,1(1X,es18.9e3))'     ) 'Std. sigma scale                 :', this%stdSigmaScale
      write(this%outFileUnit, '(3X,A,1(1X,es18.9e3))'     ) 'Global smoothing scale Silverman :', this%hSigmaScale
     end if
    end if

    ! Assign min roughness based on specified format
    select case(this%minRoughnessFormat)
    ! Assuming a Gaussian distribution
    case(0)
      select case(nDim) 
      case(1)
        this%minLimitRoughness = &
         this%minRelativeRoughness*(this%histogram%maxRawDensity**fTWO)/(this%stdCoords(this%idDim1)**fFOUR) 
      case(2)
        this%minLimitRoughness = &
         this%minRelativeRoughness*(this%histogram%maxRawDensity**fTWO)/&
         ( product(this%stdCoords**2,dim=1,mask=(this%stdCoords.ne.fZERO)) )
      case(3)
        this%minLimitRoughness = &
         this%minRelativeRoughness*(this%histogram%maxRawDensity**fTWO)/&
         ( product(this%stdCoords**(0.75),dim=1,mask=(this%stdCoords.ne.fZERO)) ) ! 3/4=0.75
      end select
    ! From minRelativeRoughness and a given length scale
    case(1)
      if ( this%minRoughnessLengthScaleAsSigma ) this%minRoughnessLengthScale = this%stdSigmaScale
      this%minLimitRoughness = & 
        this%minRelativeRoughness*(this%histogram%maxRawDensity**fTWO)/&
        (this%minRoughnessLengthScale**fFOUR)
    ! Given value
    case(2)
      continue
    ! Zero
    case(3)
      this%minLimitRoughness = fZERO
    end select

    ! Initialize database if not allocated
    if ( this%databaseOptimization ) then
      if ( .not. allocated( this%kernelDatabaseFlat ) ) then 
          call this%InitializeKernelDatabaseFlat( this%minHOverDelta(1), &
                                                  this%maxHOverDelta(1), &
                                                this%deltaHOverDelta(1), &
                                                  this%logKernelDatabase )
      end if
    end if 

    ! Assign distribution statistics as initial smoothing, Silverman (1986)
    if ( this%initialSmoothingSelection .eq. 0 ) then 
      this%initialSmoothing(:) = this%hSigmaScale
      do nd=1,3
        if ( this%dimensionMask(nd) .eq. 0 ) then 
          this%initialSmoothing(nd) = fZERO
        end if
      end do 
    end if 
    if ( this%reportToOutUnit ) then
      write( this%outFileUnit, *  )
      write( this%outFileUnit, '(1X,A)' ) 'Kernels info'
      write( this%outFileUnit, '(3X,A,3(1X,es18.9e3))') 'initialSmoothing   :', this%initialSmoothing
      flush( this%outFileUnit )
    end if

    ! Logging
    if ( this%reportToOutUnit ) then 
      if ( present( outputDataId ) .and. present( particleGroupId ) ) then 
      timeChar=''
      spcChar = ''
      write(timeChar,*)outputDataId
      write(spcChar,*)particleGroupId
      write( this%outFileUnit, * )
      write( this%outFileUnit, '(A,A,A,A)' )' Bandwidth Optimization -- Time: ', trim(adjustl(timeChar)), &
              ' -- Specie: ', trim(adjustl(spcChar))
      else
        write( this%outFileUnit, *     )
        write( this%outFileUnit, '(A)' )'|-----------------------------------------------------------|'
        write( this%outFileUnit, '(A)' )'| Optimization                                              |'
      end if 
    end if 


    ! Handler for sliced reconstruction
    if ( this%slicedReconstruction ) then
      ! If reconstruction is sliced, then nBins in 
      ! the sliced dimension should be compressed. This 
      ! variable is employed for allocating temporary grids
      ! during reconstruction (e.g. roughness)
      ! Histogram nBins remains with the original size.
      this%nBins(this%slicedDimension) = 1

      ! Loop over slices
      do sliceId=1,this%histogram%nBins(this%slicedDimension)
        ! Active bin ids for this slice
        call this%histogram%ComputeActiveBinIdsSliced(this%slicedDimension, sliceId)
        this%computeBinIds => this%histogram%activeBinIds
        this%nComputeBins  = this%histogram%nActiveBins

        if ( this%reportToOutUnit ) then 
         if ( sliceId.eq.1 ) then 
          write( this%outFileUnit, '(A)' )'|-----------------------------------------------------------|'
         end if 
        end if
        if ( this%nComputeBins .eq. 0 ) then 
         ! No bins to compute 
         if ( this%reportToOutUnit ) then 
            write( this%outFileUnit, '(1X,A,I6,A)')& 
              '  Slice ', sliceId, ' without active bins.'
            write( this%outFileUnit, '(A)' )'|-----------------------------------------------------------|'
         end if
         ! Cycle to the next slice
         cycle
        else
         ! Report
         if ( this%reportToOutUnit ) then 
          write( this%outFileUnit, '(1X,A,I6,A,I10,A)')& 
              '  Slice ', sliceId, ' with ', this%nComputeBins, ' active bins.'
         end if
        end if

        ! Density optimization over a pointer to the slice
        ! sliced given as a range to preserve matrix rank.
        select case(this%slicedDimension) 
        case(1)
          slicedDensity => densityGrid(sliceId:sliceId,:,:)
          this%histogramCounts => this%histogram%counts(sliceId:sliceId,:,:)
          if ( allocated( this%histogram%wcounts ) )&
            this%histogramWCounts => this%histogram%wcounts(sliceId:sliceId,:,:)
        case(2)
          slicedDensity => densityGrid(:,sliceId:sliceId,:)
          this%histogramCounts => this%histogram%counts(:,sliceId:sliceId,:)
          if ( allocated( this%histogram%wcounts ) )&
            this%histogramWCounts => this%histogram%wcounts(:,sliceId:sliceId,:)
        case(3)
          slicedDensity => densityGrid(:,:,sliceId:sliceId)
          this%histogramCounts => this%histogram%counts(:,:,sliceId:sliceId)
          if ( allocated( this%histogram%wcounts ) )&
            this%histogramWCounts => this%histogram%wcounts(:,:,sliceId:sliceId)
        end select

        ! Optimization
        call system_clock(clockCountStart, clockCountRate, clockCountMax)
        call this%ComputeDensityOptimization(                              &
                slicedDensity,                                             &
                nOptimizationLoops=localNOptimizationLoops,                &
                exportOptimizationVariables=locExportOptimizationVariables,&
                skipErrorConvergence=locSkipErrorConvergence,              & 
                relativeErrorConvergence=locRelativeErrorConvergence ) 
        call system_clock(clockCountStop, clockCountRate, clockCountMax)
        if ( this%reportToOutUnit ) then 
          elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
          write( this%outFileUnit, '(A)' )'|-----------------------------------------------------------|'
          write(this%outFileUnit, '(1X,A,E15.5,A)')& 
            '  Optimization time : ', elapsedTime, ' seconds'
          write( this%outFileUnit, '(A)' )'|-----------------------------------------------------------|'
        end if 

      end do 
      ! Restore nBins to its former glory
      ! It is employed downstream in function exporting data
      this%nBins = this%histogram%gridSize

    else
      ! Normal, classical reconstruction, where all grids are coincident 
      ! with the histogram grid

      ! Active bins: Only cells with particles
      call this%histogram%ComputeActiveBinIds()
      this%computeBinIds => this%histogram%activeBinIds
      this%nComputeBins  = this%histogram%nActiveBins
      this%histogramCounts => this%histogram%counts
      if ( allocated( this%histogram%wcounts ) )& 
        this%histogramWCounts => this%histogram%wcounts

      ! Optimization
      call system_clock(clockCountStart, clockCountRate, clockCountMax)
      call this%ComputeDensityOptimization(                              &
              densityGrid,                                               &
              nOptimizationLoops=localNOptimizationLoops,                &
              exportOptimizationVariables=locExportOptimizationVariables,&
              skipErrorConvergence=locSkipErrorConvergence,              & 
              relativeErrorConvergence=locRelativeErrorConvergence ) 
      call system_clock(clockCountStop, clockCountRate, clockCountMax)
      if ( this%reportToOutUnit ) then 
        elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
        write( this%outFileUnit, '(A)' )'|-----------------------------------------------------------|'
        write(this%outFileUnit, '(1X,A,E15.5,A)')& 
          '  Optimization time : ', elapsedTime, ' seconds'
        write( this%outFileUnit, '(A)' )'|-----------------------------------------------------------|'
      end if 

    end if ! slicedReconstruction

    ! Release pointers
    if ( associated( slicedDensity ) ) slicedDensity => null()
    if ( associated( this%histogramCounts ) ) this%histogramCounts => null()
    if ( associated( this%histogramWCounts ) ) this%histogramWCounts => null()

    ! drop database ?
    if ( this%databaseOptimization ) then
      if ( .not. persistKDB ) then
        call this%DropKernelDatabase()
      end if
    end if

    ! Point to the module density
    this%densityEstimateGrid => densityGrid

    ! Some corrections to relevant variables before writing to output files 
    if ( locComputeRawDensity ) then 
      ! Histogram as raw density: histogram/binvolume
      this%histogram%counts = this%histogram%counts/this%histogram%binVolume
      if ( allocated( this%histogram%wcounts ) )&
        this%histogram%wcounts = this%histogram%wcounts/this%histogram%binVolume
    end if 
    if ( locUnitVolume ) then  
      ! If unit volume, modify 
      this%densityEstimateGrid = &
      this%densityEstimateGrid*this%histogram%binVolume
    end if
    if ( locScalingFactor .ne. fZERO ) then
      ! Apply scalingFactor to density
      this%densityEstimateGrid = this%densityEstimateGrid*locScalingFactor
    end if
    if ( locScaleHistogram ) then
      ! Apply histogramScalingFactor to histogram
      this%histogram%counts = this%histogram%counts*locHistogramScalingFactor
      if ( allocated( this%histogram%wcounts ) )&
        this%histogram%wcounts = this%histogram%wcounts*locHistogramScalingFactor
    end if 

    ! Generate vector coordinates if requested. 
    if ( present( generateVectorCoordinates ) ) then 
      if ( generateVectorCoordinates ) then
        call this%GenerateVectorCoordinates()
        !do nd=1,3
        !  if ( this%dimensionMask(nd).eq.1 ) then
        !    select case(nd)
        !    case(1)
        !      if ( allocated( this%coordinatesX ) ) deallocate( this%coordinatesX )
        !      allocate( this%coordinatesX(this%nBins(nd)) )
        !      do m = 1, this%nBins(nd)
        !         idbin = m+this%deltaBinsOrigin(nd)
        !         this%coordinatesX(m) = (real(idbin,fp) + 0.5_fp)*this%binSize(nd) + this%domainOrigin(nd)
        !      end do 
        !    case(2)
        !      if ( allocated( this%coordinatesY ) ) deallocate( this%coordinatesY )
        !      allocate( this%coordinatesY(this%nBins(nd)) )
        !      do m = 1, this%nBins(nd)
        !         idbin = m+this%deltaBinsOrigin(nd)
        !         this%coordinatesY(m) = (real(idbin,fp) + 0.5_fp)*this%binSize(nd) + this%domainOrigin(nd)
        !      end do 
        !    case(3)
        !      if ( allocated( this%coordinatesZ ) ) deallocate( this%coordinatesZ )
        !      allocate( this%coordinatesZ(this%nBins(nd)) )
        !      do m = 1, this%nBins(nd)
        !         idbin = m+this%deltaBinsOrigin(nd)
        !         this%coordinatesZ(m) = (real(idbin,fp) + 0.5_fp)*this%binSize(nd) + this%domainOrigin(nd)
        !      end do 
        !    end select 
        !  end if 
        !end do 
      end if
    end if 


    ! Assign histogram density accordingly, for exporting
    ! data on a scale consistent with density
    this%histogramDensity => null()
    if ( this%histogram%isWeighted ) then 
      select case(this%histogram%effectiveWeightFormat) 
      case (0,1)
        this%histogramDensity => this%histogram%counts
      case (2,3)
        this%histogramDensity => this%histogram%wcounts
      end select
    end if     

    ! Write output files !
    locOutputDataFormat = 0
    if ( present( outputDataFormat ) ) then
      locOutputDataFormat = outputDataFormat 
    end if 
    locOutputColumnFormat = 0
    if ( present( outputColumnFormat ) ) then
      locOutputColumnFormat = outputColumnFormat 
    end if 

    select case(locOutputDataFormat) 
    case (0)
      ! Text-Plain
      if ( & 
        present( outputFileUnit )  .and. & 
        present( outputDataId )    .and. & 
        present( outputDataIdVal ) .and. & 
        present( particleGroupId ) ) then
        ! Suitable for cases idtime, time, pgroupid
        call this%ExportDensityUnit(outputFileUnit, &
          outputDataId=outputDataId, outputDataIdVal=outputDataIdVal,& 
          particleGroupId=particleGroupId, outputColumnFormat=locOutputColumnFormat )
      else if ( present( outputFileUnit ) .and. present( outputDataId ) .and. present( particleGroupId )) then
        ! Suitable for cases idtime, pgroupid
        call this%ExportDensityUnit(&
          outputFileUnit, outputDataId=outputDataId, &
                    particleGroupId=particleGroupId, &
                outputColumnFormat=locOutputColumnFormat )
      else if ( present( outputFileUnit ) ) then
        call this%ExportDensityUnit( outputFileUnit, &
            outputColumnFormat=locOutputColumnFormat )
      else if ( present( outputFileName ) ) then  
        call this%ExportDensity( outputFileName, & 
        outputColumnFormat=locOutputColumnFormat )
      end if
    case(1)
      ! Binary
      if ( & 
        present( outputFileUnit )  .and. & 
        present( outputDataId )    .and. & 
        present( outputDataIdVal ) .and. & 
        present( particleGroupId ) ) then
        ! Suitable for cases idtime, time, pgroupid
        call this%ExportDensityUnitBinary(outputFileUnit, &
          outputDataId=outputDataId, outputDataIdVal=outputDataIdVal,& 
          particleGroupId=particleGroupId, outputColumnFormat=locOutputColumnFormat )
      else if ( present( outputFileUnit ) .and. present( outputDataId ) .and. present( particleGroupId )) then
        call this%ExportDensityUnitBinary(& 
          outputFileUnit, outputDataId=outputDataId, &
                    particleGroupId=particleGroupId, &
                outputColumnFormat=locOutputColumnFormat )
      else if ( present( outputFileUnit ) ) then
        call this%ExportDensityUnitBinary( outputFileUnit, & 
                  outputColumnFormat=locOutputColumnFormat )
      else if ( present( outputFileName ) ) then  
        call this%ExportDensityBinary( outputFileName, &
              outputColumnFormat=locOutputColumnFormat )
      end if
    end select


    ! Done
    return


  end subroutine prComputeDensity 


  ! Density optimization
  subroutine prComputeDensityOptimization( this, densityEstimateGrid, nOptimizationLoops, &
                                       exportOptimizationVariables, skipErrorConvergence, &
                                                                relativeErrorConvergence  )
  !----------------------------------------------------------------------------------------
  ! Performs the optimization loop 
  ! 
  !   - Section 2.5 in Sole-Mari et al.(2019) 
  !----------------------------------------------------------------------------------------
  ! Specifications 
  !----------------------------------------------------------------------------------------
  implicit none
  ! input
  class( GridProjectedKDEType ), target:: this
  real(fp), dimension(:,:,:), intent(inout) :: densityEstimateGrid
  !real(fp), dimension(:,:,:), allocatable, intent(inout) :: densityEstimateGrid
  integer , intent(in), optional :: nOptimizationLoops
  logical , intent(in), optional :: exportOptimizationVariables
  logical , intent(in), optional :: skipErrorConvergence
  real(fp), intent(in), optional :: relativeErrorConvergence
  ! local
  ! kernels
  type( KernelMultiGaussianType )     :: kernel
  type( KernelMultiGaussianType )     :: kernelSigma
  type( KernelSecondDerivativeXType ) :: kernelSDX
  type( KernelSecondDerivativeYType ) :: kernelSDY
  type( KernelSecondDerivativeZType ) :: kernelSDZ
  ! nloops
  integer :: nOptLoops
  ! Grid cells
  type( GridCellType ), dimension(:), pointer :: activeGridCells => null()
  type( GridCellType ), pointer :: gc => null()
  ! kernelMatrix pointer
  real(fp), dimension(:,:,:), pointer :: kernelMatrix => null()
  real(fp), dimension(:,:,:), allocatable, target :: transposedKernelMatrix
  ! Utils
  integer            :: n, m, nd
  character(len=500) :: varsOutputFileName
  character(len=20)  :: loopId
  logical            :: exportVariables, skipErrorBreak
  logical            :: exportLoopError
  integer            :: errorOutputUnit
  ! Optimization error monitoring 
  real(fp), dimension(:), allocatable :: rawDensity
  real(fp), dimension(:), allocatable :: errorMetricArray
  real(fp), dimension(:), allocatable :: kernelSmoothingScaleOld
  real(fp), dimension(:), allocatable :: densityEstimateArrayOld
  real(fp) :: errorRMSE
  real(fp) :: errorRMSEOld
  real(fp) :: errorALMISEProxy 
  real(fp) :: errorALMISEProxyOld
  real(fp) :: errorALMISECumsum
  real(fp) :: errorALMISECumsumOld
  real(fp) :: errorMetricDensity
  real(fp) :: errorMetricDensityOld
  !real(fp) :: errorMetricDensityCumsum
  real(fp) :: errorMetricSmoothing
  real(fp) :: errorMetricSmoothingOld
  !real(fp) :: errorMetricSmoothingCumsum
  real(fp) :: errorMetricConvergence
  real(fp), dimension(3) :: smoothingCarrier
  ! loop n estimate
  real(fp)  :: nEstimate
  real(fp)  :: kernelSigmaScale
  real(fp)  :: kernelScale
  real(fp)  :: density
  ! health
  integer, dimension(3) :: matrixShape
  !----------------------------------------------------------------------------------------

    ! Verify that histogram matrixes were associated
    ! counts
    if ( .not. associated( this%histogramCounts ) ) then 
     if ( allocated( this%histogram%counts ) ) then
      ! Forgive and hope for the best
      this%histogramCounts => this%histogram%counts
     else
      write(*,*)'Error: histogram counts pointer is not associated and counts not allocated. Stop.'
      stop
     end if 
    end if
    if ( this%histogram%isWeighted ) then 
     ! wcounts
     select case(this%histogram%effectiveWeightFormat)
     case(2,3)
      if ( .not. associated( this%histogramWCounts ) ) then 
       if ( allocated( this%histogram%wcounts ) ) then
         ! Forgive and hope for the best
         this%histogramWCounts => this%histogram%wcounts
       else
        write(*,*)'Error: histogram weight counts pointer is not associated and wcounts not allocated. Stop.'
        stop
       end if 
      end if 
     end select
    end if 
    ! Verify shape consistency 
    matrixShape = shape(this%histogramCounts)
    do nd=1,3
      if (this%nBins(nd).ne.matrixShape(nd)) then 
       write(*,*)'Error: histogram counts matrix and nbins is inconsistent. Stop.'
       stop
      end if 
    end do
    if ( this%histogram%isWeighted ) then
     select case(this%histogram%effectiveWeightFormat)
     case(2,3) 
      matrixShape = shape(this%histogramWCounts)
      do nd=1,3
        if (this%nBins(nd).ne.matrixShape(nd)) then 
         write(*,*)'Error: histogram wcounts matrix and nbins is inconsistent. Stop.'
         stop
        end if 
      end do 
     end select
    end if

    ! Initialize vars
    exportVariables  = .false.
    skipErrorBreak   = .false.
    exportLoopError  = .false.
    errorOutputUnit  = 999 
    
    ! Pointers to null
    gc => null()
    kernelMatrix => null()

    ! Allocate arrays according to nComputebins
    call prAllocateArrays( this, this%nComputeBins,&
                           kernelSmoothing,        &
                           kernelSmoothingScale,   &
                           kernelSmoothingShape,   &
                           kernelSigmaSupportScale,&
                           curvatureBandwidth,     &
                           densityEstimateArray,   &
                           nEstimateArray,         &
                           roughnessXXArray,       &
                           roughnessYYArray,       &
                           roughnessZZArray,       &
                           netRoughnessArray,      &
                           activeGridCellsMod)
    if ( allocated(rawDensity) ) deallocate(rawDensity) 
    allocate( rawDensity(this%nComputeBins) )
    if ( allocated(densityEstimateArrayOld) ) deallocate(densityEstimateArrayOld) 
    allocate( densityEstimateArrayOld(this%nComputeBins) )
    if ( allocated(errorMetricArray) ) deallocate(errorMetricArray) 
    allocate( errorMetricArray(this%nComputeBins) )

    ! Initialize allocated arrays with zeroes
    kernelSmoothing = fZERO
    kernelSmoothingScale = fZERO
    kernelSigmaSupportScale = fZERO
    netRoughnessArray = fZERO
    nEstimateArray = fZERO
    densityEstimateArray = fZERO
    curvatureBandwidth = fZERO

    ! Initialize and process arguments
    if ( present( nOptimizationLoops ) ) then 
      nOptLoops = nOptimizationLoops
    else 
      nOptLoops = this%nOptimizationLoops
    end if 
    if ( present( exportOptimizationVariables ) ) then 
      exportVariables = exportOptimizationVariables
    end if 
    ! Error convergence
    if ( present( skipErrorConvergence ) ) then 
      skipErrorBreak = skipErrorConvergence
    end if 
    if ( present( relativeErrorConvergence ) ) then
      if ( relativeErrorConvergence .gt. fZERO ) then 
        errorMetricConvergence =  relativeErrorConvergence
      else 
        errorMetricConvergence = defaultRelativeErrorConvergence
      end if
    else 
      errorMetricConvergence = defaultRelativeErrorConvergence
    end if

    ! Initialize active grid cells
    ! and compute rawDensity
    ! Necessary ?
    rawDensity = fZERO
    !$omp parallel do schedule(dynamic,1) &
    !$omp private( n )                    &
    !$omp private( gc ) 
    do n = 1, this%nComputeBins
      gc => activeGridCellsMod(n)
      call gc%Initialize( this%computeBinIds( :, n ) )
      rawDensity(n) = this%histogramCounts(gc%id(1),gc%id(2),gc%id(3))
      !rawDensity(n) = this%histogram%counts(gc%id(1),gc%id(2),gc%id(3))
    end do
    !$omp end parallel do
    rawDensity = rawDensity/this%histogram%binVolume
    activeGridCells => activeGridCellsMod

    ! Initialize kernels
    call kernel%Initialize(      this%binSize, matrixRange=defaultKernelRange )
    call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )
    call kernelSDX%Initialize( this%binSize, matrixRange=defaultKernelSDRange )
    call kernelSDY%Initialize( this%binSize, matrixRange=defaultKernelSDRange )
    call kernelSDZ%Initialize( this%binSize, matrixRange=defaultKernelSDRange )

    ! Initial smoothing !
    ! 
    ! For a second run, consider something 
    ! to detect previous smoothing values and start 
    ! from there, although the distribution of active/inactive 
    ! bins may change for a second time 
    !
    kernelSmoothing = spread( this%initialSmoothing, 2, this%nComputeBins )
    call prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )
    kernelSigmaSupportScale = defaultInitialSigmaFactor*kernelSmoothingScale
    do nd =1, 3
      if ( this%dimensionMask(nd) .eq. 1 ) then 
        where ( kernelSmoothingScale .gt. fZERO )
          kernelSmoothingShape(nd,:) = kernelSmoothing(nd,:)/kernelSmoothingScale
        end where
      else
        ! No smoothing in compressed dimension 
        kernelSmoothing(nd,:)      = fZERO
        kernelSmoothingShape(nd,:) = fZERO
      end if 
    end do
    kernelSmoothingScaleOld = kernelSmoothingScale

    ! Initialize density grid
    densityEstimateGrid(:,:,:) = fZERO
    densityEstimateArray(:) = fZERO
    !$omp parallel do schedule( dynamic, 1 )  &
    !$omp default( none )                     &
    !$omp shared( this )                      &
    !$omp shared( activeGridCells )           & 
    !$omp shared( kernelSmoothing )           & 
    !$omp private( gc )                       & 
    !$omp private( n )                        & 
    !$omp firstprivate( kernel )              & 
    !$omp reduction( +: densityEstimateGrid )  
    do n = 1, this%nComputeBins
        
      ! Assign gc pointer 
      gc => activeGridCells(n)

      if ( any( kernelSmoothing( :, n ) .lt. fZERO ) ) cycle

      ! Set kernel 
      call this%SetKernel( gc, kernel, kernelSmoothing( :, n ) )

      ! Compute estimate
      densityEstimateGrid(                         &
            gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
            gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
            gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
        ) = densityEstimateGrid(                   &
            gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
            gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
            gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
        ) + this%histogramCounts(                  &
        !) + this%histogram%counts(                 &
            gc%id(1), gc%id(2), gc%id(3) )*gc%kernelMatrix(&
                 gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                 gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                 gc%kernelZMSpan(1):gc%kernelZMSpan(2)  &
        )
      !! Cannot be done here ! Reduction !
      ! Assign into array   
      !densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
    end do
    !$omp end parallel do 
    densityEstimateGrid = densityEstimateGrid/this%histogram%binVolume
    ! Transfer grid density to array
    do n = 1, this%nComputeBins
      gc => activeGridCells(n)
      densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
    end do
    call kernel%ResetMatrix()

    ! Error monitoring
    errorRMSE = sqrt(sum( ((densityEstimateArray - rawDensity)/real(this%histogram%nPoints,fp))**fTWO )/real(this%nComputeBins,fp))

    ! Initialize error metric 
    errorMetricArray = fZERO
    where ( kernelSmoothingScale .ne. fZERO ) 
      errorMetricArray = (nEstimateArray/( (kernelSmoothingScale**fNDim)*(fFOUR*pi)**(0.5*fNDim)) + &
      0.25*netRoughnessArray*kernelSmoothingScale**fFOUR)/(real(this%histogram%nPoints,fp)**fTWO)
    end where
    errorALMISEProxy  = sqrt(sum(errorMetricArray**fTWO)/real(this%nComputeBins,fp))
    errorALMISECumsum = errorALMISEProxy 

    ! Initialize smoothing error trackers
    errorMetricArray = fZERO
    where (kernelSmoothingScaleOld .ne. fZERO ) 
      errorMetricArray = abs(kernelSmoothingScale - & 
           kernelSmoothingScaleOld)/kernelSmoothingScaleOld
    end where
    errorMetricSmoothing = sqrt(sum(errorMetricArray**2)/this%nComputeBins)

    if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, "(a)" )       '|-----------------------------------------------------------|'
      write( this%outFileUnit, "(a,a,a,a)" ) '| Loop |', '  hHatOverDelta  |', '     ALMISE      |', '      RMSE      |'
      write( this%outFileUnit, "(a)" )       '|-----------------------------------------------------------|'
      write( this%outFileUnit, "(I6,3es18.9e3)" ) 0, & 
        sum(kernelSmoothingScale)/this%nComputeBins/this%histogram%binDistance,errorALMISEProxy, errorRMSE
      flush( this%outFileUnit ) 
    end if 

    ! For those special cases where nOptLoops .eq. 0
    if ( nOptLoops .eq. 0 ) then

      ! Export optimization variables 
      if ( (exportVariables) ) then
        write( unit=loopId, fmt=* )0
        write( unit=varsOutputFileName, fmt='(a)' )trim(adjustl(this%outputFileName))//trim(adjustl(loopId))
        call prExportOptimizationVariablesExtended( this, varsOutputFileName, & 
          densityEstimateArray, kernelSmoothing, kernelSmoothingScale,kernelSmoothingShape,  & 
          kernelSigmaSupportScale, &
          curvatureBandwidth, nEstimateArray, roughnessXXArray, &
          roughnessYYArray, roughnessZZArray, netRoughnessArray )
      end if

      ! Final density estimate for weighted histograms
      if ( this%histogram%isWeighted ) then 
        select case(this%histogram%effectiveWeightFormat)
        case(0,1)
          ! The formats where densities are transformed scaling by a uniform weight
          this%histogramCounts = this%histogramCounts*this%histogram%effectiveMass
          !this%histogram%counts = this%histogram%counts*this%histogram%effectiveMass
          densityEstimateGrid = densityEstimateGrid*this%histogram%effectiveMass
        case(2,3)
          ! The formats where histogram stored both 
          ! the count of points and cummulative weights. 
          ! A last reconstruction over the weighted histogram is performed 
          densityEstimateGrid = fZERO
          !$omp parallel do schedule( dynamic, 1 )  &
          !$omp default( none )                     &
          !$omp shared( this )                      &
          !$omp shared( activeGridCells )           & 
          !$omp shared( kernelSmoothing )           & 
          !$omp reduction( +: densityEstimateGrid ) & 
          !$omp firstprivate( kernel )              & 
          !$omp private( n )                        & 
          !$omp private( gc )                        
          do n = 1, this%nComputeBins
            ! Assign pointer 
            gc => activeGridCells(n)
            ! Set kernel
            call this%SetKernel( gc, kernel, kernelSmoothing(:,n) ) 
            ! Compute estimate
            densityEstimateGrid(                         &
                  gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
                  gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
                  gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
              ) = densityEstimateGrid(                   &
                  gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
                  gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
                  gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
              ) + this%histogramWCounts(                &  ! Notice histogram%wcounts !
              !) + this%histogram%wcounts(                &  ! Notice histogram%wcounts !
                  gc%id(1), gc%id(2), gc%id(3) )*gc%kernelMatrix(&
                       gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                       gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                       gc%kernelZMSpan(1):gc%kernelZMSpan(2)  &
              )
          end do
          !$omp end parallel do
          densityEstimateGrid = densityEstimateGrid/this%histogram%binVolume
        end select
      end if 

      ! Clean
      call kernel%Reset()
      call kernelSigma%Reset()
      call kernelSDX%Reset()
      call kernelSDY%Reset()
      call kernelSDZ%Reset()
      ! Deallocate stuff
      kernelMatrix    => null()
      gc              => null()
      activeGridCells => null()
      if( allocated( transposedKernelMatrix ) ) deallocate( transposedKernelMatrix )
      !$omp parallel do
      do n = 1, this%nComputeBins
        call activeGridCellsMod(n)%Reset()
      end do
      !$omp end parallel do
      ! Probably more to be deallocated !
      deallocate( kernelSmoothing        )
      deallocate( kernelSmoothingScale   )
      deallocate( kernelSmoothingShape   )
      deallocate( kernelSigmaSupportScale)
      deallocate( curvatureBandwidth     )
      deallocate( densityEstimateArray   )
      deallocate( nEstimateArray         )
      if (allocated(roughnessXXArray)) deallocate(roughnessXXArray)
      if (allocated(roughnessYYArray)) deallocate(roughnessYYArray)
      if (allocated(roughnessZZArray)) deallocate(roughnessZZArray)
      deallocate( netRoughnessArray      )
      deallocate( activeGridCellsMod     )
      deallocate( rawDensity             ) 
      deallocate( densityEstimateArrayOld) 
      deallocate( errorMetricArray       ) 

      ! Done
      return

    end if

    ! Initialize old error trackers
    errorALMISEProxyOld     = errorALMISEProxy
    errorALMISECumsumOld    = errorALMISECumsum
    errorRMSEOld            = errorRMSE
    errorMetricDensityOld   = fZERO
    errorMetricSmoothingOld = errorMetricSmoothing
    densityEstimateArrayOld = densityEstimateArray
    kernelSmoothingScaleOld = kernelSmoothingScale

    ! Optimization loop !
    do m = 1, nOptLoops

      ! nEstimate
      smoothingCarrier = fZERO
      nEstimate        = fZERO
      kernelSigmaScale = fZERO
      kernelScale      = fZERO
      density          = fZERO
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( densityEstimateGrid )      &
      !$omp shared( nEstimateArray )           &
      !$omp shared( kernelSigmaSupportScale )  &
      !$omp shared( kernelSmoothingScale )     &
      !$omp shared( onePlusNDimQuarter )       &
      !$omp firstprivate( kernelSigma )        &
      !$omp firstprivate( smoothingCarrier )   &
      !$omp firstprivate( nEstimate )          &
      !$omp firstprivate( kernelSigmaScale )   &
      !$omp firstprivate( kernelScale )        &
      !$omp firstprivate( density )            &
      !$omp shared( nDim )                     &
      !$omp private( n )                       &            
      !$omp private( gc )            
      do n = 1, this%nComputeBins
      
        ! Assign gc pointer
        gc => activeGridCells( n )
        density          = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
        kernelSigmaScale = kernelSigmaSupportScale(n) 
        kernelScale      = kernelSmoothingScale(n) 

        ! Set kernel sigma
        smoothingCarrier(:) =  kernelSigmaScale
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )

        ! Compute n estimate
        nEstimate = sum(&
          densityEstimateGrid(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

        ! Scale for the support kernel !
        ! - Eq. 23 in Sole-Mari et al. (2019)
        kernelSigmaScale = nEstimate**(0.5)*kernelScale**onePlusNDimQuarter/&
                       ( ( fFOUR*density )**0.25 )*this%supportDimensionConstant

        ! Bound support size
        if ( this%boundKernels ) then 
         if (kernelSigmaScale.gt.this%maxSigmaGrowth*kernelSigmaSupportScale(n)) &
           kernelSigmaScale=this%maxSigmaGrowth*kernelSigmaSupportScale(n)
         kernelSigmaScale = minval(& 
           (/&
             kernelSigmaScale,                               &
             this%maxKernelSize(this%maxSizeDimId)           & 
           /),dim=1)
         if (kernelSigmaScale.lt.this%minKernelSize(this%minSizeDimId)) &
                 kernelSigmaScale=this%minKernelSize(this%minSizeDimId)
        end if 

        ! Update kernel sigma
        smoothingCarrier(:) = kernelSigmaScale
        call this%SetKernelSigma( gc, kernelSigma, smoothingCarrier )
         
        ! Compute n estimate
        nEstimate = sum(&
          densityEstimateGrid(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

        ! To arrays
        nEstimateArray( n )          = nEstimate
        kernelSigmaSupportScale( n ) = kernelSigmaScale

      end do
      !$omp end parallel do
      call kernelSigma%ResetMatrix()

      ! Curvature bandwidths
      call this%ComputeCurvatureKernelBandwidth( densityEstimateArray, nEstimateArray, &
                                           kernelSmoothingScale, kernelSmoothingShape, & 
                                           kernelSigmaSupportScale, curvatureBandwidth )
      ! Net roughness
      call this%ComputeNetRoughnessEstimate(activeGridCells, curvatureBandwidth, &
                           roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                     netRoughnessArray, kernelSigmaSupportScale, &
                                                 kernelSDX, kernelSDY, kernelSDZ )
      ! Optimal smoothing
      call this%ComputeOptimalSmoothingAndShape( nEstimateArray, netRoughnessArray, & 
                              roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                                                   kernelSmoothing, & 
                                                              kernelSmoothingScale, & 
                                                               kernelSmoothingShape )
      ! Update density
      densityEstimateGrid = fZERO
      densityEstimateArray = fZERO
      !$omp parallel do schedule( dynamic, 1 )  &
      !$omp default( none )                     &
      !$omp shared( this )                      &
      !$omp shared( activeGridCells )           & 
      !$omp shared( kernelSmoothing )           & 
      !$omp reduction( +: densityEstimateGrid ) & 
      !$omp firstprivate( kernel )              & 
      !$omp private( n )                        & 
      !$omp private( gc )                        
      do n = 1, this%nComputeBins

        ! Assign pointer 
        gc => activeGridCells(n)

        ! Set kernel
        call this%SetKernel( gc, kernel, kernelSmoothing(:,n) ) 

        ! Compute estimate
        densityEstimateGrid(                         &
              gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
              gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
              gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
          ) = densityEstimateGrid(                   &
              gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
              gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
              gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
          ) + this%histogramCounts(                  &
          !) + this%histogram%counts(                 &
              gc%id(1), gc%id(2), gc%id(3) )*gc%kernelMatrix(&
                   gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                   gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                   gc%kernelZMSpan(1):gc%kernelZMSpan(2)  &
          )

        !! Cannot be done here ! Reduction !
        ! Assign into array   
        !densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
      end do
      !$omp end parallel do
      densityEstimateGrid = densityEstimateGrid/this%histogram%binVolume
      ! Transfer grid density to array
      do n = 1, this%nComputeBins
        gc => activeGridCells(n)
        densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
      end do
      call kernel%ResetMatrix()

      ! A proxy to error: relative density change
      errorMetricArray = fZERO
      where ( densityEstimateArrayOld .ne. fZERO )
        errorMetricArray = abs(  densityEstimateArray/maxval(densityEstimateArray) -  & 
                                   densityEstimateArrayOld/maxval(densityEstimateArrayOld) & 
                                )/(densityEstimateArrayOld/maxval(densityEstimateArrayOld) )
      end where
      errorMetricDensity = sqrt( sum(errorMetricArray**fTWO)/real(this%nComputeBins,fp) )

      ! A proxy to error: relative smoothing change
      errorMetricArray = fZERO
      where (kernelSmoothingScaleOld .ne. fZERO ) 
        errorMetricArray = abs(kernelSmoothingScale - & 
             kernelSmoothingScaleOld)/kernelSmoothingScaleOld
      end where
      errorMetricSmoothing = sqrt(sum(errorMetricArray**fTWO)/real(this%nComputeBins,fp))

      ! A proxy to error: ALMISE
      errorMetricArray = fZERO
      where ( kernelSmoothingScale .ne. fZERO ) 
        errorMetricArray = (nEstimateArray/( (kernelSmoothingScale**fNDim)*(fFOUR*pi)**(0.5*fNDim)) + &
        0.25*netRoughnessArray*kernelSmoothingScale**fFOUR)/(real(this%histogram%nPoints,fp)**fTWO)
      end where
      errorALMISEProxy = sqrt(sum(errorMetricArray**fTWO)/real(this%nComputeBins,fp))
      errorALMISECumsum = errorALMISECumsum + errorALMISEProxy

      ! A proxy to error: RMSE versus histogram density
      errorRMSE = &
        sqrt(sum(((densityEstimateArray - rawDensity)/real(this%histogram%nPoints,fp))**fTWO)/real(this%nComputeBins,fp))

      ! Error analysis:
      if ( .not. skipErrorBreak ) then
        ! ALMISE convergence
        if ( errorALMISEProxyOld.gt.fZERO ) then 
          if ( abs(errorALMISEProxy - errorALMISEProxyOld)/errorALMISEProxyOld .lt. errorMetricConvergence ) then  
            if ( this%reportToOutUnit ) then 
              write( this%outFileUnit, '(A,es13.4e2)' ) '    - ALMISE convergence :',&
                      abs(errorALMISEProxy - errorALMISEProxyOld)/errorALMISEProxyOld
            end if 
            ! Break
            exit
          end if
        end if
        ! Density convergence
        if ( errorMetricDensity .lt. errorMetricConvergence ) then  
          if ( this%reportToOutUnit ) then 
            write( this%outFileUnit, '(A,es13.4e2)' ) '    - Density convergence :', errorMetricDensity
          end if 
          ! Break
          exit
        end if
        ! Smoothing convergence
        if ( ( errorMetricSmoothing .lt. errorMetricConvergence ) ) then 
          if ( this%reportToOutUnit ) then 
            write( this%outFileUnit, '(A,es13.4e2)' ) '    - Bandwidth convergence :', errorMetricSmoothing
          end if
          ! Break
          exit
        end if
        ! Relative change in averaged ALMISE convergence
        if ( ( m.gt.1 ) ) then
          if (&
            abs( errorALMISECumsumOld/real(m-1,fp) - errorALMISECumsum/real(m,fp) )/&
                 errorALMISECumsumOld/real(m-1,fp) .lt. errorMetricConvergence ) then 
            if ( this%reportToOutUnit ) then
              write( this%outFileUnit, '(A,es13.4e2)' ) '    - ALMISE convergence :',&
               abs( errorALMISECumsumOld/real(m-1,fp) - errorALMISECumsum/real(m,fp) )/&
                    errorALMISECumsumOld/real(m-1,fp)
            end if 
            ! Break
            exit
          end if                   
        end if 
      end if

      ! Continue to next loop !
      errorALMISEProxyOld      = errorALMISEProxy
      errorALMISECumsumOld     = errorALMISECumsum
      errorRMSEOld             = errorRMSE
      errorMetricDensityOld    = errorMetricDensity
      densityEstimateArrayOld  = densityEstimateArray
      kernelSmoothingScaleOld  = kernelSmoothingScale

      ! Export optimization variables
      if ( exportVariables ) then
        write( unit=loopId, fmt=* )m
        write( unit=varsOutputFileName, fmt='(a)' )trim(adjustl(this%outputFileName))//trim(adjustl(loopId))
        call prExportOptimizationVariablesExtended( this, varsOutputFileName, & 
                 densityEstimateArray, kernelSmoothing, kernelSmoothingScale, & 
                               kernelSmoothingShape, kernelSigmaSupportScale, &
                        curvatureBandwidth, nEstimateArray, roughnessXXArray, &
                        roughnessYYArray, roughnessZZArray, netRoughnessArray )
      end if

      if ( this%reportToOutUnit ) then 
      write( this%outFileUnit, "(I6,3es18.9e3)" ) m, & 
        sum(kernelSmoothingScale)/this%nComputeBins/this%histogram%binDistance,errorALMISEProxy, errorRMSE
      flush( this%outFileUnit ) 
      end if 

    end do
    ! End optimization loop ! 

    ! Report if max loops
    if ( ((m-1).eq.nOptLoops).and.this%reportToOutUnit ) then 
      write( this%outFileUnit, '(A)' ) '    - Max loops '
    end if

    ! Final density estimate for weighted histograms
    if ( this%histogram%isWeighted ) then 
      select case(this%histogram%effectiveWeightFormat)
      case(0,1)
        ! The formats where densities are transformed scaling by a uniform weight
        this%histogramCounts = this%histogramCounts*this%histogram%effectiveMass
        !this%histogram%counts = this%histogram%counts*this%histogram%effectiveMass
        densityEstimateGrid = densityEstimateGrid*this%histogram%effectiveMass
      case(2,3)
        ! The formats where histogram stored both 
        ! the count of points and cummulative weights. 
        ! A last reconstruction over the weighted histogram is performed 
        densityEstimateGrid = fZERO
        !$omp parallel do schedule( dynamic, 1 )  &
        !$omp default( none )                     &
        !$omp shared( this )                      &
        !$omp shared( activeGridCells )           & 
        !$omp shared( kernelSmoothing )           & 
        !$omp reduction( +: densityEstimateGrid ) & 
        !$omp firstprivate( kernel )              & 
        !$omp private( n )                        & 
        !$omp private( gc )                        
        do n = 1, this%nComputeBins
          ! Assign pointer 
          gc => activeGridCells(n)
          ! Set kernel
          call this%SetKernel( gc, kernel, kernelSmoothing(:,n) ) 
          ! Compute estimate
          densityEstimateGrid(                         &
                gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
                gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
                gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
            ) = densityEstimateGrid(                   &
                gc%kernelXGSpan(1):gc%kernelXGSpan(2), &
                gc%kernelYGSpan(1):gc%kernelYGSpan(2), & 
                gc%kernelZGSpan(1):gc%kernelZGSpan(2)  & 
            ) + this%histogramWcounts(                &  ! Notice histogram%wcounts !
            !) + this%histogram%wcounts(                &  ! Notice histogram%wcounts !
                gc%id(1), gc%id(2), gc%id(3) )*gc%kernelMatrix(&
                     gc%kernelXMSpan(1):gc%kernelXMSpan(2), &
                     gc%kernelYMSpan(1):gc%kernelYMSpan(2), & 
                     gc%kernelZMSpan(1):gc%kernelZMSpan(2)  &
            )
        end do
        !$omp end parallel do
        densityEstimateGrid = densityEstimateGrid/this%histogram%binVolume
      end select
    end if 

    ! Clean
    call kernel%Reset()
    call kernelSigma%Reset()
    call kernelSDX%Reset()
    call kernelSDY%Reset()
    call kernelSDZ%Reset()

    ! Deallocate stuff
    kernelMatrix    => null()
    gc              => null()
    activeGridCells => null()
    if( allocated( transposedKernelMatrix ) ) deallocate( transposedKernelMatrix )
    !$omp parallel do
    do n = 1, this%nComputeBins
      call activeGridCellsMod(n)%Reset()
    end do
    !$omp end parallel do

    ! Probably more to be deallocated !
    deallocate( kernelSmoothing        )
    deallocate( kernelSmoothingScale   )
    deallocate( kernelSmoothingShape   )
    deallocate( kernelSigmaSupportScale)
    deallocate( curvatureBandwidth     )
    deallocate( densityEstimateArray   )
    deallocate( nEstimateArray         )
    if (allocated(roughnessXXArray)) deallocate(roughnessXXArray)
    if (allocated(roughnessYYArray)) deallocate(roughnessYYArray)
    if (allocated(roughnessZZArray)) deallocate(roughnessZZArray)
    deallocate( netRoughnessArray      )
    deallocate( activeGridCellsMod     )
    deallocate( rawDensity             ) 
    deallocate( densityEstimateArrayOld) 
    deallocate( errorMetricArray       ) 

    ! It run
    this%firstRun = .false.


  end subroutine prComputeDensityOptimization


  ! Optimization loop functions !
  subroutine prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )
    !------------------------------------------------------------------------------
    ! 
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType )          :: this
    real(fp), dimension(:,:), intent(in)   :: kernelSmoothing
    real(fp), dimension(:), intent(inout)  :: kernelSmoothingScale
    integer :: nd
    integer, dimension(:), allocatable :: dimCorrected
    !------------------------------------------------------------------------------
    
    allocate( dimCorrected(this%nComputeBins) )

    dimCorrected = nDim
    kernelSmoothingScale = fONE
    do nd = 1, 3
      if ( this%dimensionMask(nd) .eq. 1 ) then
        where( kernelSmoothing(nd,:) .ne. fZERO )  
          kernelSmoothingScale = kernelSmoothingScale*kernelSmoothing(nd,:) 
        elsewhere
          dimCorrected = dimCorrected - 1  
        end where
      end if 
    end do
    where (dimCorrected .ne. 0 )
      kernelSmoothingScale = ( kernelSmoothingScale )**( fONE/fNDim )
    elsewhere
      kernelSmoothingScale = fZERO
    end where

    deallocate( dimCorrected )

    ! Done
    return

  end subroutine prComputeKernelSmoothingScale


  ! To be deprecated
  subroutine prComputeCurvatureKernelBandwidth( this, densityEstimate, nEstimate, &
                                      kernelSmoothingScale, kernelSmoothingShape, &
                                      kernelSigmaSupportScale, curvatureBandwidth )
    !----------------------------------------------------------------------------
    ! Estimate badwidth for the curvature kernel
    !
    !   - Eqs. 25,26 in Sole-Mari et al. (2019)
    !----------------------------------------------------------------------------
    ! Specifications 
    !----------------------------------------------------------------------------
    implicit none
    ! input
    class( GridProjectedKDEType ) :: this
    real(fp), dimension(:),   intent(in)    :: nEstimate
    real(fp), dimension(:),   intent(in)    :: densityEstimate
    real(fp), dimension(:),   intent(in)    :: kernelSmoothingScale
    real(fp), dimension(:,:), intent(in)    :: kernelSmoothingShape
    real(fp), dimension(:),   intent(in)    :: kernelSigmaSupportScale
    ! out 
    real(fp), dimension(:,:), intent(inout) :: curvatureBandwidth
    ! local 
    real(fp), dimension(:),   allocatable   :: nVirtualPowerBeta
    real(fp), dimension(:,:), allocatable   :: shapeTerm
    real(fp), dimension(:)  , allocatable   :: shapeTermSum
    integer, dimension(3)                   :: shapeTermNums = 1
    integer :: n, nActiveBins, nd
    !----------------------------------------------------------------------------

    ! Allocate local arrays
    nActiveBins = size( nEstimate ) 
    allocate( shapeTerm( 3, nActiveBins  ) )
    allocate( shapeTermSum(  nActiveBins ) )
    allocate( nVirtualPowerBeta( nActiveBins ) )

    ! Compute virtual particle cloud size
    nVirtualPowerBeta = fZERO
    where ( densityEstimate .gt. fZERO )  
      nVirtualPowerBeta = ( ( sqrtEightPi*kernelSigmaSupportScale )**nDim*&
          nEstimate**fTWO/densityEstimate )**this%betaDimensionConstant
    end where

    ! Compute shape dependent terms
    curvatureBandwidth = fZERO
    shapeTerm = fZERO
    do nd = 1, 3
      if ( this%dimensionMask(nd) .eq. 1 ) then 
        shapeTermNums     = 1
        shapeTermNums(nd) = 5
        ! Compute sum for shape term
        shapeTermSum = fZERO
        do n =1,3
          if ( this%dimensionMask(n) .eq. 1 ) then
            where( kernelSmoothingShape(n,:) .ne. fZERO ) 
              shapeTermSum = shapeTermSum + shapeTermNums(n)/( kernelSmoothingShape(n,:)**2 ) 
            end where
          end if
        end do 
        where( kernelSmoothingShape(nd,:) .ne. fZERO ) 
          shapeTerm( nd, : ) = (                                           &
            ( oneOverNDimPlusFour/( kernelSmoothingShape( nd, : )**fFOUR ) )*&
                (                                                          &
                    shapeTermSum                                           &
                )                                                          &
            )**( minusOneOverNDimPlusSix )
        end where

        curvatureBandwidth( nd, : ) = &
          this%alphaDimensionConstant*nVirtualPowerBeta*shapeTerm( nd, : )*kernelSmoothingScale

        ! Bound kernel sizes
        if ( this%boundKernels ) then 
          where( curvatureBandwidth(nd,:).gt.this%maxKernelSDSize(nd) ) 
            curvatureBandwidth(nd,:) = this%maxKernelSDSize(nd)
          end where
          where( curvatureBandwidth(nd,:).lt.this%minKernelSDSize(nd) ) 
            curvatureBandwidth(nd,:) = this%minKernelSDSize(nd)
          end where
        end if

      end if
    end do 

    deallocate( shapeTerm )
    deallocate( shapeTermSum )
    deallocate( nVirtualPowerBeta )

    ! Done
    return

  end subroutine prComputeCurvatureKernelBandwidth



  subroutine prComputeCurvatureBandwidth( this, densityEstimate, nEstimate, &
                                kernelSmoothingScale, kernelSmoothingShape, &
                                kernelSigmaSupportScale, curvatureBandwidth )
    !----------------------------------------------------------------------------
    ! Estimate badwidth for the curvature kernel
    !
    !   - Eqs. 25,26 in Sole-Mari et al. (2019)
    !----------------------------------------------------------------------------
    ! Specifications 
    !----------------------------------------------------------------------------
    implicit none
    ! input
    class( GridProjectedKDEType ) :: this
    real(fp), dimension(:),   intent(in)    :: nEstimate
    real(fp), dimension(:),   intent(in)    :: densityEstimate
    real(fp), dimension(:),   intent(in)    :: kernelSmoothingScale
    real(fp), dimension(:,:), intent(in)    :: kernelSmoothingShape
    real(fp), dimension(:),   intent(in)    :: kernelSigmaSupportScale
    ! out 
    real(fp), dimension(:,:), intent(inout) :: curvatureBandwidth
    ! local 
    real(fp)               :: nVirtualPowerBeta
    real(fp)               :: shapeTerm
    real(fp)               :: shapeTermSum
    real(fp), dimension(3) :: shapeTermNums = fONE
    integer :: n, nd, ndd, did, didd
    !----------------------------------------------------------------------------

    ! Compute curvature bandwidth in parallel 
    curvatureBandwidth = fZERO
    !$omp parallel do schedule( dynamic, 1 ) &
    !$omp default( none )                    &
    !$omp shared( this )                     & 
    !$omp shared( nDim )                     &
    !$omp shared( fNDim )                    &
    !!$omp shared( sqrtEightPi )              & ! is a module parameter
    !$omp shared( oneOverNDimPlusFour )      &
    !$omp shared( minusOneOverNDimPlusSix )  &
    !$omp shared( nEstimate )                &
    !$omp shared( densityEstimate )          &
    !$omp shared( kernelSmoothingScale )     &
    !$omp shared( kernelSmoothingShape )     &
    !$omp shared( kernelSigmaSupportScale )  &
    !$omp shared( curvatureBandwidth )       &
    !$omp private( shapeTerm )               &
    !$omp private( shapeTermSum )            &
    !$omp private( shapeTermNums )           &
    !$omp private( nVirtualPowerBeta )       &
    !$omp private( nd, ndd, did, didd )      &
    !$omp private( n ) 
    do n =1,this%nComputeBins

      ! Something wrong
      if ( densityEstimate(n) .le. fZERO ) cycle

      ! N virtual power beta 
      nVirtualPowerBeta = ( ( sqrtEightPi*kernelSigmaSupportScale(n) )**fNDim*&
          nEstimate(n)**fTWO/densityEstimate(n) )**this%betaDimensionConstant

      ! Loop for dimensions
      do nd =1, nDim
        did = this%dimensions(nd)

        ! Reset shape numbers
        shapeTermNums      = fONE
        shapeTermNums(did) = fFIVE

        ! Sum for shape term, only over active dimensions
        shapeTermSum = fZERO
        do ndd=1, nDim
          didd = this%dimensions(ndd) 
          shapeTermSum = shapeTermSum + shapeTermNums(didd)/( kernelSmoothingShape(didd,n)**fTWO )
        end do
       
        ! The shape term 
        shapeTerm = (                                                    &
          ( oneOverNDimPlusFour/( kernelSmoothingShape(did,n)**fFOUR ) )*&
              (                                                          &
                  shapeTermSum                                           &
              )                                                          &
          )**( minusOneOverNDimPlusSix )
      
        ! The bandwidth for direction did 
        curvatureBandwidth( did, n ) = &
          this%alphaDimensionConstant*nVirtualPowerBeta*shapeTerm*kernelSmoothingScale(n)

        ! Bound kernel sizes
        if ( this%boundKernels ) then 
          if ( curvatureBandwidth(did,n).gt.this%maxKernelSDSize(did) ) curvatureBandwidth(did,n) = this%maxKernelSDSize(did)
          if ( curvatureBandwidth(did,n).lt.this%minKernelSDSize(did) ) curvatureBandwidth(did,n) = this%minKernelSDSize(did)
        end if

      end do

    end do 
    !$omp end parallel do 

    ! Done
    return

  end subroutine prComputeCurvatureBandwidth



  subroutine prComputeOptimalSmoothingAndShape( this, nEstimate, netRoughness, &
                         roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                                              kernelSmoothing, & 
                                                         kernelSmoothingScale, & 
                                                         kernelSmoothingShape  )
    !----------------------------------------------------------------------------
    ! Determines optimal smoothing based on the shape factors obtained 
    ! from roughesses. 
    ! 
    ! Eq. 20b in Sole-Mari et al. (2019)
    !
    !----------------------------------------------------------------------------
    ! Specifications 
    !----------------------------------------------------------------------------
    implicit none
    ! input
    class( GridProjectedKDEType) :: this
    real(fp), dimension(:)  ,         intent(in)    :: nEstimate 
    real(fp), dimension(:)  ,         intent(inout) :: netRoughness 
    real(fp), dimension(:)  , target, intent(in)    :: roughnessXXArray
    real(fp), dimension(:)  , target, intent(in)    :: roughnessYYArray
    real(fp), dimension(:)  , target, intent(in)    :: roughnessZZArray
    real(fp), dimension(:,:),         intent(inout) :: kernelSmoothing
    real(fp), dimension(:)  ,         intent(inout) :: kernelSmoothingScale
    real(fp), dimension(:,:),         intent(inout) :: kernelSmoothingShape
    ! local
    real(fp), dimension(:)  , pointer       :: roughness11Array => null()
    real(fp), dimension(:)  , pointer       :: roughness22Array => null()
    real(fp), dimension(:)  , allocatable   :: sumMainRoughnesses
    integer  :: n, nd
    logical  :: updateScale
    real(fp) :: normShape
    real(fp), parameter :: isoRoughnessFactor = 0.5
    !----------------------------------------------------------------------------

    ! Initialize smoothing scale
    kernelSmoothingScale(:) = fZERO

    ! Announce that something is wrong and leave
    if (maxval(netRoughness).eq.fZERO ) then 
      write(*,*) 'Error: the maximum value of roughness is fZERO, something wrong with discretization or kernel sizes.'
      stop
    end if

    ! If using global smoothing and reached this point 
    ! then netRoughness is safely non-zero. useGlobalSmoothing
    ! forces isotropic kernels. 
    if ( this%useGlobalSmoothing ) then 
     kernelSmoothingScale =&
       ( fNDim/( ( fFOUR*pi )**( 0.5*fNDim )*netRoughness*this%histogram%nEffective ) )**( oneOverNDimPlusFour )
    end if 

    ! Apply a correction to net roughness in case one of the 
    ! directional roughness is strongly dominating above the 
    ! others. Only applies for nDim .ge. 2. and kernels are anisotropic.

    ! If the surface is fully uniform in one dimension, 
    ! let's say, there is a non-zero curvature in x, but zero curvature in y,
    ! then it makes sense that smoothing is controlled by the curvature x. 
    ! However, if until this point computations were performed employing 
    ! anisotropic expressions, then it is likely that net roughness
    ! is close to zero in the scenario that one of the roughnesses is also 
    ! close to zero. This effect is also somehow determined by the alignment 
    ! of the distribution shape with the reference coordinates 

    if ( (nDim .ge. 2).and.(.not.this%isotropic) ) then 
      select case(nDim)
      case(2)
        select case(this%idDim1)
        case(1)
          roughness11Array => roughnessXXArray
        case(2)
          roughness11Array => roughnessYYArray
        case(3)
          roughness11Array => roughnessZZArray
        end select
        select case(this%idDim2)
        case(1)
          roughness22Array => roughnessXXArray
        case(2)
          roughness22Array => roughnessYYArray
        case(3)
          roughness22Array => roughnessZZArray
        end select
        ! It should also verify if this guys are zero, to avoid fpe
        ! At this point it was already discarded that all netRoughness were zero, 
        ! as handled by stop starting this function.
        sumMainRoughnesses = roughness11Array + roughness22Array
        where( sumMainRoughnesses.eq.fZERO )
          sumMainRoughnesses = minval(netRoughness, mask=(netRoughness.gt.fZERO))
        end where
        where(& 
          ((roughness11Array/sumMainRoughnesses).gt.this%isotropicThreshold) .or.& 
          ((roughness22Array/sumMainRoughnesses).gt.this%isotropicThreshold) )
          netRoughness = (1.0_fp-isoRoughnessFactor)*netRoughness +  isoRoughnessFactor*sumMainRoughnesses
        end where
      case(3)
        ! It should also verify if this guys are zero, to avoid fpe
        ! At this point it was already discarded that all netRoughness were zero, 
        ! as handled by stop starting this function.
        sumMainRoughnesses = roughnessXXArray + roughnessYYArray + roughnessZZArray
        where( sumMainRoughnesses.eq.fZERO )
          sumMainRoughnesses = minval(netRoughness, mask=(netRoughness.gt.fZERO))
        end where
        where(& 
          ((roughnessXXArray/sumMainRoughnesses).gt.this%isotropicThreshold) .or.& 
          ((roughnessYYArray/sumMainRoughnesses).gt.this%isotropicThreshold) .or.&
          ((roughnessZZArray/sumMainRoughnesses).gt.this%isotropicThreshold) )
          netRoughness = (1.0_fp-isoRoughnessFactor)*netRoughness +  isoRoughnessFactor*sumMainRoughnesses
        end where
      end select 
    end if 

    ! Compute kernelSmoothingScale
    if ( .not. this%useGlobalSmoothing ) then 
     if ( this%minLimitRoughness .gt. fZERO ) then  
      where ( netRoughness .gt. this%minLimitRoughness )
       kernelSmoothingScale = ( fNDim*nEstimate/( ( fFOUR*pi )**( 0.5*fNDim )*netRoughness ) )**( oneOverNDimPlusFour )
      elsewhere
       ! Estimate a scale based on minLimitRoughness
       kernelSmoothingScale = ( fNDim*nEstimate/( ( fFOUR*pi )**( 0.5*fNDim )*this%minLimitRoughness ) )**( oneOverNDimPlusFour )
      end where
     else
      ! If minLimitRoughness is zero, then compute 
      ! kernelSmoothingScale where possible and for any 
      ! case where netRoughnes is zero then assign the maximum 
      ! computed kernelSmoothingScale.
      where ( netRoughness .gt. this%minLimitRoughness )
       kernelSmoothingScale = ( fNDim*nEstimate/( ( fFOUR*pi )**( 0.5*fNDim )*netRoughness ) )**( oneOverNDimPlusFour )
      end where
      where ( netRoughness .eq. fZERO ) 
       kernelSmoothingScale = maxval(kernelSmoothingScale, mask=(kernelSmoothingScale.gt.fZERO)) 
      end where
     end if
    end if

    ! Shape determination based on roughnesses ! 
    kernelSmoothing(:,:) = fZERO
    kernelSmoothingShape(:,:) = fONE
    if ( .not. this%isotropic ) then 
      select case(nDim)
      case(1)
        continue
      case(2)
        ! At this stage, pointers to dims 1 and 2 were already assigned
        ! If bounded kernels, limit shape factor to something that avoid overcoming limit 
        if ( this%boundKernels ) then
         !$omp parallel do schedule(dynamic,1) & 
         !$omp default( none )                 & 
         !$omp shared( this )                  & 
         !$omp shared( roughness11Array )      & 
         !$omp shared( roughness22Array )      &
         !$omp shared( kernelSmoothingShape )  &
         !$omp shared( kernelSmoothingScale )  &
         !$omp shared( netRoughness )          &
         !$omp private( normShape )            &
         !$omp private( n ) 
         do n=1,this%nComputeBins
           if ( (roughness11Array(n).eq.fZERO).or.(roughness22Array(n).eq.fZERO) ) then 
             kernelSmoothingShape(this%idDim1,n) = fONE
             kernelSmoothingShape(this%idDim2,n) = fONE
             cycle
           end if
           kernelSmoothingShape(this%idDim1,n) = max( (netRoughness(n)/roughness11Array(n))**(fONE/fEIGHT), &
             this%minKernelSize(this%idDim1)/kernelSmoothingScale(n) )
           kernelSmoothingShape(this%idDim2,n) = fONE/kernelSmoothingShape(this%idDim1,n)
           normShape = sqrt(kernelSmoothingShape(this%idDim1,n)*kernelSmoothingShape(this%idDim2,n))
           ! Precaution 
           if (normShape.eq.fZERO) then 
             kernelSmoothingShape(this%idDim1,n) = fONE
             kernelSmoothingShape(this%idDim2,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(this%idDim1,n) = kernelSmoothingShape(this%idDim1,n)/normShape
           kernelSmoothingShape(this%idDim2,n) = kernelSmoothingShape(this%idDim2,n)/normShape
         end do
         !$omp end parallel do
        else
         !$omp parallel do schedule(dynamic,1) & 
         !$omp default( none )                 & 
         !$omp shared( this )                  & 
         !$omp shared( roughness11Array )      & 
         !$omp shared( roughness22Array )      &
         !$omp shared( kernelSmoothingShape )  &
         !$omp shared( kernelSmoothingScale )  &
         !$omp shared( netRoughness )          &
         !$omp private( normShape )            &
         !$omp private( n ) 
         do n=1,this%nComputeBins
           if ( (roughness11Array(n).eq.fZERO).or.(roughness22Array(n).eq.fZERO) ) then 
             kernelSmoothingShape(this%idDim1,n) = fONE
             kernelSmoothingShape(this%idDim2,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(this%idDim1,n) = (netRoughness(n)/roughness11Array(n))**(fONE/fEIGHT) 
           kernelSmoothingShape(this%idDim2,n) = fONE/kernelSmoothingShape(this%idDim1,n)
           normShape = sqrt(kernelSmoothingShape(this%idDim1,n)*kernelSmoothingShape(this%idDim2,n))
           ! Precaution 
           if (normShape.eq.fZERO) then 
             kernelSmoothingShape(this%idDim1,n) = fONE
             kernelSmoothingShape(this%idDim2,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(this%idDim1,n) = kernelSmoothingShape(this%idDim1,n)/normShape
           kernelSmoothingShape(this%idDim2,n) = kernelSmoothingShape(this%idDim2,n)/normShape
         end do
         !$omp end parallel do
        end if
      case(3)
        if ( this%boundKernels ) then
         !$omp parallel do schedule(dynamic,1) & 
         !$omp default( none )                 & 
         !$omp shared( this )                  &
         !$omp shared( roughnessXXArray )      & 
         !$omp shared( roughnessYYArray )      &
         !$omp shared( roughnessZZArray )      &
         !$omp shared( kernelSmoothingShape )  &
         !$omp shared( kernelSmoothingScale )  &
         !$omp shared( netRoughness )          &
         !$omp private( normShape )            &
         !$omp private( n ) 
         do n=1,this%nComputeBins
           if (&
             (roughnessXXArray(n).eq.fZERO).or.&
             (roughnessYYArray(n).eq.fZERO).or.&
             (roughnessZZArray(n).eq.fZERO) ) then 
             kernelSmoothingShape(:,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(1,n) = max( (netRoughness(n)/roughnessXXArray(n))**(fONE/12.0_fp), &
                                                 this%minKernelSize(1)/kernelSmoothingScale(n) )
           kernelSmoothingShape(2,n) = max( (netRoughness(n)/roughnessYYArray(n))**(fONE/12.0_fp), &
                                                 this%minKernelSize(2)/kernelSmoothingScale(n) )
           kernelSmoothingShape(3,n) = fONE/kernelSmoothingShape(1,n)/kernelSmoothingShape(2,n)
           normShape =& 
             (kernelSmoothingShape(1,n)*kernelSmoothingShape(2,n)*kernelSmoothingShape(3,n))**(0.333)
           ! Precaution 
           if (normShape.eq.fZERO) then 
             kernelSmoothingShape(:,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(1,n) = kernelSmoothingShape(1,n)/normShape
           kernelSmoothingShape(2,n) = kernelSmoothingShape(2,n)/normShape
           kernelSmoothingShape(3,n) = kernelSmoothingShape(3,n)/normShape
         end do
         !$omp end parallel do 
        else
         !$omp parallel do schedule(dynamic,1) & 
         !$omp default( none )                 & 
         !$omp shared( this )                  &
         !$omp shared( roughnessXXArray )      & 
         !$omp shared( roughnessYYArray )      &
         !$omp shared( roughnessZZArray )      &
         !$omp shared( kernelSmoothingShape )  &
         !$omp shared( kernelSmoothingScale )  &
         !$omp shared( netRoughness )          &
         !$omp private( normShape )            &
         !$omp private( n ) 
         do n=1,this%nComputeBins
           if (&
             (roughnessXXArray(n).eq.fZERO).or.&
             (roughnessYYArray(n).eq.fZERO).or.&
             (roughnessZZArray(n).eq.fZERO) ) then 
             kernelSmoothingShape(:,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(1,n) = (netRoughness(n)/roughnessXXArray(n))**(fONE/12.0_fp)
           kernelSmoothingShape(2,n) = (netRoughness(n)/roughnessYYArray(n))**(fONE/12.0_fp)
           kernelSmoothingShape(3,n) = fONE/kernelSmoothingShape(1,n)/kernelSmoothingShape(2,n)
           normShape =& 
             (kernelSmoothingShape(1,n)*kernelSmoothingShape(2,n)*kernelSmoothingShape(3,n))**(0.333)
           ! Precaution 
           if (normShape.eq.fZERO) then 
             kernelSmoothingShape(:,n) = fONE
             cycle
           end if 
           kernelSmoothingShape(1,n) = kernelSmoothingShape(1,n)/normShape
           kernelSmoothingShape(2,n) = kernelSmoothingShape(2,n)/normShape
           kernelSmoothingShape(3,n) = kernelSmoothingShape(3,n)/normShape
         end do
         !$omp end parallel do 
        end if
      end select

      ! The product of the shape factors should be one
      if ( .not. all(abs(product(kernelSmoothingShape,dim=1)-fONE).lt.0.99) ) then
        write(*,*) 'Error: the product of kernelSmoothingShape factors is not one. Verify shape calculation.'
        stop
      end if 

    end if ! if .not this%isotropic 

    ! Once shape factors are valid, compute smoothing
    updateScale = .false.
    do nd=1,3
     if( this%dimensionMask(nd).eq.0 ) cycle
     kernelSmoothing(nd,:) = kernelSmoothingShape(nd,:)*kernelSmoothingScale
     if ( this%boundKernels ) then 
      if ( any(kernelSmoothing(nd,:).gt.this%maxKernelSize(nd) ) ) then
       updateScale = .true.
       where( kernelSmoothing(nd,:).gt.this%maxKernelSize(nd) ) 
         kernelSmoothing(nd,:) = this%maxKernelSize(nd)
       end where
      end if
      if ( any(kernelSmoothing(nd,:).lt.this%minKernelSize(nd) ) ) then
       updateScale = .true.
       where( kernelSmoothing(nd,:).lt.this%minKernelSize(nd) ) 
         kernelSmoothing(nd,:) = this%minKernelSize(nd)
       end where
      end if
     end if 
    end do
    if ( updateScale ) then 
      call prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )
    end if 

    ! Done 
    return

  end subroutine prComputeOptimalSmoothingAndShape


  subroutine prComputeKernelDatabaseFlatIndexesLog( this, smoothing, flatDBIndexes, transposeKernel )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType )        :: this
    real(fp), dimension(3), intent(in)   :: smoothing
    integer, dimension(2), intent(inout) :: flatDBIndexes
    logical, intent(inout)               :: transposeKernel
    ! local 
    integer, dimension(3) :: indexes 
    integer :: nd, did
    ! A more robust function would be good 
    !------------------------------------------------------------------------------

    ! Initialize indexes 
    indexes(:) = 1
    transposeKernel = .false.

    ! Compute index value active dimensions
    do nd = 1, nDim
      did = this%dimensions(nd)
      if ( (smoothing( did ).le.fZERO) ) cycle
      indexes(did) = min(&
        max(&
          floor(&
            log( smoothing(did)/this%binSize(did)/this%minHOverDelta(did) )/this%deltaHOverDelta(did)&
          ) + 1, 1 &
        ), &
      this%nDeltaHOverDelta(did)  )
    end do 
    
    ! 1D
    if ( nDim .eq. 1 ) then 
      flatDBIndexes(1) = maxval(indexes)
      flatDBIndexes(2) = 1 
      ! Done 
      return
    end if 

    ! 2D
    ! Will work properly as long nDeltaHOverDelta
    ! has the same value for each axis. 
    ! This is linked to database initialization function.
    if ( nDim .eq. 2 ) then
      if ( indexes(this%idDim1) .lt. indexes(this%idDim2) ) then
        transposeKernel  = .true.
        flatDBIndexes(1) = indexes(this%idDim2)*( indexes(this%idDim2) - 1 )/2 + indexes(this%idDim1)
        flatDBIndexes(2) = 1
      else
        transposeKernel  = .false.
        flatDBIndexes(1) = indexes(this%idDim1)*( indexes(this%idDim1) - 1 )/2 + indexes(this%idDim2)
        flatDBIndexes(2) = 1
      end if
      ! Done 
      return
    end if 

    ! 3D 
    ! Will work properly as long nDeltaHOverDelta
    ! has the same value for each axis. 
    ! This is linked to database initialization function.
    if ( indexes(1) .lt. indexes(2) ) then
      transposeKernel  = .true.
      flatDBIndexes(1) = indexes(2)*( indexes(2) - 1 )/2 + indexes(1)
      flatDBIndexes(2) = indexes(3)
    else
      transposeKernel  = .false.
      flatDBIndexes(1) = indexes(1)*( indexes(1) - 1 )/2 + indexes(2)
      flatDBIndexes(2) = indexes(3)
    end if     

    ! Done
    return


  end subroutine prComputeKernelDatabaseFlatIndexesLog

  ! NEEDS UPDATE
  ! Kernel Database indexes, Flat
  subroutine prComputeKernelDatabaseFlatIndexesLinear( this, smoothing, flatDBIndexes, transposeKernel )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType )        :: this
    real(fp), dimension(3), intent(in)   :: smoothing
    integer, dimension(2), intent(inout) :: flatDBIndexes
    logical, intent(inout)               :: transposeKernel
    ! local 
    integer, dimension(3) :: indexes 
    integer :: nd 
    !------------------------------------------------------------------------------
    
    indexes = 0

    do nd = 1, nDim
      indexes(nd) = min(&
        max(&
          floor(&
            (smoothing(nd)/this%binSize(nd) - this%minHOverDelta(nd))/this%deltaHOverDelta(nd)&
          ) + 1, 1 &
        ), &
      this%nDeltaHOverDelta(nd)  )
    end do 
    

    ! 1D
    if ( nDim .eq. 1 ) then 
      flatDBIndexes(1) = maxval(indexes)
      flatDBIndexes(2) = 1 
      ! Done 
      return
    end if 

    ! Will work properly as long nDeltaHOverDelta
    ! has the same value for each axis. 
    ! This is linked to database initialization function.
    if ( indexes(1) < indexes(2) ) then
      transposeKernel  = .true.
      flatDBIndexes(1) = indexes(2)*( indexes(2) - 1 )/2 + indexes(1)
      flatDBIndexes(2) = indexes(3)
    else
      transposeKernel  = .false.
      flatDBIndexes(1) = indexes(1)*( indexes(1) - 1 )/2 + indexes(2)
      flatDBIndexes(2) = indexes(3)
    end if     

    ! Done
    return

  end subroutine prComputeKernelDatabaseFlatIndexesLinear
  
  ! TO BE DEPRECATED 
  ! NEEDS UPDATE, 
  ! Kernel Database indexes, 3D
  function prComputeKernelDatabaseIndexesLinear( this, smoothing ) result(indexes)
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType )      :: this
    real(fp), dimension(3), intent(in) :: smoothing
    integer, dimension(3) :: indexes 
    integer :: nd 
    !------------------------------------------------------------------------------

    indexes = 0

    do nd = 1, nDim
      indexes(nd) = min(&
        max(&
          floor(&
            (smoothing(nd)/this%binSize(nd) - this%minHOverDelta(nd))/this%deltaHOverDelta(nd)&
          ) + 1, 1 &
        ), &
      this%nDeltaHOverDelta(nd)  )
    end do 

    ! Done
    return

  end function prComputeKernelDatabaseIndexesLinear

  ! TO BE DEPRECATED 
  ! NEEDS UPDATE
  function prComputeKernelDatabaseIndexesLog( this, smoothing ) result(indexes)
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType )      :: this
    real(fp), dimension(3), intent(in) :: smoothing
    integer, dimension(3) :: indexes 
    integer :: nd 
    !------------------------------------------------------------------------------

    ! Initialize indexes 
    indexes = 1

    ! Compute indexes where smoothing > fZERO
    do nd = 1, 3
      if ( smoothing( nd ) .le. fZERO ) cycle
      indexes(nd) = min(&
        max(&
          floor(&
            log( smoothing(nd)/this%binSize(nd)/this%minHOverDelta(nd) )/this%deltaHOverDelta(nd)&
          ) + 1, 1 &
        ), &
      this%nDeltaHOverDelta(nd)  )
    end do 

    ! Done
    return

  end function prComputeKernelDatabaseIndexesLog


  function prComputeLogDatabaseIndex( this, smoothing, dimId ) result(dbindex)
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ) :: this
    real(fp), intent(in)          :: smoothing
    integer , intent(in)          :: dimId
    integer                       :: dbindex
    !------------------------------------------------------------------------------

    ! Initialize index
    dbindex = 1

    ! Compute index for smoothing .gt. fZERO 
    if ( smoothing .le. fZERO ) return
    dbindex = min(&
      max(&
        floor(&
          log( smoothing/this%binSize(dimId)/this%minHOverDelta(dimId) )/this%deltaHOverDelta(dimId)&
        ) + 1, 1 &
      ), &
    this%nDeltaHOverDelta(dimId)  )

    ! Done
    return

  end function prComputeLogDatabaseIndex


  subroutine prSetKernelFromDatabase( this, gridCell, kernel,  smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target       :: this
    type( GridCellType ), intent(inout)         :: gridCell
    class( KernelType ), target, intent(inout)  :: kernel
    real(fp), dimension(3), intent(in)          :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrix
    call kernel%ResetMatrix()

    ! Compute indexes on kernel database
    call this%ComputeKernelDatabaseFlatIndexes( smoothing,   &
      gridCell%kernelDBFlatIndexes, gridCell%transposeKernel )

    ! Copy kernel from database
    call kernel%CopyFrom(& 
      this%kernelDatabaseFlat( gridCell%kernelDBFlatIndexes(1), gridCell%kernelDBFlatIndexes(2) ) )

    if ( .not. gridCell%transposeKernel ) then
      call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
        gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
        gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan  ) 
    else
      call kernel%ComputeSpansBoundedTranspose( gridCell%id, this%nBins,     &
        gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
        gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan  )
    end if 
    gridCell%kernelMatrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelFromDatabase


  subroutine prSetKernelSigmaFromDatabase( this, gridCell, kernel, smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target        :: this
    type( GridCellType ),  intent(inout)         :: gridCell
    class( KernelType ), target, intent(inout)   :: kernel
    real(fp), dimension(3), intent(in)           :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute indexes on kernel database
    ! transposeKernelSigma will always be false as this kernel is isotropic.
    ! Regardless, index function requires to compute indexes in all 
    ! dimensions to take into account potential differences on cell sizes.
    call this%ComputeKernelDatabaseFlatIndexes( smoothing, &
      gridCell%kernelSigmaDBFlatIndexes, gridCell%transposeKernelSigma ) 

    ! Copy kernel from database
    call kernel%CopyFrom(& 
      this%kernelDatabaseFlat( gridCell%kernelSigmaDBFlatIndexes(1), gridCell%kernelSigmaDBFlatIndexes(2) ) )

    ! Determine spans
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSigmaXGSpan, gridCell%kernelSigmaYGSpan, gridCell%kernelSigmaZGSpan, & 
      gridCell%kernelSigmaXMSpan, gridCell%kernelSigmaYMSpan, gridCell%kernelSigmaZMSpan  )

    gridCell%kernelSigmaMatrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelSigmaFromDatabase


  subroutine prSetKernelSD1DFromDatabase( this, gridCell, kernel, smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target        :: this
    type( GridCellType ),  intent(inout)         :: gridCell
    class( KernelType ), target, intent(inout)   :: kernel
    real(fp), dimension(3), intent(in)           :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute indexes on kernel database
    gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

    ! Copy kernel from database
    call kernel%CopyFrom( this%kernelSDDatabase1( gridCell%kernelSDDBIndexes(this%idDim1) ) )

    ! Determine spans
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
      gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan  )

    gridCell%kernelSD1Matrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelSD1DFromDatabase


  subroutine prSetKernelSD2DFromDatabase( this, gridCell, kernel1, kernel2, smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target      :: this
    type( GridCellType ),  intent(inout)       :: gridCell
    class( KernelType ), target, intent(inout) :: kernel1
    class( KernelType ), target, intent(inout) :: kernel2
    real(fp), dimension(3), intent(in)         :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel1%ResetMatrix()
    call kernel2%ResetMatrix()

    ! Compute indexes on kernel database
    gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

    ! Copy kernel from database
    call kernel1%CopyFrom( this%kernelSDDatabase1( gridCell%kernelSDDBIndexes(this%idDim1) ) )

    ! Determine spans
    call kernel1%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
      gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan  )

    ! Copy kernel from database
    call kernel2%CopyFrom( this%kernelSDDatabase2( gridCell%kernelSDDBIndexes(this%idDim2) ) )

    ! Determine spans
    call kernel2%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
      gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan  )

    gridCell%kernelSD1Matrix => kernel1%matrix
    gridCell%kernelSD2Matrix => kernel2%matrix

    ! Done
    return

  end subroutine prSetKernelSD2DFromDatabase


  !subroutine prSetKernelSD3DFromDatabase( this, gridCell, kernel1, kernel2, kernel3, smoothing )
  !  !---------------------------------------------------------------------
  !  ! 
  !  !---------------------------------------------------------------------
  !  ! Specifications 
  !  !---------------------------------------------------------------------
  !  implicit none
  !  class( GridProjectedKDEType ), target      :: this
  !  type( GridCellType ),  intent(inout)       :: gridCell
  !  class( KernelType ), target, intent(inout) :: kernel1
  !  class( KernelType ), target, intent(inout) :: kernel2
  !  class( KernelType ), target, intent(inout) :: kernel3
  !  real(fp), dimension(3), intent(in)  :: smoothing
  !  !---------------------------------------------------------------------

  !  ! Restart kernel matrices
  !  call kernel1%ResetMatrix()
  !  call kernel2%ResetMatrix()
  !  call kernel3%ResetMatrix()

  !  ! Compute indexes on kernel database
  !  gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

  !  ! Copy kernel from database
  !  call kernel1%CopyFrom( this%kernelSDXDatabase( gridCell%kernelSDDBIndexes(1) ) )

  !  ! Determine spans bounded
  !  call kernel1%ComputeSpansBounded( gridCell%id, this%nBins, &
  !    gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
  !    gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan  )

  !  ! Copy kernel from database
  !  call kernel2%CopyFrom( this%kernelSDYDatabase( gridCell%kernelSDDBIndexes(2) ) )

  !  ! Determine spans bounded
  !  call kernel2%ComputeSpansBounded( gridCell%id, this%nBins, &
  !    gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
  !    gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan  )

  !  ! Copy kernel from database
  !  call kernel3%CopyFrom( this%kernelSDZDatabase( gridCell%kernelSDDBIndexes(3) ) )

  !  ! Determine spans bounded
  !  call kernel3%ComputeSpansBounded( gridCell%id, this%nBins, &
  !    gridCell%kernelSD3XGSpan, gridCell%kernelSD3YGSpan, gridCell%kernelSD3ZGSpan, & 
  !    gridCell%kernelSD3XMSpan, gridCell%kernelSD3YMSpan, gridCell%kernelSD3ZMSpan  )

  !  gridCell%kernelSD1Matrix => kernel1%matrix
  !  gridCell%kernelSD2Matrix => kernel2%matrix
  !  gridCell%kernelSD3Matrix => kernel3%matrix

  !  ! Done
  !  return

  !end subroutine prSetKernelSD3DFromDatabase


  subroutine prSetKernelSDFromDatabase( this, & 
    gridCell, kernel, smoothing, kernelDatabase, dimId )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target             :: this
    type( GridCellType ),               intent(inout) :: gridCell
    class( KernelType ) , target,       intent(inout) :: kernel
    real(fp)            ,               intent(in)    :: smoothing
    class( KernelType ) , dimension(:), intent(in)    :: kernelDatabase
    integer             ,               intent(in)    :: dimId
    ! local 
    integer :: dbIndex
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute indexes on kernel database
    dbIndex = prComputeLogDatabaseIndex( this, smoothing, dimId )

    ! Copy kernel from database
    call kernel%CopyFrom( kernelDatabase( dbIndex ) )

    ! Determine spans bounded
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
      gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan  )

    gridCell%kernelSDMatrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelSDFromDatabase


  subroutine prSetKernelBrute( this, gridCell, kernel,  smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target      :: this
    type( GridCellType ), intent(inout)        :: gridCell
    class( KernelType ), target, intent(inout) :: kernel
    real(fp), dimension(3), intent(in)         :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Setup kernel matrix
    call kernel%SetupMatrix( smoothing )

    ! Determine spans
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
      gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan )

    gridCell%kernelMatrix => kernel%matrix

    ! Done 
    return

  end subroutine prSetKernelBrute


  subroutine prSetKernelSigmaBrute( this, gridCell, kernel, smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target      :: this
    type( GridCellType ), intent(inout)        :: gridCell
    class( KernelType ), target, intent(inout) :: kernel
    real(fp), dimension(3), intent(in)         :: smoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Setup kernel matrix
    call kernel%SetupMatrix( smoothing )

    ! Determine spans
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSigmaXGSpan, gridCell%kernelSigmaYGSpan, gridCell%kernelSigmaZGSpan, & 
      gridCell%kernelSigmaXMSpan, gridCell%kernelSigmaYMSpan, gridCell%kernelSigmaZMSpan  )

    gridCell%kernelSigmaMatrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelSigmaBrute


  subroutine prSetKernelSD1DBrute( this, gridCell, kernel, smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target      :: this
    type( GridCellType ),  intent(inout)       :: gridCell
    class( KernelType ), target, intent(inout) :: kernel
    real(fp), dimension(3), intent(in)         :: smoothing
    !-----------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute matrix
    call kernel%SetupMatrix( smoothing )

    ! Determine spans
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
      gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan  )

    gridCell%kernelSD1Matrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelSD1DBrute


  subroutine prSetKernelSD2DBrute( this, gridCell, kernel1, kernel2, smoothing )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target      :: this
    type( GridCellType ),  intent(inout)       :: gridCell
    class( KernelType ), target, intent(inout) :: kernel1
    class( KernelType ), target, intent(inout) :: kernel2
    real(fp), dimension(3), intent(in)         :: smoothing
    real(fp), dimension(3)                     :: locSmoothing
    !---------------------------------------------------------------------

    ! Restart kernel matrices
    call kernel1%ResetMatrix()
    call kernel2%ResetMatrix()

    ! Fill smoothing, curvature kernel is isotropic
    locSmoothing = fZERO
    locSmoothing(this%idDim1) = smoothing(this%idDim1)
    locSmoothing(this%idDim2) = smoothing(this%idDim1)

    ! Compute matrix
    call kernel1%SetupMatrix( locSmoothing )

    ! Determine spans
    call kernel1%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
      gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan  ) 

    ! Fill smoothing, curvature kernel is isotropic
    locSmoothing = fZERO
    locSmoothing(this%idDim1) = smoothing(this%idDim2)
    locSmoothing(this%idDim2) = smoothing(this%idDim2)

    ! Compute matrix
    call kernel2%SetupMatrix( locSmoothing )

    ! Determine spans
    call kernel2%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
      gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan  ) 

    ! For boundaries
    gridCell%kernelSD1Matrix => kernel1%matrix
    gridCell%kernelSD2Matrix => kernel2%matrix

    ! Done
    return

  end subroutine prSetKernelSD2DBrute


  !subroutine prSetKernelSD3DBrute( this, gridCell, kernel1, kernel2, kernel3, smoothing )
  !  !---------------------------------------------------------------------
  !  ! 
  !  !---------------------------------------------------------------------
  !  ! Specifications 
  !  !---------------------------------------------------------------------
  !  implicit none
  !  class( GridProjectedKDEType ), target      :: this
  !  type( GridCellType ),  intent(inout)       :: gridCell
  !  class( KernelType ), target, intent(inout) :: kernel1
  !  class( KernelType ), target, intent(inout) :: kernel2
  !  class( KernelType ), target, intent(inout) :: kernel3
  !  real(fp), dimension(3), intent(in)  :: smoothing
  !  !-----------------------------------------------------------

  !  ! Restart kernel matrices
  !  call kernel1%ResetMatrix()
  !  call kernel2%ResetMatrix()
  !  call kernel3%ResetMatrix()

  !  ! Compute matrix
  !  call kernel1%SetupMatrix( (/smoothing(1),smoothing(1),smoothing(1)/) )

  !  ! Determine spans
  !  call kernel1%ComputeSpansBounded( gridCell%id, this%nBins, &
  !    gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
  !    gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan  ) 

  !  ! Compute matrix
  !  call kernel2%SetupMatrix( (/smoothing(2),smoothing(2),smoothing(2)/) )

  !  ! Determine spans
  !  call kernel2%ComputeSpansBounded( gridCell%id, this%nBins, &
  !    gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
  !    gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan  ) 

  !  ! Compute matrix
  !  call kernel3%SetupMatrix( (/smoothing(3),smoothing(3),smoothing(3)/) )

  !  ! Determine spans
  !  call kernel3%ComputeSpansBounded( gridCell%id, this%nBins, &
  !    gridCell%kernelSD3XGSpan, gridCell%kernelSD3YGSpan, gridCell%kernelSD3ZGSpan, & 
  !    gridCell%kernelSD3XMSpan, gridCell%kernelSD3YMSpan, gridCell%kernelSD3ZMSpan  ) 

  !  ! For boundaries
  !  gridCell%kernelSD1Matrix => kernel1%matrix
  !  gridCell%kernelSD2Matrix => kernel2%matrix
  !  gridCell%kernelSD3Matrix => kernel3%matrix

  !  ! Done
  !  return

  !end subroutine prSetKernelSD3DBrute


  subroutine prSetKernelSDBrute( this, &
    gridCell, kernel, smoothing, kernelDatabase, dimId )
    !---------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------
    ! Specifications 
    !---------------------------------------------------------------------
    implicit none
    class( GridProjectedKDEType ), target           :: this
    type( GridCellType ),          intent(inout)    :: gridCell
    class( KernelType ) , target, intent(inout)     :: kernel
    real(fp)            ,  intent(in)               :: smoothing
    ! For compatibilty with interface, not used
    class( KernelType ) , dimension(:), intent(in)  :: kernelDatabase
    integer             ,               intent(in)  :: dimId
    !-----------------------------------------------------------

    ! Restart kernel matrices
    call kernel%ResetMatrix()

    ! Compute matrix
    call kernel%SetupMatrix( (/smoothing,smoothing,smoothing/) )

    ! Determine spans
    call kernel%ComputeSpansBounded( gridCell%id, this%nBins, &
      gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
      gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan  ) 

    gridCell%kernelSDMatrix => kernel%matrix

    ! Done
    return

  end subroutine prSetKernelSDBrute


  ! Utils kernel database
  function prGenerateLogSpaceData( initPoint, endPoint, nPoints ) result( output )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    real(fp), intent(in) :: initPoint, endPoint
    integer , intent(in) :: nPoints
    real(fp)             :: deltaExponent, deltaAccumulated
    real(fp)             :: logInit, logEnd, logBase
    real(fp), dimension(:), allocatable :: output
    integer :: n
    !-----------------------------------------

    allocate( output( nPoints ) )

    logBase = 10.0_fp
    deltaAccumulated = fZERO
    ! Some sanity to verify init smaller than end

    logInit        = log10( initPoint )
    logEnd         = log10( endPoint  )
    deltaExponent  = ( logEnd - logInit )/( nPoints - 1 )  

    do n = 1, nPoints
      output(n) = logBase**( logInit + deltaAccumulated )
      deltaAccumulated = deltaAccumulated + deltaExponent
    end do
   
    ! Done 
    return

  end function prGenerateLogSpaceData

  ! Utils output files !

  subroutine prExportDensityUnit( this, outputUnit, outputDataId, outputDataIdVal, & 
                                               particleGroupId, outputColumnFormat )
    !------------------------------------------------------------------------------
    ! Export methods reporting cell indexes with respect to domain grid 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class(GridProjectedKDEType), target :: this
    integer, intent(in) :: outputUnit
    integer, optional, intent(in) :: outputDataId
    real(fp),optional, intent(in) :: outputDataIdVal
    integer, optional, intent(in) :: particleGroupId
    integer, optional, intent(in) :: outputColumnFormat
    integer :: ix, iy, iz
    integer :: dataId
    integer :: columnFormat
    integer :: idbinx, idbiny, idbinz
    real(fp), dimension(:,:,:), pointer :: histogramData => null()
    !------------------------------------------------------------------------------

    columnFormat = 0
    if ( present( outputColumnFormat ) ) then 
      select case(outputColumnFormat)
      case(1)
        ! bin ids, cell coordinates and density data
        columnFormat = outputColumnFormat
      case(2)
        ! cell coordinates and density data
        columnFormat = outputColumnFormat
      case default
        ! bin ids and density data
        columnFormat = 0
      end select
    end if  

    ! If the pointer is associated, it meas it was a 
    ! weighted histogram. This is mostly for mpath, 
    ! in order to export histogram data consistent 
    ! with density units, scaled.
    histogramData => null()
    if ( this%histogram%isWeighted ) then 
      if  ( associated( this%histogramDensity ) ) then 
        histogramData => this%histogramDensity
      else
        histogramData => this%histogram%counts
      end if 
    else
      histogramData => this%histogram%counts
    end if 

    ! idtime, time, pgroup
    if ( present( outputDataId ) .and. present( outputDataIdVal ) & 
                               .and. present( particleGroupId ) ) then
      ! Following column-major nesting
      select case(columnFormat)
      case(0)
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              write(outputUnit,"(I8,es18.9e3,4I8,2es18.9e3)") &
              outputDataId, outputDataIdVal, particleGroupId, &
              ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      case(1)
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit,"(I8,es18.9e3,4I8,3es18.9e3,2es18.9e3)")&
                                           outputDataId, outputDataIdVal, particleGroupId, &
                                                                   idbinx, idbiny, idbinz, & 
                        (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                        (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                        (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
                 this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      case(2)
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit,"(I8,es18.9e3,I8,3es18.9e3,2es18.9e3)")&
                                           outputDataId, outputDataIdVal, particleGroupId, &
                        (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                        (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                        (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
                 this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      end select
    else if ( present( outputDataId ) .and. present( particleGroupId ) ) then
      ! Following column-major nesting
      select case(columnFormat)
      case(0)
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              write(outputUnit,"(5I8,2es18.9e3)") outputDataId, particleGroupId, &
              ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      case(1)
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit,"(5I8,3es18.9e3,2es18.9e3)") outputDataId, particleGroupId, &
                                                                   idbinx, idbiny, idbinz, & 
                        (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                        (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                        (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
                 this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      case(2)
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit,"(2I8,3es18.9e3,2es18.9e3)") outputDataId, particleGroupId, &
                        (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                        (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                        (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
                 this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      end select
    else if ( present( outputDataId ) .or. present( particleGroupId ) ) then
      if( present(outputDataId) ) then
        dataId = outputDataId
      else
        dataId = particleGroupId
      end if
      select case(columnFormat)
      case(0)
        ! Following column-major nesting
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              write(outputUnit,"(4I8,2es18.9e3)") dataId, &
              ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      case(1)
        ! Following column-major nesting
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit,"(4I8,3es18.9e3,2es18.9e3)") dataId, idbinx, idbiny, idbinz, & 
                         (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                         (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                         (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
                  this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz ) 
            end do
          end do
        end do
      case(2)
        ! Following column-major nesting
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit,"(1I8,3es18.9e3,2es18.9e3)") dataId, &
              (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
              (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
              (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz ) 
            end do
          end do
        end do
      end select
    else
      select case(columnFormat)
      case(0)
        ! Following column-major nesting
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              write(outputUnit,"(3I8,2es18.9e3)") &
              ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz ) 
            end do
          end do
        end do
      case(1)
        ! Following column-major nesting
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit,"(3I8,3es18.9e3,2es18.9e3)") idbinx, idbiny, idbinz, & 
                 (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                 (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                 (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
                 this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz ) 
            end do
          end do
        end do
      case(2)
        ! Following column-major nesting
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit,"(3es18.9e3,2es18.9e3)") & 
                 (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                 (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                 (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
                 this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz ) 
            end do
          end do
        end do
      end select
    end if 


  end subroutine prExportDensityUnit


  subroutine prExportDensityUnitBinary( this, outputUnit, outputDataId, outputDataIdVal, &
                                                     particleGroupId, outputColumnFormat )
    !------------------------------------------------------------------------------
    ! Export methods reporting cell indexes with respect to domain grid 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class(GridProjectedKDEType), target :: this
    integer, intent(in) :: outputUnit
    integer, optional, intent(in) :: outputDataId
    real(fp),optional, intent(in) :: outputDataIdVal
    integer, optional, intent(in) :: particleGroupId
    integer, optional, intent(in) :: outputColumnFormat
    integer :: ix, iy, iz
    integer :: dataId
    integer :: columnFormat
    integer :: idbinx, idbiny, idbinz
    real(fp), dimension(:,:,:), pointer :: histogramData => null()
    !------------------------------------------------------------------------------

    columnFormat = 0
    if ( present( outputColumnFormat ) ) then 
      select case(outputColumnFormat)
      case(1)
        ! bin ids, cell coordinates and density data
        columnFormat = outputColumnFormat
      case(2)
        ! cell coordinates and density data
        columnFormat = outputColumnFormat
      case default
        ! bin ids and density data
        columnFormat = 0
      end select
    end if  


    ! If the pointer is associated, it meas it was a 
    ! weighted histogram. This is mostly for mpath, 
    ! in order to export histogram data consistent 
    ! with density units, scaled.
    histogramData => null()
    if ( this%histogram%isWeighted ) then 
      if  ( associated( this%histogramDensity ) ) then 
        histogramData => this%histogramDensity
      else
        histogramData => this%histogram%counts
      end if 
    else
      histogramData => this%histogram%counts
    end if 

    ! idtime, time, pgroup
    if ( present( outputDataId ) .and. present( outputDataIdVal ) & 
                               .and. present( particleGroupId ) ) then
      ! Following column-major nesting
      select case(columnFormat)
      case(0)
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              write(outputUnit) outputDataId, outputDataIdVal, particleGroupId, &
              ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      case(1)
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit) outputDataId, outputDataIdVal, particleGroupId,  &
                                                         idbinx, idbiny, idbinz, & 
              (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
              (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
              (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      case(2)
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit) outputDataId, outputDataIdVal, particleGroupId,  &
              (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
              (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
              (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      end select
    else if ( present( outputDataId ) .and. present( particleGroupId ) ) then
      ! Following column-major nesting
      select case(columnFormat)
      case(0)
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              write(outputUnit) outputDataId, particleGroupId, &
              ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      case(1)
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit) outputDataId, particleGroupId, &
                                                                   idbinx, idbiny, idbinz, & 
                        (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                        (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                        (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
                 this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      case(2)
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit) outputDataId, particleGroupId, &
                        (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                        (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                        (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
                 this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      end select
    else if ( present( outputDataId ) .or. present( particleGroupId ) ) then
      if( present(outputDataId) ) then
        dataId = outputDataId
      else
        dataId = particleGroupId
      end if
      select case(columnFormat)
      case(0)
        ! Following column-major nesting
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              write(outputUnit) dataId, &
              ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz )
            end do
          end do
        end do
      case(1)
        ! Following column-major nesting
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit) dataId, idbinx, idbiny, idbinz, & 
                         (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                         (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                         (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
                  this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz ) 
            end do
          end do
        end do
      case(2)
        ! Following column-major nesting
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit) dataId, &
              (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
              (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
              (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz ) 
            end do
          end do
        end do
      end select
    else
      select case(columnFormat)
      case(0)
        ! Following column-major nesting
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              write(outputUnit) &
              ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz ) 
            end do
          end do
        end do
      case(1)
        ! Following column-major nesting
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit) idbinx, idbiny, idbinz, & 
                 (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                 (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                 (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
                 this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz ) 
            end do
          end do
        end do
      case(2)
        ! Following column-major nesting
        do iz = 1, this%nBins(3)
          do iy = 1, this%nBins(2)
            do ix = 1, this%nBins(1)
              if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
              idbinx = ix+this%deltaBinsOrigin(1)
              idbiny = iy+this%deltaBinsOrigin(2)
              idbinz = iz+this%deltaBinsOrigin(3)
              write(outputUnit) & 
                 (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                 (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                 (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
                 this%densityEstimateGrid( ix, iy, iz ), histogramData( ix, iy, iz ) 
            end do
          end do
        end do
      end select
    end if 


  end subroutine prExportDensityUnitBinary




  subroutine prExportOptimizationVariablesExtended( this, outputFileName, &
             densityEstimateArray, kernelSmoothing, kernelSmoothingScale, & 
                           kernelSmoothingShape, kernelSigmaSupportScale, &
       curvatureBandwidth, nEstimate, roughnessXXArray, roughnessYYArray, &
                                           roughnessZZArray, netRoughness )
    !----------------------------------------------------------------------
    ! 
    !----------------------------------------------------------------------
    ! Specifications 
    !----------------------------------------------------------------------
    implicit none 
    class(GridProjectedKDEType) :: this
    character(len=500), intent(in) :: outputFileName
    real(fp), dimension(:)  ,intent(in) :: densityEstimateArray
    real(fp), dimension(:,:),intent(in) :: kernelSmoothing 
    real(fp), dimension(:)  ,intent(in) :: kernelSmoothingScale
    real(fp), dimension(:,:),intent(in) :: kernelSmoothingShape
    real(fp), dimension(:)  ,intent(in) :: kernelSigmaSupportScale
    real(fp), dimension(:,:),intent(in) :: curvatureBandwidth
    real(fp), dimension(:)  ,intent(in) :: nEstimate
    real(fp), dimension(:)  ,intent(in) :: netRoughness
    real(fp), dimension(:) ,intent(in), allocatable :: roughnessXXArray   
    real(fp), dimension(:) ,intent(in), allocatable :: roughnessYYArray
    real(fp), dimension(:) ,intent(in), allocatable :: roughnessZZArray
    integer :: ix, iy, iz, n
    integer :: outputUnit = 555
    !---------------------------------------------------------------------

    ! Write the output file name
    ! Add some default
    open( outputUnit, file=outputFileName, status='replace' )

    do n = 1, this%nComputeBins
      ix = this%computeBinIds( 1, n ) + this%deltaBinsOrigin(1) 
      iy = this%computeBinIds( 2, n ) + this%deltaBinsOrigin(2)
      iz = this%computeBinIds( 3, n ) + this%deltaBinsOrigin(3)
      write(outputUnit,"(3I6,17es18.9e3)") ix, iy, iz,& 
        densityEstimateArray( n ),& 
        kernelSmoothing(1,n), kernelSmoothing(2,n), kernelSmoothing(3,n),& 
        kernelSmoothingShape(1,n), kernelSmoothingShape(2,n), kernelSmoothingShape(3,n),& 
        kernelSmoothingScale(n), kernelSigmaSupportScale(n), &
        curvatureBandwidth(1,n), curvatureBandwidth(2,n), curvatureBandwidth(3,n), &
        !nEstimate(n), roughnessXXArray(n), roughnessYYArray(n), &
        !roughnessZZArray(n), netRoughness(n)
        nEstimate(n), roughnessXXArray(n), 0d0, 0d0 , &
         netRoughness(n)
    end do

    ! Finished
    close(outputUnit)


  end subroutine prExportOptimizationVariablesExtended


  subroutine prExportDensity( this, outputFileName, outputColumnFormat )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class(GridProjectedKDEType)   :: this
    character(len=*), intent(in)  :: outputFileName
    integer, intent(in), optional :: outputColumnFormat
    integer :: ix, iy, iz
    integer :: outputUnit = 555
    integer :: columnFormat
    integer :: idbinx, idbiny, idbinz
    !------------------------------------------------------------------------------

    columnFormat = 0
    if ( present( outputColumnFormat ) ) then 
      select case(outputColumnFormat)
      case(1)
        ! bin ids, cell coordinates and density data
        columnFormat = outputColumnFormat
      case(2)
        ! cell coordinates and density data
        columnFormat = outputColumnFormat
      case default
        ! bin ids and density data
        columnFormat = 0
      end select
    end if  

    ! Write the output file name
    ! Add some default
    open( outputUnit, file=outputFileName, status='replace', access='sequential', form='formatted' )

    select case(columnFormat)
    case(0)
      ! Following column-major nesting
      do iz = 1, this%nBins(3)
        do iy = 1, this%nBins(2)
          do ix = 1, this%nBins(1)
            if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
            ! cellids, density, histogram
            write(outputUnit,"(3I8,2es18.9e3)")& 
              ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), this%histogram%counts( ix, iy, iz ) 
          end do
        end do
      end do
    case(1)
      do iz = 1, this%nBins(3)
        do iy = 1, this%nBins(2)
          do ix = 1, this%nBins(1)
            if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
            idbinx = ix+this%deltaBinsOrigin(1)
            idbiny = iy+this%deltaBinsOrigin(2)
            idbinz = iz+this%deltaBinsOrigin(3)
            write(outputUnit,"(3I8,3es18.9e3,2es18.9e3)") idbinx, idbiny, idbinz, & 
               (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
               (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
               (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
               this%densityEstimateGrid( ix, iy, iz ), this%histogram%counts( ix, iy, iz )
          end do
        end do
      end do
    case(2)
      do iz = 1, this%nBins(3)
        do iy = 1, this%nBins(2)
          do ix = 1, this%nBins(1)
            if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
            idbinx = ix+this%deltaBinsOrigin(1)
            idbiny = iy+this%deltaBinsOrigin(2)
            idbinz = iz+this%deltaBinsOrigin(3)
            write(outputUnit,"(3es18.9e3,2es18.9e3)") &
                      (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                      (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                      (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
               this%densityEstimateGrid( ix, iy, iz ), this%histogram%counts( ix, iy, iz )
          end do
        end do
      end do
    end select


    ! Finished
    close(outputUnit)


  end subroutine prExportDensity


  subroutine prExportDensityBinary( this, outputFileName, outputColumnFormat )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class(GridProjectedKDEType)   :: this
    character(len=*), intent(in)  :: outputFileName
    integer, intent(in), optional :: outputColumnFormat
    integer :: ix, iy, iz
    integer :: outputUnit = 555
    integer :: columnFormat
    integer :: idbinx, idbiny, idbinz
    !------------------------------------------------------------------------------

    columnFormat = 0
    if ( present( outputColumnFormat ) ) then 
      select case(outputColumnFormat)
      case(1)
        ! bin ids, cell coordinates and density data
        columnFormat = outputColumnFormat
      case(2)
        ! cell coordinates and density data
        columnFormat = outputColumnFormat
      case default
        ! bin ids and density data
        columnFormat = 0
      end select
    end if  

    ! Write the output file name
    ! Add some default
    open( outputUnit, file=outputFileName, status='replace', access='stream', form='unformatted' )

    select case(columnFormat)
    case(0)
      ! Following column-major nesting
      do iz = 1, this%nBins(3)
        do iy = 1, this%nBins(2)
          do ix = 1, this%nBins(1)
            if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
            ! cellids, density, histogram
            write(outputUnit)& 
              ix+this%deltaBinsOrigin(1), iy+this%deltaBinsOrigin(2), iz+this%deltaBinsOrigin(3), &
              this%densityEstimateGrid( ix, iy, iz ), this%histogram%counts( ix, iy, iz ) 
          end do
        end do
      end do
    case(1)
      do iz = 1, this%nBins(3)
        do iy = 1, this%nBins(2)
          do ix = 1, this%nBins(1)
            if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
            idbinx = ix+this%deltaBinsOrigin(1)
            idbiny = iy+this%deltaBinsOrigin(2)
            idbinz = iz+this%deltaBinsOrigin(3)
            write(outputUnit) idbinx, idbiny, idbinz, & 
               (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
               (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
               (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
               this%densityEstimateGrid( ix, iy, iz ), this%histogram%counts( ix, iy, iz )
          end do
        end do
      end do
    case(2)
      do iz = 1, this%nBins(3)
        do iy = 1, this%nBins(2)
          do ix = 1, this%nBins(1)
            if ( this%densityEstimateGrid( ix, iy, iz ) .le. fZERO ) cycle
            idbinx = ix+this%deltaBinsOrigin(1)
            idbiny = iy+this%deltaBinsOrigin(2)
            idbinz = iz+this%deltaBinsOrigin(3)
            write(outputUnit) &
                      (real(idbinx,fp) + 0.5_fp)*this%binSize(1) + this%domainOrigin(1), & 
                      (real(idbiny,fp) + 0.5_fp)*this%binSize(2) + this%domainOrigin(2), &
                      (real(idbinz,fp) + 0.5_fp)*this%binSize(3) + this%domainOrigin(3), &
               this%densityEstimateGrid( ix, iy, iz ), this%histogram%counts( ix, iy, iz )
          end do
        end do
      end do
    end select


    ! Finished
    close(outputUnit)


  end subroutine prExportDensityBinary



end module GridProjectedKDEModule
