module GridProjectedKDEModule
    !------------------------------------------------------------------------------
    ! 
    ! 
    !------------------------------------------------------------------------------
    use HistogramModule, only : HistogramType
    use KernelMultiGaussianModule, only : KernelMultiGaussianType, &
                                      KernelSecondDerivativeXType, &
                                      KernelSecondDerivativeYType, &
                                      KernelSecondDerivativeZType, &
                                      KernelType
    use GridCellModule, only : GridCellType
    use omp_lib
    implicit none
    !------------------------------------------------------------------------------


    ! Default configuration
    integer, parameter :: defaultKernelRange   = 3
    integer, parameter :: defaultKernelSDRange = 4
    integer, parameter :: defaultNOptLoops     = 2

    logical, parameter :: defaultFlatKernelDatabase   = .true.
    logical, parameter :: defaultDatabaseOptimization = .true.
    logical, parameter :: defaultLogKernelDatabase    = .true.

    doubleprecision, parameter :: defaultMaxHOverLambda   = 30
    doubleprecision, parameter :: defaultMinHOverLambda   = 1
    doubleprecision, parameter :: defaultDeltaHOverLambda = 0.25
    doubleprecision, parameter :: defaultDensityRelativeConvergence = 0.01
    doubleprecision, parameter :: defaultRelaxedDensityRelativeConvergence = 0.05

    logical, parameter ::  defaultBruteOptimization       = .false. 
    logical, parameter ::  defaultAnisotropicSigmaSupport = .false.

    ! Optimization
    doubleprecision :: defaultInitialSmoothingFactor = 2d0
    doubleprecision :: defaultDensityScale           = 1d0
    !doubleprecision :: defaultMinLimitRoughness      = 1d-4
    !doubleprecision :: defaultMaxLimitRoughness      = 1d4
    doubleprecision :: defaultMinLimitRoughness      = 1d-40
    doubleprecision :: defaultMaxLimitRoughness      = 1d40
    doubleprecision :: defaultMaxSmoothingGrowth     = 10d0
    doubleprecision :: defaultMaxKernelShape         = 10d0
    doubleprecision :: defaultMinKernelShape         = 5d-1


    ! Numerical Parameters
    doubleprecision, parameter :: pi           = 4.d0*atan(1.d0)
    doubleprecision, parameter :: sqrtEightPi  = sqrt(8.d0*4.d0*atan(1.d0))


    ! Module parameters defined after initialization
    integer  :: nDim
    integer, dimension(3) :: dimensionMask = (/1,1,1/)


    ! Set default access to private
    private


    ! Grids
    doubleprecision, dimension(:,:,:), allocatable :: nEstimateGrid ! MOVE IT?

    ! Arrays (MOVE THEM?)
    doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothing
    doubleprecision, dimension(:)    , allocatable :: kernelSmoothingScale
    doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothingShape
    doubleprecision, dimension(:,:)  , allocatable :: kernelSigmaSupport
    doubleprecision, dimension(:)    , allocatable :: kernelSigmaSupportScale
    doubleprecision, dimension(:,:)  , allocatable :: curvatureBandwidth
    doubleprecision, dimension(:,:)  , allocatable :: relativeSmoothingChange
    doubleprecision, dimension(:)    , allocatable :: relativeDensityChange
    doubleprecision, dimension(:)    , allocatable :: densityEstimateArray 
    doubleprecision, dimension(:)    , allocatable :: nEstimateArray
    doubleprecision, dimension(:)    , allocatable , target :: roughnessXXArray
    doubleprecision, dimension(:)    , allocatable , target :: roughnessYYArray
    doubleprecision, dimension(:)    , allocatable , target :: roughnessZZArray
    doubleprecision, dimension(:)    , allocatable :: netRoughnessArray
    
    type( GridCellType ), dimension(:), allocatable, target :: activeGridCellsMod
    
    
    ! Main object
    type, public :: GridProjectedKDEType
    
        ! Properties
        type( HistogramType ) :: histogram
        type( KernelMultiGaussianType ), dimension(:,:,:), allocatable :: kernelDatabase
        type( KernelMultiGaussianType ), dimension(:,:)  , allocatable :: kernelDatabaseFlat
        type( KernelSecondDerivativeXType ), dimension(:), allocatable :: kernelSDXDatabase
        type( KernelSecondDerivativeYType ), dimension(:), allocatable :: kernelSDYDatabase
        type( KernelSecondDerivativeZType ), dimension(:), allocatable :: kernelSDZDatabase


        ! For 2d-3d mapping
        integer :: idDim1, idDim2
        class( KernelType ), dimension(:), pointer :: kernelSDDatabase1
        class( KernelType ), dimension(:), pointer :: kernelSDDatabase2

        ! Initialization
        doubleprecision, dimension(3)   :: binSize
        doubleprecision, dimension(3)   :: domainSize
        doubleprecision, dimension(3)   :: domainOrigin
        doubleprecision, dimension(3)   :: initialSmoothing
        integer        , dimension(3)   :: nBins
        integer        , dimension(3)   :: dimensionMask
    
        ! Variables
        ! Consider replacing some for a common grid, 
        ! replacing values when necessary
        doubleprecision, dimension(:)    , allocatable :: densityEstimate
        doubleprecision, dimension(:,:,:), allocatable :: densityEstimateGrid
        doubleprecision, dimension(:,:,:), allocatable :: rawDensityEstimateGrid
        !doubleprecision, dimension(:,:,:), allocatable :: nEstimateGrid
        doubleprecision, dimension(:,:)  , allocatable :: kernelSmoothing
        doubleprecision, dimension(:,:)  , allocatable :: kernelSigmaSupport
        doubleprecision, dimension(:,:)  , allocatable :: curvatureBandwidth
        
        ! Kernel database params 
        doubleprecision, dimension(3) :: deltaHOverLambda
        doubleprecision, dimension(3) :: minHOverLambda
        doubleprecision, dimension(3) :: maxHOverLambda
        integer, dimension(3)         :: nDeltaHOverLambda ! Computed at kernel databases
        logical                       :: logKernelDatabase
        logical                       :: databaseOptimization 
        logical                       :: flatKernelDatabase
        
        ! Optimization
        logical :: bruteOptimization 
        logical :: anisotropicSigmaSupport 
        integer :: nOptimizationLoops
        doubleprecision :: densityRelativeConvergence
        doubleprecision :: minLimitRoughness
        doubleprecision :: maxLimitRoughness
        doubleprecision :: maxSmoothingGrowth
        doubleprecision :: densityScale
        doubleprecision :: minKernelShape
        doubleprecision :: maxKernelShape
        logical :: firstRun = .true. 
        integer, allocatable, dimension(:,:) :: outputBinIds
        
        doubleprecision, dimension(3) :: averageKernelSmoothing = 0d0


        ! Report to outUnit
        logical :: reportToOutUnit = .false.
        integer :: outFileUnit
        character(len=200) :: outFileName


        ! Module constants
        doubleprecision :: supportDimensionConstant
        doubleprecision :: alphaDimensionConstant
        doubleprecision :: betaDimensionConstant
        
        ! Bins to compute
        integer, dimension(:,:), pointer :: computeBinIds
        integer                          :: nComputeBins = 0
        character( len=300 )             :: outputFileName 
        
        ! Interface
        procedure( ComputeIndexes )    , pass, pointer  :: ComputeKernelDatabaseIndexes      => null()
        procedure( ComputeFlatIndexes ), pass, pointer  :: ComputeKernelDatabaseFlatIndexes  => null()
        procedure( ComputeNetRoughness ), pass, pointer :: ComputeNetRoughnessEstimate       => null()
        procedure( SetKernelInterface )  , pass, pointer  :: SetKernel => null()
        procedure( SetKernelInterface )  , pass, pointer  :: SetKernelSigma => null()
        procedure( SetKernelInterface )  , pass, pointer  :: SetKernelSD    => null()
        procedure( SetKernelInterface2D ), pass, pointer  :: SetKernelSD2D  => null()
        procedure( SetKernelInterface3D ), pass, pointer  :: SetKernelSD3D  => null()
            
    contains
    
        ! Procedures
        procedure :: Initialize                      => prInitialize 
        procedure :: Reset                           => prReset 
        procedure :: InitializeModuleConstants       => prInitializeModuleConstants
        procedure :: InitializeNetRoughnessFunction  => prInitializeNetRoughnessFunction
        procedure :: InitializeKernelDatabaseFlat    => prInitializeKernelDatabaseFlat
        procedure :: DropKernelDatabase              => prDropKernelDatabase
        procedure :: ComputeDensity                  => prComputeDensity
        procedure :: ComputeDensityOptimization      => prComputeDensityOptimization
        procedure :: ComputeSupportScale             => prComputeSupportScale
        procedure :: ComputeCurvatureKernelBandwidth => prComputeCurvatureKernelBandwidth
        procedure :: ComputeOptimalSmoothingAndShape => prComputeOptimalSmoothingAndShape
        procedure :: ExportDensity                   => prExportDensity
        procedure :: ExportDensityUnit               => prExportDensityUnit
        procedure :: GenerateLogSpaceData            => prGenerateLogSpaceData
        procedure :: ComputeXYTranspose              => prComputeXYTranspose
    
    end type GridProjectedKDEType
            
            
    ! Interfaces
    abstract interface
    
        ! ComputeIndexes
        function ComputeIndexes( this, smoothing ) result(indexes)
            import GridProjectedKDEType
            implicit none
            class( GridProjectedKDEType )             :: this
            doubleprecision, dimension(3), intent(in) :: smoothing
            integer, dimension(3) :: indexes 
            integer :: nd 
        end function ComputeIndexes
    
    
        ! ComputeFlatIndexes
        subroutine ComputeFlatIndexes( this, smoothing, flatDBIndexes, transposeKernel )
            import GridProjectedKDEType
            implicit none
            class( GridProjectedKDEType )             :: this
            doubleprecision, dimension(3), intent(in) :: smoothing
            integer, dimension(2), intent(inout)      :: flatDBIndexes
            logical, intent(inout)                    :: transposeKernel
            integer, dimension(3) :: indexes 
            integer :: nd
        end subroutine ComputeFlatIndexes
    
    
        ! NetRoughness
        subroutine ComputeNetRoughness( this, activeGridCells, curvatureBandwidth, &
                             roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                            netRoughnessArray, kernelSigmaSupport, &
                                                   kernelSDX, kernelSDY, kernelSDZ ) 
            import GridProjectedKDEType
            import GridCellType
            import KernelSecondDerivativeXType
            import KernelSecondDerivativeYType
            import KernelSecondDerivativeZType
            implicit none 
            class( GridProjectedKDEType ), target :: this
            type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
            doubleprecision, dimension(:,:), intent(in)            :: curvatureBandwidth
            doubleprecision, dimension(:), intent(inout), target   :: roughnessXXArray
            doubleprecision, dimension(:), intent(inout), target   :: roughnessYYArray
            doubleprecision, dimension(:), intent(inout), target   :: roughnessZZArray
            doubleprecision, dimension(:), intent(inout)           :: netRoughnessArray
            doubleprecision, dimension(:,:), intent(in)            :: kernelSigmaSupport
            type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
            type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
            type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
        end subroutine ComputeNetRoughness
    

        ! SetKernelInterface
        subroutine SetKernelInterface( this, gridCell, kernel, smoothing )
            import GridProjectedKDEType
            import GridCellType
            import KernelType
            implicit none 
            class( GridProjectedKDEType ), target                  :: this
            type( GridCellType ), intent(inout)                    :: gridCell
            class( KernelType ), target, intent(inout)             :: kernel
            doubleprecision, dimension(3), intent(in)              :: smoothing
        end subroutine SetKernelInterface


        ! SetKernelInterface2D
        subroutine SetKernelInterface2D( this, gridCell, kernel1, kernel2, smoothing )
            import GridProjectedKDEType
            import GridCellType
            import KernelType
            implicit none 
            class( GridProjectedKDEType ), target                  :: this
            type( GridCellType ), intent(inout)                    :: gridCell
            class( KernelType ), target, intent(inout)             :: kernel1
            class( KernelType ), target, intent(inout)             :: kernel2
            doubleprecision, dimension(3), intent(in)              :: smoothing
        end subroutine SetKernelInterface2D
    

        ! SetKernelInterface3D
        subroutine SetKernelInterface3D( this, gridCell, kernel1, kernel2, kernel3, smoothing )
            import GridProjectedKDEType
            import GridCellType
            import KernelType
            implicit none 
            class( GridProjectedKDEType ), target                  :: this
            type( GridCellType ), intent(inout)                    :: gridCell
            class( KernelType ), target, intent(inout)             :: kernel1
            class( KernelType ), target, intent(inout)             :: kernel2
            class( KernelType ), target, intent(inout)             :: kernel3
            doubleprecision, dimension(3), intent(in)              :: smoothing
        end subroutine SetKernelInterface3D

    end interface


    ! Module contains
    contains


    ! Subroutines
    subroutine prInitialize( this, domainSize, binSize, initialSmoothing, &
                                databaseOptimization, flatKernelDatabase, &
                                          minHOverLambda, maxHOverLambda, &
                                     deltaHOverLambda, logKernelDatabase, &
                              bruteOptimization, anisotropicSigmaSupport, &
                                        nOptimizationLoops, domainOrigin, &
                                densityRelativeConvergence, densityScale, &
                          maxRoughness, minRoughness, maxSmoothingGrowth, &
                                          minKernelShape, maxKernelShape, & 
                                                              outFileName )
      !----------------------------------------------------------------------------
      !
      !----------------------------------------------------------------------------
      ! Specifications 
      !----------------------------------------------------------------------------
      implicit none
      class( GridProjectedKDEType ) :: this
      ! Reconstruction grid parameters
      doubleprecision, dimension(3), intent(in) :: domainSize
      doubleprecision, dimension(3), intent(in) :: binSize
      doubleprecision, dimension(3), intent(in), optional :: domainOrigin
      doubleprecision, dimension(3), intent(in), optional :: initialSmoothing
      integer, intent(in), optional :: nOptimizationLoops 
      ! Kernel database parameters
      logical, intent(in), optional :: databaseOptimization, flatKernelDatabase
      doubleprecision, intent(in), optional :: minHOverLambda, maxHOverLambda
      doubleprecision, intent(in), optional :: deltaHOverLambda
      doubleprecision, intent(in), optional :: densityRelativeConvergence  
      doubleprecision, intent(in), optional :: densityScale, maxSmoothingGrowth
      doubleprecision, intent(in), optional :: minRoughness, maxRoughness
      doubleprecision, intent(in), optional :: minKernelShape, maxKernelShape
      logical, intent(in), optional :: logKernelDatabase
      ! Brute optimization, no kernel database
      logical, intent(in), optional :: bruteOptimization, anisotropicSigmaSupport
      ! General use, indexes
      integer :: n
      ! Limit roughness
      doubleprecision ::  minHRoughness, maxHRoughness, nDensityScale
      ! The analog to a listUnit, reports
      character(len=200), intent(in), optional :: outFileName
      integer :: isThisFileOpen
      ! Time monitoring
      !----------------------------------------------------------------------------


      ! Enable reporting to outUnit if given 
      if( present( outFileName ) ) then
        isThisFileOpen = -1
        inquire( file=outFileName, number=isThisFileOpen )
        if ( isThisFileOpen .gt. 0 ) then 
          this%reportToOutUnit = .true.
          this%outFileUnit = isThisFileOpen
          this%outFileName = outFileName
          write( this%outFileUnit, * )
          write( this%outFileUnit, '(A)' ) '------------------------------'
          write( this%outFileUnit, '(A)' ) ' GPKDE module is initializing '
        end if
      else if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, * )
        write( this%outFileUnit, '(A)' ) '------------------------------'
        write( this%outFileUnit, '(A)' ) ' GPKDE module is initializing '
      end if 

      ! Stop if all bin sizes are zero
      if ( all( binSize .lt. 0d0 ) ) then 
        write(*,*) 'Error while initializing GPKDE, all binSizes are .lt. 0d0. Stop.'
        stop 
      end if 

      ! Initialize reconstruction grid 
      where( binSize .ne. 0d0 ) 
        this%nBins = ceiling( domainSize/binSize )
      elsewhere
        this%nBins = 1
      end where

      ! Stop if any nBins .lt. 1
      if ( any( this%nBins .lt. 1 ) ) then 
        write(*,*) 'Error while initializing GPKDE, some nBins .lt. 1. Stop.'
        stop 
      end if
      this%binSize    = binSize
      this%domainSize = domainSize

      ! domainOrigin
      if ( present( domainOrigin ) ) then 
        this%domainOrigin = domainOrigin
      else 
        this%domainOrigin = (/0,0,0/)
      end if

      ! Depending on nBins, is the number of dimensions 
      ! of the  GPDKE reconstruction process. If any nBins is 1, 
      ! then that dimension is compressed. e.g. nBins = (10,1,20),
      ! then it is a 2D reconstruction process where dimensions
      ! x and z define the 2D plane. This is not necessarily the 
      ! same for the computation of Histograms, where determination 
      ! of a particle inside the grid is related to the 
      ! binSize. If a given binSize is zero, then histogram computation 
      ! does not consider this dimension. If nBins .eq. 1 and binSize .gt. 0
      ! then dimension is considered as valid, and compared against the
      ! origin.

      ! Initialize module dimensions
      call prInitializeModuleDimensions( this, nDim, dimensionMask ) 

      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, *) '  Given binSize      :', this%binSize
        write( this%outFileUnit, *) '  Given domainSize   :', this%domainSize
        write( this%outFileUnit, *) '  Given domainOrigin :', this%domainOrigin
        write( this%outFileUnit, *) '  Computed nBins     :', this%nBins
        write( this%outFileUnit, *) '  Dimensionality for reconstruction is determined from nBins '
        write( this%outFileUnit, *) '  Will perform reconstruction in ', nDim, ' dimensions.'
      end if  
      
      ! Initialize module constants, uses nDim
      call this%InitializeModuleConstants()

      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(A)' ) ' GPKDE initializing Histogram '
      end if

      ! Initialize histogram
      call this%histogram%Initialize( &
            this%nBins, this%binSize, &
         dimensionMask=dimensionMask, & 
       domainOrigin=this%domainOrigin )

      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, *) '  Histogram determines dimensions to be analyzed based on binSizes '
        write( this%outFileUnit, *) '  Will compute Histogram considering ', this%histogram%nDim, ' dimensions '
      end if  
      
      ! Process optional arguments
      ! Kernel database 
      if ( present( databaseOptimization ) ) then 
        this%databaseOptimization = databaseOptimization
      else 
        this%databaseOptimization = defaultDatabaseOptimization
      end if 
      ! flatKernelDatabase
      if ( present( flatKernelDatabase ) ) then 
        this%flatKernelDatabase = flatKernelDatabase
      else 
        this%flatKernelDatabase = defaultFlatKernelDatabase
      end if
      ! Process kernel database discretization parameters 
      if ( present( maxHOverLambda ) ) then 
        this%maxHOverLambda = maxHOverLambda
      else 
        this%maxHOverLambda = defaultMaxHOverLambda
      end if
      if ( present( minHOverLambda ) ) then 
        this%minHOverLambda = minHOverLambda
      else 
        this%minHOverLambda = defaultMinHOverLambda
      end if
      if ( present( deltaHOverLambda ) ) then 
        this%deltaHOverLambda = deltaHOverLambda
      else 
        this%deltaHOverLambda = defaultDeltaHOverLambda
      end if
      if ( present( densityRelativeConvergence ) ) then 
        this%densityRelativeConvergence = densityRelativeConvergence
      else 
        this%densityRelativeConvergence = defaultDensityRelativeConvergence
      end if
      if ( present( minRoughness ) ) then 
        this%minLimitRoughness = minRoughness
      else 
        this%minLimitRoughness = defaultMinLimitRoughness
      end if
      if ( present( maxRoughness ) ) then 
        this%maxLimitRoughness = maxRoughness
      else 
        this%maxLimitRoughness = defaultMaxLimitRoughness
      end if
      if ( present( logKernelDatabase ) ) then 
        this%logKernelDatabase = logKernelDatabase
      else 
        this%logKernelDatabase = defaultLogKernelDatabase
      end if
      ! bruteOptimization, not used
      if ( present( bruteOptimization ) ) then 
        this%bruteOptimization = bruteOptimization
      else 
        this%bruteOptimization = defaultBruteOptimization       
      end if
      ! Not implemented 
      if ( present( anisotropicSigmaSupport ) ) then 
        this%anisotropicSigmaSupport = anisotropicSigmaSupport
      else 
        this%anisotropicSigmaSupport = defaultAnisotropicSigmaSupport
      end if
      ! nOptimizationLoops
      if ( present( nOptimizationLoops ) ) then 
        this%nOptimizationLoops = nOptimizationLoops
      else 
        this%nOptimizationLoops = defaultNOptLoops
      end if 
      ! Initialize smoothing,
      ! could be a vector for active bins 
      if ( present( initialSmoothing ) ) then
        this%initialSmoothing = initialSmoothing
      else
        ! The initial estimate could be improved
        this%initialSmoothing = defaultInitialSmoothingFactor*this%histogram%binDistance
      end if 
      ! Fix to be consistent with dimensions 
      do n =1,3
        if ( dimensionMask(n) .eq. 0 ) then 
          this%initialSmoothing(n) = 0d0
        end if 
      end do
      ! Smoothing growth
      if ( present( maxSmoothingGrowth ) ) then 
        this%maxSmoothingGrowth = maxSmoothingGrowth
      else 
        this%maxSmoothingGrowth = defaultMaxSmoothingGrowth
      end if
      ! Limit kernel shapes
      if ( present( minKernelShape ) ) then 
        this%minKernelShape = minKernelShape
      else 
        this%minKernelShape = defaultMinKernelShape
      end if
      if ( present( maxKernelShape ) ) then 
        this%maxKernelShape = maxKernelShape
      else 
        this%maxKernelShape = defaultMaxKernelShape
      end if
      ! Density scale
      if ( present( densityScale ) ) then 
        this%densityScale = densityScale
      else 
        this%densityScale = defaultDensityScale
      end if
      
      ! Limit roughnesses based on limit smoothing
      ! Default, nDensityScale = 1d0
      ! Consider only left as user parameter
      nDensityScale = this%densityScale
      maxHRoughness = maxval( this%maxHOverLambda*this%binSize )
      !this%minLimitRoughness = maxval( (/this%minLimitRoughness, &
      !    nDim*nDensityScale/( ( maxHRoughness**(nDim + 4d0) )*(4d0*pi)**(0.5*nDim) ) /) )
      minHRoughness = minval( this%minHOverLambda*this%binSize )
      !this%maxLimitRoughness = minval( (/ nDim*nDensityScale/( ( minHRoughness**(nDim + 4d0) )*(4d0*pi)**(0.5*nDim) ), &
      !    this%maxLimitRoughness /) )


      ! Initialize kernel database (ISSUES in 2D)
      if ( this%databaseOptimization ) then
        if ( this%flatKernelDatabase ) then
          ! Initialize kernel database 
          ! Needs review for 2D
          call this%InitializeKernelDatabaseFlat( this%minHOverLambda(1), &
                                                  this%maxHOverLambda(1), &
                                                this%deltaHOverLambda(1), &
                                                  this%logKernelDatabase  )
        end if
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

      ! Allocate matrixes for density
      if ( allocated( this%densityEstimateGrid ) ) deallocate( this%densityEstimateGrid )
      allocate( this%densityEstimateGrid(this%nBins(1), this%nBins(2), this%nBins(3)) )
      if ( allocated( nEstimateGrid ) ) deallocate( nEstimateGrid )
      allocate( nEstimateGrid(this%nBins(1), this%nBins(2), this%nBins(3)) )

      ! Report intialization
      if ( this%reportToOutUnit ) then 
        write( this%outFileUnit, '(A)' ) ' GPKDE module is initialized  '
        write( this%outFileUnit, '(A)' ) '------------------------------'
        write( this%outFileUnit,  *    ) 
      end if


    end subroutine prInitialize


    subroutine prReset( this )
      !------------------------------------------------------------------------------
      ! Incomplete ! 
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      implicit none
      class( GridProjectedKDEType ) :: this
      !------------------------------------------------------------------------------

      call this%histogram%Reset()
      !call this%kernel%Reset()


      ! MAYBE HERE
      !deallocate( kernelSmoothing )
      !deallocate( kernelSigmaSupport )

      !deallocate( densityEstimateActiveBins )
      !deallocate( nEstimateActiveBins )

      !deallocate( densityGridEstimate )
      !deallocate( nGridEstimate )

    end subroutine prReset


    subroutine prAllocateArrays( nComputeBins,  &
                                inkernelSmoothing,&
                                inkernelSmoothingScale,&
                                inkernelSmoothingShape,&
                                inkernelSigmaSupport,&
                                inkernelSigmaSupportScale,&
                                incurvatureBandwidth,&
                                indensityEstimateArray ,&
                                innEstimateArray,&
                                inroughnessXXArray,&
                                inroughnessYYArray,&
                                inroughnessZZArray,&
                                innetRoughnessArray,&
                                activeGridCellsIn   )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        integer, intent(in) :: nComputeBins
        doubleprecision, dimension(:,:), allocatable, intent(out) :: inkernelSmoothing
        doubleprecision, dimension(:)  , allocatable, intent(out) :: inkernelSmoothingScale
        doubleprecision, dimension(:,:), allocatable, intent(out) :: inkernelSmoothingShape
        doubleprecision, dimension(:,:), allocatable, intent(out) :: inkernelSigmaSupport
        doubleprecision, dimension(:)  , allocatable, intent(out) :: inkernelSigmaSupportScale
        doubleprecision, dimension(:,:), allocatable, intent(out) :: incurvatureBandwidth
        doubleprecision, dimension(:)  , allocatable, intent(out) :: indensityEstimateArray 
        doubleprecision, dimension(:)  , allocatable, intent(out) :: innEstimateArray
        doubleprecision, dimension(:)  , allocatable, intent(out) :: inroughnessXXArray
        doubleprecision, dimension(:)  , allocatable, intent(out) :: inroughnessYYArray
        doubleprecision, dimension(:)  , allocatable, intent(out) :: inroughnessZZArray
        doubleprecision, dimension(:)  , allocatable, intent(out) :: innetRoughnessArray
        doubleprecision, dimension(:,:), allocatable :: lockernelSmoothing
        doubleprecision, dimension(:)  , allocatable :: lockernelSmoothingScale
        doubleprecision, dimension(:,:), allocatable :: lockernelSmoothingShape
        doubleprecision, dimension(:,:), allocatable :: lockernelSigmaSupport
        doubleprecision, dimension(:)  , allocatable :: lockernelSigmaSupportScale
        doubleprecision, dimension(:,:), allocatable :: loccurvatureBandwidth
        doubleprecision, dimension(:)  , allocatable :: locdensityEstimateArray 
        doubleprecision, dimension(:)  , allocatable :: locnEstimateArray
        doubleprecision, dimension(:)  , allocatable :: locroughnessXXArray
        doubleprecision, dimension(:)  , allocatable :: locroughnessYYArray
        doubleprecision, dimension(:)  , allocatable :: locroughnessZZArray
        doubleprecision, dimension(:)  , allocatable :: locnetRoughnessArray
        type( GridCellType ), dimension(:), allocatable, intent(out) :: activeGridCellsIn
        type( GridCellType ), dimension(:), allocatable :: activeGridCellsLocal
        !------------------------------------------------------------------------------


        !! Allocate arrays
        allocate(         lockernelSmoothing( 3, nComputeBins ) )
        allocate(      lockernelSigmaSupport( 3, nComputeBins ) )
        allocate(    lockernelSmoothingShape( 3, nComputeBins ) )
        allocate(      loccurvatureBandwidth( 3, nComputeBins ) )
        !allocate(         lockernelSmoothing( nDim, nComputeBins ) )
        !allocate(      lockernelSigmaSupport( nDim, nComputeBins ) )
        !allocate(    lockernelSmoothingShape( nDim, nComputeBins ) )
        !allocate(      loccurvatureBandwidth( nDim, nComputeBins ) )
        allocate(          lockernelSmoothingScale( nComputeBins ) )
        allocate(       lockernelSigmaSupportScale( nComputeBins ) )
        allocate(          locdensityEstimateArray( nComputeBins ) )
        allocate(                locnEstimateArray( nComputeBins ) )
        allocate(              locroughnessXXArray( nComputeBins ) )  
        allocate(              locroughnessYYArray( nComputeBins ) )
        allocate(              locroughnessZZArray( nComputeBins ) )
        allocate(             locnetRoughnessArray( nComputeBins ) )
        allocate(             activeGridCellsLocal( nComputeBins ) )

        call move_alloc(        activeGridCellsLocal,       activeGridCellsIn  )
        call move_alloc(          lockernelSmoothing,        inkernelSmoothing )
        call move_alloc(       lockernelSigmaSupport,     inkernelSigmaSupport )
        call move_alloc(     lockernelSmoothingShape,   inkernelSmoothingShape )
        call move_alloc(       loccurvatureBandwidth,     incurvatureBandwidth )
        call move_alloc(     lockernelSmoothingScale,   inkernelSmoothingScale )
        call move_alloc(  lockernelSigmaSupportScale,inkernelSigmaSupportScale )
        call move_alloc(     locdensityEstimateArray,   indensityEstimateArray )
        call move_alloc(           locnEstimateArray,         innEstimateArray )
        call move_alloc(         locroughnessXXArray,       inroughnessXXArray )  
        call move_alloc(         locroughnessYYArray,       inroughnessYYArray )
        call move_alloc(         locroughnessZZArray,       inroughnessZZArray )
        call move_alloc(        locnetRoughnessArray,      innetRoughnessArray )



    end subroutine prAllocateArrays


    subroutine prInitializeModuleDimensions( this, nDim, dimensionMask )
      !------------------------------------------------------------------------------
      ! 
      !
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      class( GridProjectedKDEType ), target :: this 
      integer, intent(inout)        :: nDim
      integer, dimension(3), intent(inout) :: dimensionMask
      integer :: n, nd, currentDim
      !------------------------------------------------------------------------------

      ! Determine dimensions based on number of bins
      do n = 1,3
        if (this%nBins(n) .eq. 1) dimensionMask(n) = 0 
      end do 
      nDim = sum(dimensionMask)
      this%dimensionMask = dimensionMask

      if ( nDim .le. 0 ) then 
        write(*,*) 'Error while initializing GPKDE dimensions. nDim .le. 0. Stop.'
        stop
      end if 

      ! Identify directions
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

      return

    end subroutine prInitializeModuleDimensions 


    subroutine prInitializeModuleConstants( this )
      !------------------------------------------------------------------------------
      ! 
      !
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      class( GridProjectedKDEType ) :: this 
      !------------------------------------------------------------------------------

      ! Compute constants
      this%supportDimensionConstant = ( ( nDim + 2d0 )*( 8d0*pi )**( 0.5*nDim ) )**( 0.25 )

      this%alphaDimensionConstant = ( ( 1d0 + 2d0**(0.5*nDim + 2d0) )/( 3d0*2d0**( 4d0/( nDim + 4d0 ) ) ) )**(&
          1d0/(nDim + 6d0) )*( nDim + 2d0 )**( 1d0/(nDim + 4d0) )/( ( nDim + 4d0 )**( 1d0/(nDim + 6d0) ) )

      this%betaDimensionConstant  = 2d0/( nDim + 4d0)/( nDim + 6d0 ) 

      return

    end subroutine prInitializeModuleConstants 



    ! NET ROUGHNESS
    ! net roughness
    ! initialize
    subroutine prInitializeNetRoughnessFunction( this, nDim )
      !------------------------------------------------------------------------------
      ! 
      !
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      class( GridProjectedKDEType ), target :: this 
      integer, intent(in)        :: nDim
      !------------------------------------------------------------------------------

      if ( nDim .eq. 1 ) then 
        ! Assign interfaces
        this%ComputeNetRoughnessEstimate => prComputeNetRoughness1D
        if ( this%databaseOptimization ) then 
          this%SetKernelSD => prSetKernelSD1DFromDatabase
        else
          this%SetKernelSD => prSetKernelSD1DBrute
        end if 
      end if

      if ( nDim .eq. 2 ) then
        ! Assign interface
        this%ComputeNetRoughnessEstimate => prComputeNetRoughness2D
        if ( this%databaseOptimization ) then 
          this%SetKernelSD2D => prSetKernelSD2DFromDatabase
        else
          this%SetKernelSD2D => prSetKernelSD2DBrute
        end if 
      end if

      if ( nDim .eq. 3 ) then 
        ! Assign interface
        this%ComputeNetRoughnessEstimate => prComputeNetRoughness3D
        if ( this%databaseOptimization ) then 
            this%SetKernelSD3D => prSetKernelSD3DFromDatabase
        else
            this%SetKernelSD3D => prSetKernelSD3DBrute
        end if 
      end if

      return

    end subroutine prInitializeNetRoughnessFunction 



    ! NET ROUGHNESS
    ! net roughness
    ! 1D
    subroutine prComputeNetRoughness1D( this, activeGridCells, curvatureBandwidth, &
                             roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                            netRoughnessArray, kernelSigmaSupport, &
                                                   kernelSDX, kernelSDY, kernelSDZ ) 
        !------------------------------------------------------------------------------
        ! Net roughness in 1D
        ! 
        !  - Eq. 13a in Sole-Mari et al. (2019)
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        !input
        class( GridProjectedKDEType ), target                   :: this
        type( GridCellType ), dimension(:), intent(in), target  :: activeGridCells
        doubleprecision, dimension(:,:), intent(in)             :: curvatureBandwidth
        doubleprecision, dimension(:,:), intent(in)             :: kernelSigmaSupport
        type( KernelSecondDerivativeXType ), intent(inout)      :: kernelSDX
        type( KernelSecondDerivativeYType ), intent(inout)      :: kernelSDY
        type( KernelSecondDerivativeZType ), intent(inout)      :: kernelSDZ
        ! out
        doubleprecision, dimension(:), intent(inout), target :: roughnessXXArray
        doubleprecision, dimension(:), intent(inout), target :: roughnessYYArray
        doubleprecision, dimension(:), intent(inout), target :: roughnessZZArray
        doubleprecision, dimension(:), intent(inout)         :: netRoughnessArray
        ! local 
        type( GridCellType ), pointer :: gc => null()
        doubleprecision, dimension(:), pointer :: roughness11Array
        doubleprecision, dimension(:,:,:), allocatable, target :: curvature1
        doubleprecision, dimension(:,:,:), allocatable, target :: curvature11
        doubleprecision, dimension(:,:,:), allocatable, target :: roughness11
        integer :: n
        integer :: iX, iY, iZ
        type( KernelMultiGaussianType ) :: kernelSigma
        !------------------------------------------------------------------------------
        allocate( curvature1(  this%nBins(1), this%nBins(2), this%nBins(3) )) 
        allocate( curvature11( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        allocate( roughness11( this%nBins(1), this%nBins(2), this%nBins(3) )) 
        !------------------------------------------------------------------------------
   

        ! Initialize
        curvature1  = 0d0
        roughness11 = 0d0

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
            !$omp private( gc )                       
            do n = 1, this%nComputeBins
    
                ! Assign gc pointer 
                gc => activeGridCells(n)
  
                if ( ( curvatureBandwidth(this%idDim1,n) .lt. 0d0 ) .or. & 
                     ( curvatureBandwidth(this%idDim1,n) /= curvatureBandwidth(this%idDim1,n) ) ) cycle

                ! Set kernel
                call this%SetKernelSD( gc, kernelSDX, curvatureBandwidth(:,n) )

                ! Compute curvature
                curvature1( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) = curvature1( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
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
            !$omp private( gc )                       
            do n = 1, this%nComputeBins
    
                ! Assign gc pointer 
                gc => activeGridCells(n)
  
                if ( ( curvatureBandwidth(this%idDim1,n) .lt. 0d0 ) .or. & 
                     ( curvatureBandwidth(this%idDim1,n) /= curvatureBandwidth(this%idDim1,n) ) ) cycle

                ! Set kernel
                call this%SetKernelSD( gc, kernelSDY, curvatureBandwidth(:,n) )

                ! Compute curvature
                curvature1( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) = curvature1( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
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
            !$omp private( gc )                       
            do n = 1, this%nComputeBins
    
                ! Assign gc pointer 
                gc => activeGridCells(n)
  
                if ( ( curvatureBandwidth(this%idDim1,n) .lt. 0d0 ) .or. & 
                     ( curvatureBandwidth(this%idDim1,n) /= curvatureBandwidth(this%idDim1,n) ) ) cycle

                ! Set kernel
                call this%SetKernelSD( gc, kernelSDZ, curvatureBandwidth(:,n) )

                ! Compute curvature
                curvature1( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) = curvature1( &
                        gc%kernelSDXGSpan(1):gc%kernelSDXGSpan(2), &
                        gc%kernelSDYGSpan(1):gc%kernelSDYGSpan(2), & 
                        gc%kernelSDZGSpan(1):gc%kernelSDZGSpan(2)  & 
                    ) + this%histogram%counts(                             &
                         gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(&
                                gc%kernelSDXMSpan(1):gc%kernelSDXMSpan(2), &
                                gc%kernelSDYMSpan(1):gc%kernelSDYMSpan(2), & 
                                gc%kernelSDZMSpan(1):gc%kernelSDZMSpan(2)  & 
                         )

            end do
            !$omp end parallel do

        end select
        curvature1 = curvature1/this%histogram%binVolume
        ! Matrix from curvature kernels is lambda**2*KernelVMatrix
        curvature1 = curvature1/( this%binSize(this%idDim1)**2 )


        ! Product curvatures, roughness
        curvature11 = curvature1*curvature1


        ! kernelSigma was already computed ? 
        call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange   )


        ! 11
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( activeGridCells )          &
        !$omp shared( curvature11 )              &
        !$omp shared( roughness11 )              & 
        !$omp shared( kernelSigmaSupport )       & 
        !$omp firstprivate( kernelSigma )        &
        !$omp private( gc )
        do n = 1, this%nComputeBins

            ! Assign pointer 
            gc => activeGridCells(n)

            call this%SetKernelSigma( gc, kernelSigma, kernelSigmaSupport( :, n ) )

            ! Compute roughness grid estimates
            roughness11( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
                curvature11(&
                    gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                    gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                    gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
                )*gc%kernelSigmaMatrix(&
                    gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                    gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                    gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

        end do
        !$omp end parallel do 


        ! Net roughness
        roughness11Array  = 0d0 
        netRoughnessArray = 0d0
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( activeGridCells )          &
        !$omp shared( roughness11 )              &
        !$omp shared( roughness11Array )         &
        !$omp shared( netRoughnessArray )        &
        !$omp private( gc )                      & 
        !$omp private( iX, iY, iZ )  
        do n = 1, this%nComputeBins

            ! Assign pointer 
            gc => activeGridCells(n)

            if ( gc%skipKernelSigma ) cycle

            iX = gc%id(1)
            iY = gc%id(2)
            iZ = gc%id(3)

            ! Assign info for needed arrays 
            roughness11Array( n ) = roughness11(iX,iY,iZ)

            ! Compute net roughness
            ! 1D
            netRoughnessArray( n ) = roughness11(iX,iY,iZ)

        end do
        !$omp end parallel do
        

        ! Deallocate
        deallocate( curvature1  ) 
        deallocate( curvature11 ) 
        deallocate( roughness11 ) 


        return


    end subroutine prComputeNetRoughness1D


    ! NET ROUGHNESS
    ! net roughness
    ! 2D
    subroutine prComputeNetRoughness2D( this, activeGridCells, curvatureBandwidth, &
                             roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                            netRoughnessArray, kernelSigmaSupport, &
                                                   kernelSDX, kernelSDY, kernelSDZ ) 
      !------------------------------------------------------------------------------
      ! Net roughness in 2D
      ! 
      !  - Eq. 13b in Sole-Mari et al. (2019)
      ! 
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      ! input
      class( GridProjectedKDEType ), target                  :: this
      type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
      doubleprecision, dimension(:,:), intent(in)            :: curvatureBandwidth
      doubleprecision, dimension(:,:), intent(in)            :: kernelSigmaSupport
      type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
      type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
      type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
      ! out
      doubleprecision, dimension(:), intent(inout), target :: roughnessXXArray
      doubleprecision, dimension(:), intent(inout), target :: roughnessYYArray
      doubleprecision, dimension(:), intent(inout), target :: roughnessZZArray
      doubleprecision, dimension(:), intent(inout)         :: netRoughnessArray
      ! local
      type( GridCellType ), pointer :: gc => null()
      doubleprecision, dimension(:), pointer         :: roughness11Array
      doubleprecision, dimension(:), pointer         :: roughness22Array
      doubleprecision, dimension(:,:,:), allocatable :: curvature1
      doubleprecision, dimension(:,:,:), allocatable :: curvature2
      doubleprecision, dimension(:,:,:), allocatable :: curvature11
      doubleprecision, dimension(:,:,:), allocatable :: curvature22
      doubleprecision, dimension(:,:,:), allocatable :: curvature12
      doubleprecision, dimension(:,:,:), allocatable :: roughness11
      doubleprecision, dimension(:,:,:), allocatable :: roughness22
      doubleprecision, dimension(:,:,:), allocatable :: roughness12
      integer :: n
      integer :: iX, iY, iZ
      type( KernelMultiGaussianType ) :: kernelSigma
      !------------------------------------------------------------------------------
      allocate( curvature1( this%nBins(1), this%nBins(2), this%nBins(3)  )) 
      allocate( curvature2( this%nBins(1), this%nBins(2), this%nBins(3)  )) 
      allocate( curvature11( this%nBins(1), this%nBins(2), this%nBins(3) )) 
      allocate( curvature22( this%nBins(1), this%nBins(2), this%nBins(3) )) 
      allocate( curvature12( this%nBins(1), this%nBins(2), this%nBins(3) )) 
      allocate( roughness11( this%nBins(1), this%nBins(2), this%nBins(3) )) 
      allocate( roughness22( this%nBins(1), this%nBins(2), this%nBins(3) )) 
      allocate( roughness12( this%nBins(1), this%nBins(2), this%nBins(3) )) 
      !------------------------------------------------------------------------------


      ! Initialize 
      curvature1  = 0d0
      curvature2  = 0d0
      roughness11 = 0d0
      roughness22 = 0d0
      roughness12 = 0d0

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
        !$omp private( gc )                       
        do n = 1, this%nComputeBins
    
            ! Assign gc pointer 
            gc => activeGridCells(n)
  
            if ( ( any( curvatureBandwidth(:,n) .lt. 0d0 ) ) .or. & 
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
                ) + this%histogram%counts(                               &
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
                ) + this%histogram%counts(                               &
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
        !$omp private( gc )                       
        do n = 1, this%nComputeBins
    
          ! Assign gc pointer 
          gc => activeGridCells(n)
  
          if ( ( any( curvatureBandwidth( :, n ) .lt. 0d0 ) ) .or. & 
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
              ) + this%histogram%counts(                               &
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
              ) + this%histogram%counts(                               &
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
        !$omp private( gc )                       
        do n = 1, this%nComputeBins
    
          ! Assign gc pointer 
          gc => activeGridCells(n)
  
          if ( ( any( curvatureBandwidth( :, n ) .lt. 0d0 ) ) .or. & 
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
              ) + this%histogram%counts(                               &
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
              ) + this%histogram%counts(                               &
                     gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD2Matrix(&
                          gc%kernelSD2XMSpan(1):gc%kernelSD2XMSpan(2), &
                          gc%kernelSD2YMSpan(1):gc%kernelSD2YMSpan(2), & 
                          gc%kernelSD2ZMSpan(1):gc%kernelSD2ZMSpan(2)  & 
                     )

        end do
        !$omp end parallel do
      end if 


      curvature1 = curvature1/this%histogram%binVolume
      curvature2 = curvature2/this%histogram%binVolume
      ! Matrix from curvature kernels is lambda**2*KernelVMatrix
      curvature1 = curvature1/( this%binSize(this%idDim1)**2 )
      curvature2 = curvature2/( this%binSize(this%idDim2)**2 )

      ! Product curvatures, roughness
      curvature11 = curvature1*curvature1
      curvature22 = curvature2*curvature2
      curvature12 = curvature1*curvature2

      ! Initialize kernelSigma
      call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )

      ! 11
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvature11 )              &
      !$omp shared( roughness11 )              & 
      !$omp firstprivate( kernelSigma )        &
      !$omp shared( kernelSigmaSupport )       &
      !$omp private( gc )
      do n = 1, this%nComputeBins

          ! Assign pointer 
          gc => activeGridCells(n)

          ! Notice this stage: kernelSigma is used in following loops
          call this%SetKernelSigma( gc, kernelSigma, kernelSigmaSupport( :, n ) )

          ! Compute roughness grid estimates
          roughness11( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
              curvature11(&
                  gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                  gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                  gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
              )*gc%kernelSigmaMatrix(&
                  gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                  gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                  gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

      end do
      !$omp end parallel do 


      ! 22
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvature22 )              &
      !$omp shared( roughness22 )              & 
      !$omp private( gc )
      do n = 1, this%nComputeBins

          ! Assign pointer 
          gc => activeGridCells(n)

          ! Compute roughness grid estimates
          roughness22( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
              curvature22(&
                  gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                  gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                  gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
              )*gc%kernelSigmaMatrix(&
                  gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                  gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                  gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

      end do
      !$omp end parallel do 


      ! 12
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( curvature12 )              &
      !$omp shared( roughness12 )              & 
      !$omp private( gc )
      do n = 1, this%nComputeBins

          ! Assign pointer 
          gc => activeGridCells(n)

          ! Compute roughness grid estimates
          roughness12( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
              curvature12(&
                  gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                  gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                  gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
              )*gc%kernelSigmaMatrix(&
                  gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                  gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                  gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

      end do
      !$omp end parallel do 


      ! Net roughness
      roughness11Array  = 0d0 
      roughness22Array  = 0d0 
      netRoughnessArray = 0d0
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( roughness11, roughness22 ) &
      !$omp shared( roughness12 )              & 
      !$omp shared( roughness11Array )         &
      !$omp shared( roughness22Array )         &
      !$omp shared( netRoughnessArray )        &
      !$omp private( gc )                      & 
      !$omp private( iX, iY, iZ )  
      do n = 1, this%nComputeBins

          ! Assign pointer 
          gc => activeGridCells(n)

          iX = gc%id(1)
          iY = gc%id(2)
          iZ = gc%id(3)

          ! Assign info for needed arrays 
          roughness11Array( n ) = roughness11(iX,iY,iZ)
          roughness22Array( n ) = roughness22(iX,iY,iZ)

          ! Compute net roughness
          ! 2D
          ! ISOTROPIC
          !netRoughnessArray( n ) = roughness11(iX,iY,iZ) + 2*roughness12(iX,iY,iZ) + roughness22(iX,iY,iZ) 
          ! ANISOTROPIC
          !netRoughnessArray( n ) = 2*sqrt( roughness11(iX,iY,iZ)*roughness22(iX,iY,iZ) ) + 2*roughness12(iX,iY,iZ)
          ! Combined ( avoid smoothing blow ups )
          netRoughnessArray( n ) = 2*roughness12(iX,iY,iZ) + maxval( (/     &
                  roughness11(iX,iY,iZ) + roughness22(iX,iY,iZ),           &
                  2*sqrt( roughness11(iX,iY,iZ)*roughness22(iX,iY,iZ) ) /) )

      end do
      !$omp end parallel do
      

      ! Deallocate
      deallocate( curvature1  ) 
      deallocate( curvature2  ) 
      deallocate( curvature11 ) 
      deallocate( curvature22 ) 
      deallocate( curvature12 ) 
      deallocate( roughness11 ) 
      deallocate( roughness22 ) 
      deallocate( roughness12 ) 


      return


    end subroutine prComputeNetRoughness2D


    ! NET ROUGHNESS
    ! net roughness
    ! 3D
    subroutine prComputeNetRoughness3D( this, activeGridCells, curvatureBandwidth, &
                             roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                            netRoughnessArray, kernelSigmaSupport, &
                                                   kernelSDX, kernelSDY, kernelSDZ ) 
      !------------------------------------------------------------------------------
      ! Net roughness in 3D
      !
      !  - Eq. 13c in Sole-Mari et al. (2019)
      !
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      !input
      class( GridProjectedKDEType ), target :: this
      type( GridCellType ), dimension(:), intent(in), target :: activeGridCells
      doubleprecision, dimension(:,:), intent(in)            :: curvatureBandwidth
      doubleprecision, dimension(:,:), intent(in)            :: kernelSigmaSupport
      type( KernelSecondDerivativeXType ), intent(inout)     :: kernelSDX
      type( KernelSecondDerivativeYType ), intent(inout)     :: kernelSDY
      type( KernelSecondDerivativeZType ), intent(inout)     :: kernelSDZ
      ! out
      doubleprecision, dimension(:), intent(inout), target   :: roughnessXXArray
      doubleprecision, dimension(:), intent(inout), target   :: roughnessYYArray
      doubleprecision, dimension(:), intent(inout), target   :: roughnessZZArray
      doubleprecision, dimension(:), intent(inout)   :: netRoughnessArray
      ! local 
      type( GridCellType ), pointer :: gc => null()
      doubleprecision, dimension(:,:,:), pointer :: curvature
      doubleprecision, dimension(:,:,:), pointer :: roughness
      doubleprecision, dimension(:,:,:), allocatable, target ::   curvatureX
      doubleprecision, dimension(:,:,:), allocatable, target ::   curvatureY
      doubleprecision, dimension(:,:,:), allocatable, target ::   curvatureZ
      doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureXX
      doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureXY
      doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureXZ
      doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureYY
      doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureYZ
      doubleprecision, dimension(:,:,:), allocatable, target ::  curvatureZZ
      doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessXX
      doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessXY
      doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessXZ
      doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessYY
      doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessYZ
      doubleprecision, dimension(:,:,:), allocatable, target ::  roughnessZZ
      integer :: n, nr  
      integer :: iX, iY, iZ 
      !------------------------------------------------------------------------------
      allocate(          curvatureX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(          curvatureY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(          curvatureZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(         curvatureXX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(         curvatureXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(         curvatureXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(         curvatureYY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(         curvatureYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(         curvatureZZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(         roughnessXX( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(         roughnessXY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(         roughnessXZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(         roughnessYY( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(         roughnessYZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      allocate(         roughnessZZ( this%nBins(1), this%nBins(2), this%nBins(3) ) )
      !------------------------------------------------------------------------------


      ! Initialize 
      curvatureX  = 0d0
      curvatureY  = 0d0
      curvatureZ  = 0d0
      roughnessXX = 0d0 
      roughnessYY = 0d0
      roughnessZZ = 0d0
      roughnessXY = 0d0
      roughnessXZ = 0d0
      roughnessYZ = 0d0


      ! Curvatures, kappa
      !$omp parallel do schedule( dynamic, 1 )           & 
      !$omp default( none )                              &
      !$omp shared( this )                               &
      !$omp shared( activeGridCells )                    &
      !$omp reduction( +:curvatureX )                    &
      !$omp reduction( +:curvatureY )                    &
      !$omp reduction( +:curvatureZ )                    &
      !$omp shared( curvatureBandwidth )                 &
      !$omp firstprivate( kernelSDX )                    &
      !$omp firstprivate( kernelSDY )                    &
      !$omp firstprivate( kernelSDZ )                    &
      !$omp private( gc )                       
      do n = 1, this%nComputeBins
    
        ! Assign gc pointer 
        gc => activeGridCells(n)
  
        if ( ( any( curvatureBandwidth( :, n ) .lt. 0d0 ) ) .or. & 
             ( any( curvatureBandwidth( :, n ) /= curvatureBandwidth( :, n ) ) ) ) cycle

        ! Set kernels        
        call this%SetKernelSD3D( gc, kernelSDX, kernelSDY, kernelSDZ, curvatureBandwidth( :, n ) )

        ! Compute curvature
        curvatureX( &
                gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
            ) = curvatureX( &
                gc%kernelSD1XGSpan(1):gc%kernelSD1XGSpan(2), &
                gc%kernelSD1YGSpan(1):gc%kernelSD1YGSpan(2), & 
                gc%kernelSD1ZGSpan(1):gc%kernelSD1ZGSpan(2)  & 
            ) + this%histogram%counts(                               &
                gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD1Matrix(   &
                        gc%kernelSD1XMSpan(1):gc%kernelSD1XMSpan(2), &
                        gc%kernelSD1YMSpan(1):gc%kernelSD1YMSpan(2), & 
                        gc%kernelSD1ZMSpan(1):gc%kernelSD1ZMSpan(2)  & 
                )

        ! Compute curvature
        curvatureY( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) = curvatureY( &
                gc%kernelSD2XGSpan(1):gc%kernelSD2XGSpan(2), &
                gc%kernelSD2YGSpan(1):gc%kernelSD2YGSpan(2), & 
                gc%kernelSD2ZGSpan(1):gc%kernelSD2ZGSpan(2)  & 
            ) + this%histogram%counts(                               &
                   gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD2Matrix(&
                        gc%kernelSD2XMSpan(1):gc%kernelSD2XMSpan(2), &
                        gc%kernelSD2YMSpan(1):gc%kernelSD2YMSpan(2), & 
                        gc%kernelSD2ZMSpan(1):gc%kernelSD2ZMSpan(2)  & 
                   )
        
        ! Compute curvature 
        curvatureZ( &
                gc%kernelSD3XGSpan(1):gc%kernelSD3XGSpan(2), &
                gc%kernelSD3YGSpan(1):gc%kernelSD3YGSpan(2), & 
                gc%kernelSD3ZGSpan(1):gc%kernelSD3ZGSpan(2)  & 
            ) = curvatureZ( &
                gc%kernelSD3XGSpan(1):gc%kernelSD3XGSpan(2), &
                gc%kernelSD3YGSpan(1):gc%kernelSD3YGSpan(2), & 
                gc%kernelSD3ZGSpan(1):gc%kernelSD3ZGSpan(2)  & 
            ) + this%histogram%counts(                             &
                gc%id(1), gc%id(2), gc%id(3) )*gc%kernelSD3Matrix(&
                        gc%kernelSD3XMSpan(1):gc%kernelSD3XMSpan(2), &
                        gc%kernelSD3YMSpan(1):gc%kernelSD3YMSpan(2), & 
                        gc%kernelSD3ZMSpan(1):gc%kernelSD3ZMSpan(2)  & 
                )
    
      end do
      !$omp end parallel do
      curvatureX = curvatureX/this%histogram%binVolume
      curvatureY = curvatureY/this%histogram%binVolume
      curvatureZ = curvatureZ/this%histogram%binVolume
      ! Matrix from curvature kernels is lambda**2*KernelVMatrix
      curvatureX = curvatureX/( this%binSize(1)**2 )
      curvatureY = curvatureY/( this%binSize(2)**2 )
      curvatureZ = curvatureZ/( this%binSize(3)**2 )

      ! Compute curvatures product
      curvatureXX = curvatureX*curvatureX
      curvatureYY = curvatureY*curvatureY
      curvatureZZ = curvatureZ*curvatureZ
      curvatureXY = curvatureX*curvatureY
      curvatureXZ = curvatureX*curvatureZ
      curvatureYZ = curvatureY*curvatureZ

      ! Compute roughnesses
      do nr = 1, 6
        select case(nr) 
          case (1)
            roughness => roughnessXX
            curvature => curvatureXX
          case (2)
            roughness => roughnessYY
            curvature => curvatureYY
          case (3)
            roughness => roughnessZZ
            curvature => curvatureZZ
          case (4)
            roughness => roughnessXY
            curvature => curvatureXY
          case (5)
            roughness => roughnessXZ
            curvature => curvatureXZ
          case (6)
            roughness => roughnessYZ
            curvature => curvatureYZ
        end select

        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( activeGridCells )          &
        !$omp shared( curvature )                &
        !$omp shared( roughness )                & 
        !$omp private( gc )
        do n = 1, this%nComputeBins

          ! Assign pointer 
          gc => activeGridCells(n)

          ! Compute roughness grid estimates
          roughness( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
              curvature(&
                  gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                  gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                  gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
              )*gc%kernelSigmaMatrix(&
                  gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                  gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                  gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2) )) 

        end do
        !$omp end parallel do 
      end do

      ! Net roughness
      roughnessXXArray  = 0d0 
      roughnessYYArray  = 0d0 
      roughnessZZArray  = 0d0 
      netRoughnessArray = 0d0
      !$omp parallel do schedule( dynamic, 1 ) &
      !$omp default( none )                    &
      !$omp shared( this )                     &
      !$omp shared( activeGridCells )          &
      !$omp shared( roughnessXX, roughnessYY ) &
      !$omp shared( roughnessZZ, roughnessXY ) &
      !$omp shared( roughnessXZ, roughnessYZ ) &
      !$omp shared( roughnessXXArray )         &
      !$omp shared( roughnessYYArray )         &
      !$omp shared( roughnessZZArray )         &
      !$omp shared( netRoughnessArray )        &
      !$omp private( gc )                      & 
      !$omp private( iX, iY, iZ )  
      do n = 1, this%nComputeBins

        ! Assign pointer 
        gc => activeGridCells(n)

        iX = gc%id(1)
        iY = gc%id(2)
        iZ = gc%id(3)

        ! Assign info for needed arrays 
        roughnessXXArray( n ) = roughnessXX(iX,iY,iZ)
        roughnessYYArray( n ) = roughnessYY(iX,iY,iZ)
        roughnessZZArray( n ) = roughnessZZ(iX,iY,iZ)

        ! Compute net roughness
        ! 3D
        ! ISOTROPIC
        netRoughnessArray( n ) = roughnessXX(iX,iY,iZ) + 2*roughnessXY(iX,iY,iZ) + 2*roughnessXZ(iX,iY,iZ) + &
                                       roughnessYY(iX,iY,iZ) + 2*roughnessYZ(iX,iY,iZ) + roughnessZZ(iX,iY,iZ)

        ! ANISOTROPIC
        !netRoughnessArray( n ) = 3*( roughnessXX(iX,iY,iZ)*roughnessYY(iX,iY,iZ)*roughnessZZ(iX,iY,iZ) )**(1d0/3d0) + &
        !    2*roughnessYZ(iX,iY,iZ)*( roughnessXX(iX,iY,iZ)**2/roughnessYY(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6d0) + &
        !    2*roughnessXZ(iX,iY,iZ)*( roughnessYY(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessZZ(iX,iY,iZ) )**(1d0/6d0) + &
        !    2*roughnessXY(iX,iY,iZ)*( roughnessZZ(iX,iY,iZ)**2/roughnessXX(iX,iY,iZ)/roughnessYY(iX,iY,iZ) )**(1d0/6d0)

        ! NEED COMBINED !
      end do
      !$omp end parallel do
      
      ! Deallocate
      deallocate(          curvatureX )
      deallocate(          curvatureY )
      deallocate(          curvatureZ )
      deallocate(         curvatureXX )
      deallocate(         curvatureXY )
      deallocate(         curvatureXZ )
      deallocate(         curvatureYY )
      deallocate(         curvatureYZ )
      deallocate(         curvatureZZ )
      deallocate(         roughnessXX )
      deallocate(         roughnessXY )
      deallocate(         roughnessXZ )
      deallocate(         roughnessYY )
      deallocate(         roughnessYZ )
      deallocate(         roughnessZZ )


      ! Done
      return


    end subroutine prComputeNetRoughness3D



    subroutine prInitializeKernelDatabaseFlat( this,  &
                      minHOverLambda, maxHOverLambda, &
                 deltaHOverLambda, logKernelDatabase, &
                           kernelRange, kernelSDRange )
      !------------------------------------------------------------------------------
      ! 
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      implicit none
      class( GridProjectedKDEType ), target :: this
      ! input
      doubleprecision,   intent(in) :: minHOverLambda
      doubleprecision,   intent(in) :: maxHOverLambda
      doubleprecision,   intent(in) :: deltaHOverLambda
      logical, intent(in), optional :: logKernelDatabase
      integer, intent(in), optional :: kernelRange
      integer, intent(in), optional :: kernelSDRange
      ! local
      doubleprecision, dimension(3) :: inputSmoothing
      doubleprecision, dimension(:), allocatable :: hOverLambda
      integer :: nDelta
      integer :: i, n, m, o, dbi
      logical :: localLogDatabase
      integer :: localKernelRange
      integer :: localKernelSDRange

      ! Mem debug
      doubleprecision :: kernelMatrixMemory
      doubleprecision :: kernelDBMemory    
      doubleprecision :: kernelSDDBMemory  

      ! Time monitoring
      integer         :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
      doubleprecision :: elapsedTime
      !------------------------------------------------------------------------------

      ! Sanity check for input parameters

      ! Default parameters
      ! logDatabase as false
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
        ! Verify this 
        nDelta      = ceiling( log10( maxHOverLambda/minHOverLambda )/log10( 1 + deltaHOverLambda ) ) + 1
        allocate( hOverLambda( nDelta ) )
        hOverLambda = this%GenerateLogSpaceData( minHOverLambda, maxHOverLambda, nDelta )

        ! Assign indexes interfaces
        this%ComputeKernelDatabaseFlatIndexes => prComputeKernelDatabaseFlatIndexesLog
        this%ComputeKernelDatabaseIndexes     => prComputeKernelDatabaseIndexesLog ! meanwhile for SD's

        this%deltaHOverLambda = log( hOverLambda(2)/hOverLambda(1) ) ! Fix this inconsistency, is overwritten
      else 
        ! LINEAR FORM
        nDelta      = floor( ( maxHOverLambda - minHOverLambda )/deltaHOverLambda )
        allocate( hOverLambda( nDelta ) )
        hOverLambda = [ (minHOverLambda + i*deltaHOverLambda, i=0, nDelta ) ]

        ! Assign indexes interface
        this%ComputeKernelDatabaseFlatIndexes => prComputeKernelDatabaseFlatIndexesLinear
        this%ComputeKernelDatabaseIndexes     => prComputeKernelDatabaseIndexesLinear ! meanwhile for SD's
      end if 

      ! Assign to the object
      ! Temporarilly the same value for each axis
      this%nDeltaHOverLambda   = nDelta

      ! Depending on the number of dimensions
      ! is the required kernel database.
      select case(nDim)
      !1D
      case(1)
        ! Allocate kernel databases
        allocate( this%kernelDatabaseFlat( nDelta, 1 ) )

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
        kernelMatrixMemory = 0d0
        kernelDBMemory = 0d0
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( hOverLambda )              &
        !$omp shared( nDelta )                   &
        !$omp shared( localKernelRange )         &
        !$omp reduction( +:kernelDBMemory )      &
        !$omp private( kernelMatrixMemory )      &
        !$omp private( inputSmoothing )
        do n = 1, nDelta
          inputSmoothing(:) = 0
          inputSmoothing( this%idDim1 ) = hOverLambda(n)
          call this%kernelDatabaseFlat( n, 1 )%Initialize( &
            this%binSize, matrixRange=localKernelRange, dimensionMask=this%dimensionMask )
          call this%kernelDatabaseFlat( n, 1 )%SetupMatrix( inputSmoothing*this%binSize )
          kernelMatrixMemory = sizeof( this%kernelDatabaseFlat( n, 1 )%matrix )/1d6
          kernelDBMemory     = kernelDBMemory + kernelMatrixMemory
        end do
        !$omp end parallel do
        ! TOC
        call system_clock(clockCountStop, clockCountRate, clockCountMax)
        elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)
        
        ! TIC
        call system_clock(clockCountStart, clockCountRate, clockCountMax)
        ! Second derivatives
        kernelMatrixMemory = 0d0
        kernelSDDBMemory = 0d0
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( hOverLambda )              &
        !$omp shared( nDelta )                   &
        !$omp shared( localKernelSDRange )       &
        !$omp reduction( +:kernelSDDBMemory )    &
        !$omp private( kernelMatrixMemory )      &
        !$omp private( inputSmoothing )
        do n = 1, nDelta
          inputSmoothing(:) = 0
          inputSmoothing( this%idDim1 ) = hOverLambda(n)
          ! 1 
          call this%kernelSDDatabase1( n )%Initialize(& 
              this%binSize, matrixRange=localKernelSDRange, dimensionMask=this%dimensionMask )
          call this%kernelSDDatabase1( n )%SetupMatrix( inputSmoothing*this%binSize )
          kernelMatrixMemory = sizeof( this%kernelSDDatabase1( n )%matrix )/1d6
          kernelSDDBMemory   = kernelSDDBMemory + kernelMatrixMemory
        end do
        !$omp end parallel do
        ! TOC
        call system_clock(clockCountStop, clockCountRate, clockCountMax)
        elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

      ! 2D
      case(2)
        allocate( this%kernelDatabaseFlat( nDelta*( nDelta + 1 )/2, 1 ) )

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
        kernelMatrixMemory = 0d0
        kernelDBMemory = 0d0
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( hOverLambda )              &
        !$omp shared( nDelta )                   &
        !$omp shared( localKernelRange )         &
        !$omp private( m, dbi )                  &
        !$omp reduction( +:kernelDBMemory )      &
        !$omp private( kernelMatrixMemory )      &
        !$omp private( inputSmoothing )
        do n = 1, nDelta
          do m = 1, min( n, nDelta )
            dbi = n*( n - 1 )/2 + m
            inputSmoothing(:) = 0d0
            inputSmoothing( this%idDim1 ) =  hOverLambda(n)
            inputSmoothing( this%idDim2 ) =  hOverLambda(m)
            call this%kernelDatabaseFlat( dbi, 1 )%Initialize( & 
                this%binSize, matrixRange=localKernelRange )
            call this%kernelDatabaseFlat( dbi, 1 )%SetupMatrix( inputSmoothing*this%binSize )
            kernelMatrixMemory = sizeof( this%kernelDatabaseFlat( dbi, 1 )%matrix )/1d6
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
        kernelMatrixMemory = 0d0
        kernelSDDBMemory = 0d0
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( hOverLambda )              &
        !$omp shared( nDelta )                   &
        !$omp shared( localKernelSDRange )       &
        !$omp reduction( +:kernelSDDBMemory )    &
        !$omp private( kernelMatrixMemory )      &
        !$omp private( inputSmoothing )
        do n = 1, nDelta
          inputSmoothing(:) = 0
          inputSmoothing( this%idDim1 ) = hOverLambda(n)
          inputSmoothing( this%idDim2 ) = hOverLambda(n)
          ! 1 
          call this%kernelSDDatabase1( n )%Initialize(& 
              this%binSize, matrixRange=localKernelSDRange, dimensionMask=this%dimensionMask )
          call this%kernelSDDatabase1( n )%SetupMatrix( inputSmoothing*this%binSize )
          kernelMatrixMemory = sizeof( this%kernelSDDatabase1( n )%matrix )/1d6
          kernelSDDBMemory   = kernelSDDBMemory + kernelMatrixMemory
          ! 2 
          call this%kernelSDDatabase2( n )%Initialize(& 
              this%binSize, matrixRange=localKernelSDRange, dimensionMask=this%dimensionMask )
          call this%kernelSDDatabase2( n )%SetupMatrix( inputSmoothing*this%binSize )
          kernelMatrixMemory = sizeof( this%kernelSDDatabase1( n )%matrix )/1d6
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

        ! TIC
        call system_clock(clockCountStart, clockCountRate, clockCountMax)
        ! Kernel database
        kernelMatrixMemory = 0d0
        kernelDBMemory = 0d0
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( hOverLambda )              &
        !$omp shared( nDelta )                   &
        !$omp shared( localKernelRange )         &
        !$omp private( n, m, dbi )               &
        !$omp reduction( +:kernelDBMemory )      &
        !$omp private( kernelMatrixMemory )      &
        !$omp private( inputSmoothing )
        do o = 1, nDelta
          do n = 1, nDelta
            do m = 1, min( n, nDelta )
              dbi = n*( n - 1 )/2 + m
              inputSmoothing = (/ hOverLambda(n), hOverLambda(m), hOverLambda(o) /) 
              call this%kernelDatabaseFlat( dbi, o )%Initialize( & 
                      this%binSize, matrixRange=localKernelRange )
              call this%kernelDatabaseFlat( dbi, o )%SetupMatrix( inputSmoothing*this%binSize )
              kernelMatrixMemory = sizeof( this%kernelDatabaseFlat( dbi, o )%matrix )/1d6
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
        kernelMatrixMemory = 0d0
        kernelSDDBMemory = 0d0
        !$omp parallel do schedule( dynamic, 1 ) &
        !$omp default( none )                    &
        !$omp shared( this )                     &
        !$omp shared( hOverLambda )              &
        !$omp shared( nDelta )                   &
        !$omp shared( localKernelSDRange )       &
        !$omp reduction( +:kernelSDDBMemory )    &
        !$omp private( kernelMatrixMemory )      &
        !$omp private( inputSmoothing )
        do n = 1, nDelta
          inputSmoothing = (/ hOverLambda(n), hOverLambda(n), hOverLambda(n) /)
          ! X 
          call this%kernelSDXDatabase( n )%Initialize(& 
              this%binSize, matrixRange=localKernelSDRange )
          call this%kernelSDXDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )
          kernelMatrixMemory = sizeof( this%kernelSDXDatabase( n )%matrix )/1d6
          kernelSDDBMemory = kernelSDDBMemory + kernelMatrixMemory
          ! Y
          call this%kernelSDYDatabase( n )%Initialize(& 
              this%binSize, matrixRange=localKernelSDRange )
          call this%kernelSDYDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )
          kernelMatrixMemory = sizeof( this%kernelSDYDatabase( n )%matrix )/1d6
          kernelSDDBMemory = kernelSDDBMemory + kernelMatrixMemory
          ! Z
          call this%kernelSDZDatabase( n )%Initialize(& 
              this%binSize, matrixRange=localKernelSDRange )
          call this%kernelSDZDatabase( n )%SetupMatrix( inputSmoothing*this%binSize )
          kernelMatrixMemory = sizeof( this%kernelSDZDatabase( n )%matrix )/1d6
          kernelSDDBMemory = kernelSDDBMemory + kernelMatrixMemory
        end do
        !$omp end parallel do
        ! TOC
        call system_clock(clockCountStop, clockCountRate, clockCountMax)
        elapsedTime = dble(clockCountStop - clockCountStart) / dble(clockCountRate)

      end select

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

        ! Dropping database does not mean 
        ! that parameters are resetted
        !this%deltaHOverLambda  = 0d0 
        !this%nDeltaHOverLambda = 0d0
        !this%minHOverLambda    = 0d0
        !this%maxHOverLambda    = 0d0

        return


    end subroutine prDropKernelDatabase


    ! Density computation manager 
    subroutine prComputeDensity( this, dataPoints, nOptimizationLoops, &
        outputFileName, outputFileUnit, outputDataId, particleGroupId, &
                persistentKernelDatabase, exportOptimizationVariables, & 
                                     skipErrorConvergence, unitVolume, &
                                scalingFactor, histogramScalingFactor, &
                                                    computeRawDensity, &
                weightedHistogram, weights, onlyHistogram, exactPoint  ) 
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        ! input
        class( GridProjectedKDEType ), target               :: this
        doubleprecision, dimension(:,:), intent(in)         :: dataPoints
        integer, intent(in), optional                       :: nOptimizationLoops
        character(len=*), intent(in), optional              :: outputFileName
        integer, intent(in), optional                       :: outputFileUnit
        integer, intent(in), optional                       :: outputDataId
        integer, intent(in), optional                       :: particleGroupId
        logical, intent(in), optional                       :: persistentKernelDatabase
        logical, intent(in), optional                       :: exportOptimizationVariables
        logical, intent(in), optional                       :: skipErrorConvergence
        logical, intent(in), optional                       :: unitVolume
        doubleprecision, intent(in), optional               :: scalingFactor
        doubleprecision, intent(in), optional               :: histogramScalingFactor
        logical, intent(in), optional                       :: computeRawDensity
        logical, intent(in), optional                       :: weightedHistogram
        logical, intent(in), optional                       :: onlyHistogram
        logical, intent(in), optional                       :: exactPoint
        doubleprecision, dimension(:), intent(in), optional :: weights
        ! local 
        logical               :: persistKDB
        logical               :: locExportOptimizationVariables
        logical               :: locSkipErrorConvergence
        logical               :: locUnitVolume
        doubleprecision       :: locScalingFactor
        logical               :: locScaleHistogram
        logical               :: locComputeRawDensity
        doubleprecision       :: locHistogramScalingFactor
        logical               :: locWeightedHistogram
        logical               :: locOnlyHistogram 
        logical               :: locExactPoint 
        integer               :: localNOptimizationLoops
        integer, dimension(2) :: dataPointsShape
        character(len=16)     :: timeChar
        character(len=16)     :: spcChar

        ! DEV (DEPRECATE)
        logical :: useBoundingBox = .false.
        type( KernelMultiGaussianType ) :: filterKernel
        integer, dimension(2) :: xGridSpan, yGridSpan, zGridSpan
        ! Kernel span are not used but required by filterKernel%ComputeGridSpans function
        integer, dimension(2) :: xKernelSpan, yKernelSpan, zKernelSpan
        integer :: n
        integer :: bcount = 1
        logical, dimension(:), allocatable :: computeThisBin
        !! Time monitoring
        !integer         :: clockCountStart, clockCountStop, clockCountRate, clockCountMax
        !doubleprecision :: elapsedTime
        !------------------------------------------------------------------------------

        ! Initialize optional arguments
        persistKDB = .true.
        locExportOptimizationVariables =.false.
        locSkipErrorConvergence =.false.
        locUnitVolume =.false.
        locScalingFactor = 1d0
        locScaleHistogram = .false.
        locComputeRawDensity = .false.
        locHistogramScalingFactor = 1d0
        locWeightedHistogram = .false.
        locOnlyHistogram = .false.
        locExactPoint = .false.
        localNOptimizationLoops = this%nOptimizationLoops

        ! Process them
        if ( present( nOptimizationLoops ) ) then 
          localNOptimizationLoops = nOptimizationLoops
        end if 
        if ( present( outputFileName ) ) then 
          this%outputFileName = outputFileName
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
        if ( (locWeightedHistogram).and.(.not.present(weights)) ) then 
          write(*,*) 'ERROR: weightedHistogram requires weights and were not given. Stop.'
          stop
        end if 

        dataPointsShape = shape(dataPoints)
        if ( (locWeightedHistogram).and.(size(weights).ne.dataPointsShape(1)) ) then 
          write(*,*) 'ERROR: given weights are not the same length than datapoints. Stop.'
          stop
        end if

        if ( locWeightedHistogram ) then 
          ! Cummulative histogram-like quantities
          call this%histogram%ComputeCountsWeighted( dataPoints, weights, locExactPoint )
        else
          ! Histogram quantities
          call this%histogram%ComputeCounts( dataPoints, locExactPoint )
        end if

        ! If only histogram, leave
        if ( locOnlyHistogram ) then 
          if ( locComputeRawDensity ) then 
            this%histogram%counts = this%histogram%counts/this%histogram%binVolume
          end if 
          return
        end if 


        ! SOON TO BE DEPRECATED !
        ! Bounding box or active bins
        ! Not in use !
        if ( useBoundingBox ) then 

            ! Compute bounding box
            call this%histogram%ComputeBoundingBox()
        
            ! Initialize filterKernel
            call filterKernel%Initialize(  this%binSize, matrixRange=defaultKernelRange )

            ! This could be a factor times the initial smoothing, needs elegance
            call filterKernel%SetupMatrix( 0.5*this%initialSmoothing ) 
            
            ! Allocate the identifier 
            allocate( computeThisBin( this%histogram%nBBoxBins ) )
            computeThisBin = .false.

            ! Now loop over the cells within the bounding box,
            ! and count how many bins will be computed.
            !$omp parallel do                                &
            !$omp firstprivate( filterKernel )               &
            !$omp private( xGridSpan, yGridSpan, zGridSpan ) &  
            !$omp private( xKernelSpan, yKernelSpan, zKernelSpan ) 
            do n = 1, this%histogram%nBBoxBins

                ! Determine spans
                call filterKernel%ComputeGridSpans(&
                    this%histogram%boundingBoxBinIds( :, n ), this%nBins, &
                                         xGridSpan, yGridSpan, zGridSpan, & 
                                   xKernelSpan, yKernelSpan, zKernelSpan, &
                                   this%dimensionMask  ) 

                if ( any( this%histogram%counts(    &
                        xGridSpan(1):xGridSpan(2),  &
                        yGridSpan(1):yGridSpan(2),  & 
                        zGridSpan(1):zGridSpan(2) ) .gt. 0 ) ) then ! Cell active

                   computeThisBin( n ) = .true.

                end if

            end do
            !$omp end parallel do 

            ! Count how many and allocate
            this%nComputeBins = count( computeThisBin )
            allocate( this%computeBinIds( nDim, this%nComputeBins ) )

            ! Fill computeBinIds
            do n = 1, this%histogram%nBBoxBins
                if ( computeThisBin( n ) ) then 
                    this%computeBinIds( :, bcount ) = this%histogram%boundingBoxBinIds( :, n )
                    bcount = bcount + 1
                end if 
            end do 

            deallocate( computeThisBin )

        else
          ! Active bins: Only cells with particles
          call this%histogram%ComputeActiveBinIds()

          this%computeBinIds => this%histogram%activeBinIds
          this%nComputeBins  = this%histogram%nActiveBins
              
          if ( this%nComputeBins .eq. 0 ) then 
           ! No bins to compute 
           if ( this%reportToOutUnit ) then
            write(this%outFileUnit, *  )
            write(this%outFileUnit, '(A)' ) 'WARNING: GPKDE module  '
            write(this%outFileUnit, '(A)' ) 'NO bins to compute. Check origin coordinates or particles. Leaving ComputeDensity.'
            write(this%outFileUnit, *  )
           end if
           ! Leaving  
           return
          end if 

        end if 


        if ( this%reportToOutUnit ) then 
          if ( present( outputDataId ) .and. present( particleGroupId ) ) then 
          timeChar=''
          spcChar = ''
          write(timeChar,*)outputDataId
          write(spcChar,*)particleGroupId
          write( this%outFileUnit, '(A,A,A,A)' )' GPKDE Optimization -- Time: ', trim(adjustl(timeChar)), &
                  ' -- Specie: ', trim(adjustl(spcChar))
          else
            write( this%outFileUnit, '(A)' )' GPKDE Optimization ---------------------------'
          end if 
        end if 

        ! Density optimization 
        if ( this%databaseOptimization ) then
          ! Initialize database if not allocated
          if ( .not. allocated( this%kernelDatabaseFlat ) ) then 
              call this%InitializeKernelDatabaseFlat( this%minHOverLambda(1), &
                                                      this%maxHOverLambda(1), &
                                                    this%deltaHOverLambda(1), &
                                                      this%logKernelDatabase  )
          end if
          ! Compute density
          call this%ComputeDensityOptimization(                              &
                  this%densityEstimateGrid,                                  &
                  nOptimizationLoops=localNOptimizationLoops,                &
                  exportOptimizationVariables=locExportOptimizationVariables,&
                  skipErrorConvergence=locSkipErrorConvergence ) 
          ! Drop database ?
          if ( .not. persistKDB ) then
              call this%DropKernelDatabase()
          end if
        else
          ! Brute force optimization
          call this%ComputeDensityOptimization(                              &
                  this%densityEstimateGrid,                                  &
                  nOptimizationLoops=localNOptimizationLoops,                &
                  exportOptimizationVariables=locExportOptimizationVariables,&
                  skipErrorConvergence=locSkipErrorConvergence ) 
        end if 

        ! Some corrections to relevant variables before writing to output files !

        if ( locComputeRawDensity ) then 
          ! Compute the rawDensityEstimate: histogram/binvolume
          this%histogram%counts = this%histogram%counts/this%histogram%binVolume
          !if (allocated( this%rawDensityEstimateGrid )) deallocate( this%rawDensityEstimateGrid )
          !this%rawDensityEstimateGrid = this%histogram%counts/this%histogram%binVolume
          !if ( locUnitVolume ) then  
          !  ! If unit volume, modify 
          !  this%rawDensityEstimateGrid = &
          !  this%rawDensityEstimateGrid*this%histogram%binVolume
          !end if
        end if 

        if ( locUnitVolume ) then  
            ! If unit volume, modify 
            this%densityEstimateGrid = &
            this%densityEstimateGrid*this%histogram%binVolume
        end if

        if ( locScalingFactor .ne. 0d0 ) then
            ! Apply scalingFactor to density
            this%densityEstimateGrid = this%densityEstimateGrid*locScalingFactor
        end if

        if ( locScaleHistogram ) then
            ! Apply histogramScalingFactor to histogram
            this%histogram%counts = this%histogram%counts*locHistogramScalingFactor
        end if 


        ! Write output files !
        if ( present( outputFileUnit ) .and. present( outputDataId ) .and. present( particleGroupId )) then
            call this%ExportDensityUnit( outputFileUnit, outputDataId, particleGroupId )
        else if ( present( outputFileName ) ) then  
            call this%ExportDensity( outputFileName )
        end if
      

        ! Done
        return


    end subroutine prComputeDensity 
    

    ! Density optimization
    subroutine prComputeDensityOptimization( this, densityEstimateGrid, nOptimizationLoops, &
                                          exportOptimizationVariables, skipErrorConvergence )
      !------------------------------------------------------------------------------
      ! Performs the optimization loop 
      ! 
      !   - Section 2.5 in Sole-Mari et al.(2019) 
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      implicit none
      class( GridProjectedKDEType ), target:: this
      doubleprecision, dimension(:,:,:), intent(inout) :: densityEstimateGrid

      ! kernels
      type( KernelMultiGaussianType )     :: kernel
      type( KernelMultiGaussianType )     :: kernelSigma
      type( KernelSecondDerivativeXType ) :: kernelSDX
      type( KernelSecondDerivativeYType ) :: kernelSDY
      type( KernelSecondDerivativeZType ) :: kernelSDZ

      ! Optimization loops
      integer, intent(in), optional :: nOptimizationLoops
      logical, intent(in), optional :: exportOptimizationVariables
      logical, intent(in), optional :: skipErrorConvergence
      integer                       :: nOptLoops

      ! Grid cells
      type( GridCellType ), dimension(:), pointer :: activeGridCells => null()
      type( GridCellType ), pointer :: gc => null()

      ! kernelMatrix pointer
      doubleprecision, dimension(:,:,:), pointer :: kernelMatrix => null()
      doubleprecision, dimension(:,:,:), allocatable, target :: transposedKernelMatrix

      ! Utils
      integer            :: n, m, nd
      integer            :: convergenceCount 
      integer            :: softConvergenceCount
      integer            :: zeroDensityCount    
      character(len=500) :: varsOutputFileName
      character(len=500) :: errorOutputFileName
      character(len=20)  :: loopId
      character(len=20)  :: auxChar
      logical            :: exportVariables, skipErrorBreak
      logical            :: exportLoopError
      integer            :: errorOutputUnit

      ! Optimization error monitoring 
      doubleprecision :: errorRMSE
      doubleprecision :: errorRMSEOld
      doubleprecision :: errorALMISEProxy 
      doubleprecision :: errorALMISEProxyOld
      doubleprecision, dimension(:), allocatable     :: squareDensityDiff
      doubleprecision, dimension(:), allocatable     :: rawDensity
      doubleprecision, dimension(:), allocatable     :: errorMetricArray
      doubleprecision, dimension(:), allocatable     :: relativeDensityChange
      doubleprecision, dimension(:), allocatable     :: relativeSmoothingChange
      doubleprecision, dimension(:), allocatable     :: relativeRoughnessChange
      doubleprecision, dimension(:), allocatable     :: kernelSmoothingScaleOld
      doubleprecision, dimension(:,:), allocatable   :: kernelSmoothingOld
      doubleprecision, dimension(:), allocatable     :: densityEstimateArrayOld
      doubleprecision, dimension(:), allocatable     :: netRoughnessArrayOld
      doubleprecision, dimension(:,:,:), allocatable :: densityGridOld
      doubleprecision :: errorMetric
      doubleprecision :: errorMetricOld
      doubleprecision :: errorMetricSmoothing
      doubleprecision :: errorMetricSmoothingOld
      doubleprecision :: errorMetricConvergence
      integer         :: nDensityConvergence
      integer         :: nDensityConvergenceOld
      integer         :: nRoughnessConvergence
      integer         :: nRoughnessConvergenceOld
      integer         :: nSmoothingConvergence
      integer         :: nSmoothingConvergenceOld
      doubleprecision :: nFractionDensity
      doubleprecision :: nFractionRoughness
      doubleprecision :: nFractionSmoothing
      !character(len=16) :: helpChar
      !------------------------------------------------------------------------------

      ! Initialize vars
      convergenceCount = 0
      softConvergenceCount = 0
      zeroDensityCount = 0
      exportVariables  = .false.
      skipErrorBreak   = .false.
      exportLoopError  = .false.
      errorOutputUnit  = 999 
      
      ! Pointers to null
      gc => null()
      kernelMatrix => null()

      ! Allocate arrays according to nComputebins
      call prAllocateArrays( this%nComputeBins,      &
                             kernelSmoothing,        &
                             kernelSmoothingScale,   &
                             kernelSmoothingShape,   &
                             kernelSigmaSupport,     &
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
      if ( allocated(relativeDensityChange) ) deallocate(relativeDensityChange) 
      allocate( relativeDensityChange(this%nComputeBins) )
      if ( allocated(densityEstimateArrayOld) ) deallocate(densityEstimateArrayOld) 
      allocate( densityEstimateArrayOld(this%nComputeBins) )
      if ( allocated(errorMetricArray) ) deallocate(errorMetricArray) 
      allocate( errorMetricArray(this%nComputeBins) )
      if ( allocated(kernelSmoothingOld) ) deallocate(kernelSmoothingOld) 
      allocate( kernelSmoothingOld(3,this%nComputeBins) )
      if ( allocated(netRoughnessArrayOld) ) deallocate(netRoughnessArrayOld) 
      allocate( netRoughnessArrayOld(this%nComputeBins) )
      if ( allocated(relativeRoughnessChange) ) deallocate(relativeRoughnessChange) 
      allocate( relativeRoughnessChange(this%nComputeBins) )
      if ( allocated(relativeSmoothingChange) ) deallocate(relativeSmoothingChange) 
      allocate( relativeSmoothingChange(this%nComputeBins) )

      errorMetricConvergence   = this%densityRelativeConvergence
      nDensityConvergence      = 0d0
      nDensityConvergenceOld   = 0d0
      nRoughnessConvergence    = 0d0
      nRoughnessConvergenceOld = 0d0
      ! Process arguments
      if ( present( nOptimizationLoops ) ) then 
        nOptLoops = nOptimizationLoops
      else 
        nOptLoops = this%nOptimizationLoops
      end if 

      if ( present( skipErrorConvergence ) ) then 
        skipErrorBreak = skipErrorConvergence
      end if 
      if ( present( exportOptimizationVariables ) ) then 
        exportVariables = exportOptimizationVariables
      end if 

      ! Initialize active grid cells
      ! and compute rawDensity
      rawDensity = 0d0
      !$omp parallel do schedule(dynamic,1) &
      !$omp private(gc) 
      do n = 1, this%nComputeBins
        gc => activeGridCellsMod(n)
        call gc%Initialize( this%computeBinIds( :, n ) )
        rawDensity(n) = this%histogram%counts(gc%id(1),gc%id(2),gc%id(3))
      end do
      !$omp end parallel do
      rawDensity = rawDensity/this%histogram%binVolume
      activeGridCells => activeGridCellsMod

      ! Initialize kernels
      call kernel%Initialize(      this%binSize, matrixRange=defaultKernelRange )
      call kernelSigma%Initialize( this%binSize, matrixRange=defaultKernelRange )
      call kernelSDX%Initialize(   this%binSize, matrixRange=defaultKernelSDRange )
      call kernelSDY%Initialize(   this%binSize, matrixRange=defaultKernelSDRange )
      call kernelSDZ%Initialize(   this%binSize, matrixRange=defaultKernelSDRange )

      ! Something to detect a new loop/second run

      ! Define initial smoothing array
      ! initialSmoothing or kernelSmoothing could
      ! be constructed from results of previous optimization
      ! AFTER THE FIRST TIME 
      !if( .not. this%firstRun ) then 
      !    kernelSmoothing = spread( this%averageKernelSmoothing, 2, this%nComputeBins )
      !else
          kernelSmoothing = spread( this%initialSmoothing, 2, this%nComputeBins )
      !end if
      call prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )
      kernelSigmaSupportScale = 3d0*kernelSmoothingScale
      kernelSigmaSupport      = spread( kernelSigmaSupportScale, 1, 3 )
      do nd =1, 3
        if ( this%dimensionMask(nd) .eq. 1 ) then 
          where ( kernelSmoothingScale .gt. 0d0 )
            kernelSmoothingShape(nd,:) = kernelSmoothing(nd,:)/kernelSmoothingScale
          end where
        else
          ! No smoothing in compressed dimension 
          kernelSmoothing(nd,:)      = 0
          kernelSmoothingShape(nd,:) = 0
          kernelSigmaSupport(nd,:)   = 0
        end if 
      end do


      ! Initialize nEstimate
      nEstimateGrid = 0d0
      nEstimateArray = 0d0
      !$omp parallel do schedule( dynamic, 1 )                    &
      !$omp default( none )                                       &
      !$omp shared( this )                                        &
      !$omp shared( activeGridCells )                             &
      !$omp shared( nEstimateGrid, nEstimateArray )               &
      !$omp shared( kernelSigmaSupport, kernelSigmaSupportScale ) &
      !$omp firstprivate( kernelSigma )                           &
      !$omp private( gc )            
      do n = 1, this%nComputeBins

        ! Assign gc pointer
        gc => activeGridCells( n )

        if (  kernelSigmaSupportScale( n ) .lt. 0d0 ) cycle ! yes ?

        ! Set kernel sigma
        call this%SetKernelSigma( gc, kernelSigma, kernelSigmaSupport( :, n ) )

        ! Compute estimate, using rawDensity
        ! Notice division by volume after the loop
        nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
          this%histogram%counts(&
              gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
              gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
              gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
          )*gc%kernelSigmaMatrix(&
              gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
              gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
              gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

        ! Assign into array     
        nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
      end do
      !$omp end parallel do
      nEstimateGrid  = nEstimateGrid/this%histogram%binVolume
      nEstimateArray = nEstimateArray/this%histogram%binVolume

      ! Initialize variables
      !curvatureBandwidth = kernelSmoothing
      curvatureBandwidth = kernelSigmaSupport
      kernelSmoothingOld = kernelSmoothing
      kernelSmoothingScaleOld = kernelSmoothingScale
      ! Initial roughness
      call this%ComputeNetRoughnessEstimate(activeGridCells, curvatureBandwidth, &
                           roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                          netRoughnessArray, kernelSigmaSupport, &
                                                 kernelSDX, kernelSDY, kernelSDZ )

      ! Optimal smoothing
      call this%ComputeOptimalSmoothingAndShape( nEstimateArray, netRoughnessArray, & 
                              roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                               kernelSmoothing, kernelSmoothingOld, & 
                                     kernelSmoothingScale, kernelSmoothingScaleOld, & 
                                                               kernelSmoothingShape )
      ! Update smoothing scale
      call prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )

      ! Initialize density grid
      densityEstimateGrid = 0d0
      densityEstimateArray = 0d0
      !$omp parallel do schedule( dynamic, 1 )  &
      !$omp default( none )                     &
      !$omp shared( this )                      &
      !$omp shared( activeGridCells )           & 
      !$omp shared( kernelSmoothing )           & 
      !$omp reduction( +: densityEstimateGrid ) & 
      !$omp firstprivate( kernel )              & 
      !$omp private( gc )                        
      do n = 1, this%nComputeBins
          
          ! Assign gc pointer 
          gc => activeGridCells(n)

          if ( any( kernelSmoothing( :, n ) .lt. 0d0 ) ) cycle

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
            ) + this%histogram%counts(                 &
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

      ! Error monitoring
      squareDensityDiff = (densityEstimateArray - rawDensity)**2
      errorRMSE         = sqrt(sum( squareDensityDiff )/this%nComputeBins)

      ! Initialize error metric 
      errorMetricArray = 0d0
      where ( kernelSmoothingScale .ne. 0d0 ) 
        errorMetricArray = nEstimateArray/( (kernelSmoothingScale**nDim)*(4d0*pi)**(0.5*nDim)) + &
        0.25*netRoughnessArray*kernelSmoothingScale**4d0
      end where
      errorALMISEProxy = sqrt(sum(errorMetricArray**2)/this%nComputeBins)

      ! Initialize smoothing error trackers
      relativeSmoothingChange = 0d0
      where (kernelSmoothingScaleOld .ne. 0d0 ) 
        relativeSmoothingChange = abs(kernelSmoothingScale - & 
             kernelSmoothingScaleOld)/kernelSmoothingScaleOld
      end where
      nSmoothingConvergence = 0
      do n = 1, this%nComputeBins
        if ( relativeSmoothingChange(n) < errorMetricConvergence ) then 
            nSmoothingConvergence = nSmoothingConvergence + 1
        end if
      end do
      errorMetricSmoothing = sqrt(sum(relativeSmoothingChange**2)/this%nComputeBins)
      nFractionSmoothing = real(nSmoothingConvergence)/real(this%nComputeBins)


      ! Write error variables to output file
      if ( (exportLoopError) ) then
        ! Name of the loop-errors file
        write( unit=auxChar, fmt=* )999
        write( unit=errorOutputFileName, fmt='(a)' )trim(adjustl(this%outputFileName))//trim(adjustl(auxChar))
        ! Open output file name
        open( errorOutputUnit, file=errorOutputFileName, status='replace' )
        ! Write initial records
        call prWriteErrorMetricsRecord( this, errorOutputUnit, 0, &
                                  0d0, 0d0, errorMetricSmoothing, &  
                                    0d0, 0d0, nFractionSmoothing, & 
                                 errorRMSE, errorALMISEProxy, 0d0 )
      end if 


      ! Initialize old error trackers
      errorALMISEProxyOld     = errorALMISEProxy
      errorRMSEOld            = errorRMSE
      errorMetricOld          = 1d10 ! something big
      errorMetricSmoothingOld = errorMetricSmoothing
      densityEstimateArrayOld = densityEstimateArray
      kernelSmoothingScaleOld = kernelSmoothingScale
      kernelSmoothingOld      = kernelSmoothing
      netRoughnessArrayOld    = netRoughnessArray
      densityGridOld          = densityEstimateGrid


      ! Optimization loop !
      do m = 1, nOptLoops

          ! nEstimate
          nEstimateGrid = 0d0 
          nEstimateArray = 0d0
          !$omp parallel do schedule( dynamic, 1 )                    &
          !$omp default( none )                                       &
          !$omp shared( this )                                        &
          !$omp shared( activeGridCells )                             &
          !$omp shared( densityEstimateGrid )                         &
          !$omp shared( nEstimateGrid, nEstimateArray )               &
          !$omp shared( kernelSigmaSupport, kernelSigmaSupportScale ) &
          !$omp firstprivate( kernelSigma )                           &
          !$omp private( gc )            
          do n = 1, this%nComputeBins

            ! Assign gc pointer
            gc => activeGridCells( n )

            if (  kernelSigmaSupportScale( n ) .lt. 0d0 ) cycle ! yes ?

            ! Set kernel sigma
            call this%SetKernelSigma( gc, kernelSigma, kernelSigmaSupport( :, n ) )

            ! Compute estimate
            nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
              densityEstimateGrid(&
                  gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                  gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                  gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
              )*gc%kernelSigmaMatrix(&
                  gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                  gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                  gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

            ! Assign into array     
            nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
          end do
          !$omp end parallel do

          ! Export optimization variables 
          if ( (exportVariables) .and. (m.eq.1) ) then
            ! Export initial optimization variables
            write( unit=loopId, fmt=* )m-1
            write( unit=varsOutputFileName, fmt='(a)' )trim(adjustl(this%outputFileName))//trim(adjustl(loopId))
            call prExportOptimizationVariablesExtended( this, varsOutputFileName, & 
              densityEstimateArray, kernelSmoothing, kernelSmoothingScale,kernelSmoothingShape,  & 
              kernelSigmaSupportScale, &
              curvatureBandwidth, nEstimateArray, roughnessXXArray, &
              roughnessYYArray, roughnessZZArray, netRoughnessArray )
          end if 
       
          ! Compute support scale 
          call this%ComputeSupportScale( kernelSmoothingScale, densityEstimateArray, & 
                                             nEstimateArray, kernelSigmaSupportScale )

          ! Spread the support scale as isotropic 
          ! And deactivate compressed dimension
          kernelSigmaSupport = spread( kernelSigmaSupportScale, 1, 3 )
          do nd =1, 3
            if ( this%dimensionMask(nd) .eq. 0 ) then 
              ! No smoothing in compressed dimension 
              kernelSigmaSupport(nd,:)   = 0d0
            end if 
          end do

          ! Update nEstimate
          nEstimateGrid  = 0d0
          nEstimateArray = 0d0
          !$omp parallel do schedule( dynamic, 1 )                    &
          !$omp default( none )                                       &
          !$omp shared( this )                                        &
          !$omp shared( activeGridCells )                             &
          !$omp shared( densityEstimateGrid )                         &
          !$omp shared( nEstimateGrid, nEstimateArray )               &
          !$omp shared( kernelSigmaSupport, kernelSigmaSupportScale ) &
          !$omp firstprivate( kernelSigma )                           &
          !$omp private( gc )
          do n = 1, this%nComputeBins

            ! Assign gc pointer 
            gc => activeGridCells(n)

            ! Set kernel sigma
            call this%SetKernelSigma( gc, kernelSigma, kernelSigmaSupport( :, n ) )

            ! Compute estimate
            nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) ) = sum(&
              densityEstimateGrid(&
                  gc%kernelSigmaXGSpan(1):gc%kernelSigmaXGSpan(2), &
                  gc%kernelSigmaYGSpan(1):gc%kernelSigmaYGSpan(2), & 
                  gc%kernelSigmaZGSpan(1):gc%kernelSigmaZGSpan(2)  & 
              )*gc%kernelSigmaMatrix(&
                  gc%kernelSigmaXMSpan(1):gc%kernelSigmaXMSpan(2), &
                  gc%kernelSigmaYMSpan(1):gc%kernelSigmaYMSpan(2), & 
                  gc%kernelSigmaZMSpan(1):gc%kernelSigmaZMSpan(2)) )

            ! Assign into array     
            nEstimateArray( n ) = nEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
          end do
          !$omp end parallel do 

          ! Curvature bandwidths
          call this%ComputeCurvatureKernelBandwidth( densityEstimateArray, nEstimateArray, &
                              kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape, & 
                                               kernelSigmaSupportScale, curvatureBandwidth )

          ! Net roughness
          call this%ComputeNetRoughnessEstimate(activeGridCells, curvatureBandwidth, &
                               roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                              netRoughnessArray, kernelSigmaSupport, &
                                                     kernelSDX, kernelSDY, kernelSDZ )
   
          ! Optimal smoothing
          call this%ComputeOptimalSmoothingAndShape( nEstimateArray, netRoughnessArray, & 
                                  roughnessXXArray, roughnessYYArray, roughnessZZArray, &
                                                   kernelSmoothing, kernelSmoothingOld, & 
                                         kernelSmoothingScale, kernelSmoothingScaleOld, & 
                                                                   kernelSmoothingShape )
          ! Update smoothing scale
          call prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )

          ! Update density
          densityEstimateGrid = 0d0
          densityEstimateArray = 0d0
          !$omp parallel do schedule( dynamic, 1 )  &
          !$omp default( none )                     &
          !$omp shared( this )                      &
          !$omp shared( activeGridCells )           & 
          !$omp shared( kernelSmoothing )           & 
          !$omp reduction( +: densityEstimateGrid ) & 
          !$omp firstprivate( kernel )              & 
          !$omp private( gc )                        
          do n = 1, this%nComputeBins

            ! Any smoothing < 0 or NaN, skip
            if ( ( any( kernelSmoothing( :, n ) .lt. 0d0 ) ) .or.             &
                ( any( kernelSmoothing( :, n ) /= kernelSmoothing( :, n ) ) ) ) then
              cycle
            end if

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
              ) + this%histogram%counts(                 &
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

          ! A proxy to error: relative density change
          relativeDensityChange = 0d0
          where ( densityEstimateArrayOld .ne. 0d0 )
            relativeDensityChange = abs(densityEstimateArray/maxval(densityEstimateArray) -   & 
                                      densityEstimateArrayOld/maxval(densityEstimateArrayOld) & 
                                      )/(densityEstimateArrayOld/maxval(densityEstimateArrayOld) )
          end where
          errorMetric = sqrt( sum(relativeDensityChange**2)/this%nComputeBins )
          ! A proxy to error: relative roughness change 
          relativeRoughnessChange = 0d0
          where (netRoughnessArrayOld .ne. 0d0 ) 
            relativeRoughnessChange = abs(netRoughnessArray/maxval(netRoughnessArray) - &
                  netRoughnessArrayOld/maxval(netRoughnessArrayOld) & 
                  )/(netRoughnessArrayOld/maxval(netRoughnessArrayOld))
          end where
          ! A proxy to error: relative smoothing change
          relativeSmoothingChange = 0d0
          where (kernelSmoothingScaleOld .ne. 0d0 ) 
            relativeSmoothingChange = abs(kernelSmoothingScale - & 
                 kernelSmoothingScaleOld)/kernelSmoothingScaleOld
          end where
          errorMetricSmoothing = sqrt(sum(relativeSmoothingChange**2)/this%nComputeBins)

          ! Counters
          nDensityConvergence = 0
          nRoughnessConvergence = 0
          nSmoothingConvergence = 0
          !$omp parallel do schedule(static)      &    
          !$omp default(none)                     &
          !$omp shared(this)                      &
          !$omp shared(errorMetricConvergence)    &
          !$omp shared(relativeDensityChange)     &
          !$omp shared(relativeRoughnessChange)   &
          !$omp shared(relativeSmoothingChange)   &
          !$omp reduction(+:nDensityConvergence)  &
          !$omp reduction(+:nRoughnessConvergence)&
          !$omp reduction(+:nSmoothingConvergence)
          do n = 1, this%nComputeBins
            if ( relativeDensityChange(n) < errorMetricConvergence ) then 
              nDensityConvergence = nDensityConvergence + 1
            end if
            if ( relativeRoughnessChange(n) < errorMetricConvergence ) then 
              nRoughnessConvergence = nRoughnessConvergence + 1
            end if
            if ( relativeSmoothingChange(n) < errorMetricConvergence ) then 
              nSmoothingConvergence = nSmoothingConvergence + 1
            end if
          end do
          !$omp end parallel do
          nFractionDensity   = real(nDensityConvergence)/real(this%nComputeBins)
          nFractionRoughness = real(nRoughnessConvergence)/real(this%nComputeBins)
          nFractionSmoothing = real(nSmoothingConvergence)/real(this%nComputeBins)

          ! A proxy to error: an estimate of ALIMISE
          errorMetricArray = 0d0
          where ( kernelSmoothingScale .ne. 0d0 ) 
              errorMetricArray = nEstimateArray/( (kernelSmoothingScale**nDim)*(4d0*pi)**(0.5*nDim)) + &
              0.25*netRoughnessArray*kernelSmoothingScale**4d0
          end where
          errorALMISEProxy = sqrt(sum(errorMetricArray**2)/this%nComputeBins)

          ! A proxy to error: RMSE versus histogram density
          squareDensityDiff = (densityEstimateArray - rawDensity)**2
          errorRMSE         = sqrt(sum( squareDensityDiff )/this%nComputeBins)

          ! Error analysis:
          if ( .not. skipErrorBreak ) then
            ! NOTE: a new criteria could consider the 
            ! densityGrid in order to include in the error 
            ! estimate those cells without particles/mass.
            if (  ( errorMetric .lt. errorMetricConvergence ) ) then
              ! Criteria:
              ! Break optimization loop if 
              ! relative density change lower than 
              ! a given convergence
              this%averageKernelSmoothing = sum( kernelSmoothing, dim=2 )/this%nComputeBIns
              !print *, '!! DENSITY CONVERGENCE !!'
              if ( this%reportToOutUnit ) then 
              write( this%outFileUnit, '(A,es13.4e2)' ) '    - Density convergence ', errorMetric
              end if 
              ! Break
              exit
            end if 
            if ( ( errorMetricSmoothing .lt. errorMetricConvergence ) ) then 
              ! Criteria:
              ! Break optimization loop if 
              ! relative smoothing change lower than 
              ! a given convergence
              this%averageKernelSmoothing = sum( kernelSmoothing, dim=2 )/this%nComputeBIns
              !print *, '!! SMOOTHING CONVERGENCE !!'
              if ( this%reportToOutUnit ) then 
              write( this%outFileUnit, '(A,es13.4e2)' ) '    - Bandwidth convergence ', errorMetricSmoothing
              end if
              ! Break
              exit
            end if 
            if ( (errorALMISEProxy .gt. errorALMISEProxyOld ) .and. & 
                    (errorRMSE .gt. errorRMSEOld ) .and.               &
                    (errorMetric .lt. defaultRelaxedDensityRelativeConvergence) ) then 
              ! Criteria
              ! If both RMSE and ALMISE are increasing
              ! and density convergence is below the relaxed limit, 
              ! return previous density and leave
              densityEstimateGrid = densityGridOld
              ! Transfer grid density to array
              do n = 1, this%nComputeBins
                ! Assign gc pointer 
                gc => activeGridCells(n)
                densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
              end do
              this%averageKernelSmoothing = sum( kernelSmoothing, dim=2 )/this%nComputeBIns
              !print *, '!! INCREASED RMSE/ALMISE AND SOFT CONVERGENCE !!'
              if ( this%reportToOutUnit ) then 
              write( this%outFileUnit, '(A,es13.4e2)' ) '    - Relaxed density convergence ', errorMetricOld
              end if 
              ! Break
              exit
            end if 
            !if ( & 
            !  (errorRMSE .lt. errorRMSEOld) .and. ( errorALMISEProxy .gt. errorALMISEProxyOld ) .and. & 
            !  (errorMetric .lt. defaultRelaxedDensityRelativeConvergence )   ) then
            !  ! Criteria
            !  ! If the RMSE versus the histogram decreases, 
            !  ! but the ALMISE increases it is probably and indication 
            !  ! of high resolution in the particle model. 
            !  ! This has been observed for high number of particles 
            !  ! and error indicators oscillate.
            !  if ( this%reportToOutUnit ) then 
            !  write( this%outFileUnit, '(A,es10.4e2)' ) '    - Relaxed density convergence ', errorMetric
            !  end if 
            !  !print *, '!! INCREASED ALMISE DECREASE RMSE AND SOFT CONVERGENCE !!'
            !  ! Break
            !  exit
            !end if 
            !if ( (errorMetric .gt. errorMetricOld) .and. & 
            !  (errorMetric .lt. defaultRelaxedDensityRelativeConvergence) ) then
            !  ! Criteria
            !  ! If the relative change in density increased with 
            !  ! respect to previous value, return previous density and leave
            !  ! This could be complemented with maximum density analysis, 
            !  ! requires mass. 
            !  densityEstimateGrid = densityGridOld
            !  ! Transfer grid density to array
            !  do n = 1, this%nComputeBins
            !    gc => activeGridCells(n)
            !    densityEstimateArray( n ) = densityEstimateGrid( gc%id(1), gc%id(2), gc%id(3) )
            !  end do
            !  this%averageKernelSmoothing = sum( kernelSmoothing, dim=2 )/this%nComputeBIns
            !  !print *, '!! INCREASED RELATIVE DENSITY CHANGE AND SOFT CONVERGENCE !!'
            !  if ( this%reportToOutUnit ) then 
            !  write( this%outFileUnit, '(A,es10.4e2)' ) '    - Relaxed density convergence ', errorMetricOld
            !  end if 
            !  ! Break
            !  exit
            !end if
            if ( nFractionDensity .gt. 0.98 ) then 
              ! Criteria
              ! If the fraction of cells that is presenting changes in 
              ! density below the convergence criteria is close to the total 
              ! number of cells, exit.
              !print *, '!! NFRACTIONDENSITY !!'
              if ( this%reportToOutUnit ) then 
              write( this%outFileUnit, '(A)' ) '    - nFractionDensity '
              end if 
              ! Break
              exit
            end if  
            if ( nFractionSmoothing .gt. 0.98 ) then 
              ! Criteria
              ! If the fraction of cells that is presenting changes in 
              ! smoothing below the convergence criteria is close to the total 
              ! number of cells, exit.
              !print *, '!! NFRACTIONSMOOTHING !!'
              if ( this%reportToOutUnit ) then 
              write( this%outFileUnit, '(A)' ) '    - nFractionSmoothing '
              end if 
              ! Break
              exit
            end if  
          end if


          if ( exportLoopError ) then 
            ! Write initial records
            call prWriteErrorMetricsRecord( this, errorOutputUnit, m,      &
                  sqrt(sum(relativeDensityChange**2)/this%nComputeBins),   &  
                  sqrt(sum(relativeRoughnessChange**2)/this%nComputeBins), &  
                                                     errorMetricSmoothing, &  
                 nFractionDensity, nFractionRoughness, nFractionSmoothing, & 
                                        errorRMSE, errorALMISEProxy, 0d0 )
          end if 

          ! Continue to next loop !
          errorALMISEProxyOld      = errorALMISEProxy
          errorRMSEOld             = errorRMSE
          errorMetricOld           = errorMetric
          kernelSmoothingOld       = kernelSmoothing
          netRoughnessArrayOld     = netRoughnessArray
          densityEstimateArrayOld  = densityEstimateArray
          densityGridOld           = densityEstimateGrid
          nDensityConvergenceOld   = nDensityConvergence
          nRoughnessConvergenceOld = nRoughnessConvergence
          nSmoothingConvergenceOld = nSmoothingConvergence
          call prComputeKernelSmoothingScale( this, kernelSmoothingOld, kernelSmoothingScaleOld )

          ! Export optimization variables
          if ( exportVariables ) then
            write( unit=loopId, fmt=* )m
            write( unit=varsOutputFileName, fmt='(a)' )trim(adjustl(this%outputFileName))//trim(adjustl(loopId))
            !call prExportOptimizationVariables( this, varsOutputFileName, & 
            !    densityEstimateArray, kernelSmoothing, kernelSigmaSupportScale, &
            !    curvatureBandwidth, nEstimateArray, netRoughnessArray )
            call prExportOptimizationVariablesExtendedError( this, varsOutputFileName, & 
                densityEstimateArray, kernelSmoothing, kernelSmoothingScale,kernelSmoothingShape,  & 
                kernelSigmaSupportScale, &
                curvatureBandwidth, nEstimateArray, roughnessXXArray, &
                roughnessYYArray, roughnessZZArray, netRoughnessArray, &
                         relativeDensityChange, relativeRoughnessChange )

          end if 

      end do
      ! End optimization loop ! 

      if ( ((m-1).eq.nOptLoops).and.this%reportToOutUnit ) then 
        write( this%outFileUnit, '(A)' ) '    - Max loops '
      end if

      ! Finalize error output unit
      if ( exportLoopError ) then
        close( errorOutputUnit )
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
      
      ! It run
      this%firstRun = .false.


    end subroutine prComputeDensityOptimization



    subroutine prComputeKernelSmoothingScale( this, kernelSmoothing, kernelSmoothingScale )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(:,:), intent(in)   :: kernelSmoothing
        doubleprecision, dimension(:), intent(inout)  :: kernelSmoothingScale
        integer :: nd
        integer, dimension(:), allocatable :: dimCorrected
        !------------------------------------------------------------------------------
        
        allocate( dimCorrected(this%nComputeBins) )
        dimCorrected = nDim
        kernelSmoothingScale = 1d0
        do nd = 1, 3
            if ( this%dimensionMask(nd) .eq. 1 ) then
                where( kernelSmoothing(nd,:) .ne. 0d0 )  
                    kernelSmoothingScale = kernelSmoothingScale*kernelSmoothing(nd,:) 
                elsewhere
                    dimCorrected = dimCorrected - 1  
                end where
            end if 
        end do
        where (dimCorrected .ne. 0 )
            kernelSmoothingScale = ( kernelSmoothingScale )**( 1d0/nDim )
        elsewhere
            kernelSmoothingScale = 0d0
        end where
        deallocate( dimCorrected )

        ! Done
        return

    end subroutine prComputeKernelSmoothingScale



    ! Optimization loop functions
    subroutine prComputeSupportScale( this, kernelSmoothingScale, densityEstimate, &
                                              nEstimate, kernelSigmaSupportScale )
      !------------------------------------------------------------------------------
      ! Smoothing scale for the support kernel
      ! 
      !   - Eq. 23 in Sole-Mari et al. (2019) 
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      implicit none
      class( GridProjectedKDEType ) :: this
      doubleprecision, dimension(:), intent(in)    :: kernelSmoothingScale
      doubleprecision, dimension(:), intent(in)    :: densityEstimate
      doubleprecision, dimension(:), intent(in)    :: nEstimate
      doubleprecision, dimension(:), intent(inout) :: kernelSigmaSupportScale
      !------------------------------------------------------------------------------

      ! Reset array
      kernelSigmaSupportScale = 0d0

      ! Compute
      where ( densityEstimate .ne. 0d0 ) 
        kernelSigmaSupportScale = nEstimate**(0.5)*kernelSmoothingScale**( 1d0 + 0.25*nDim )/&
                       ( ( 4d0*densityEstimate )**0.25 )*this%supportDimensionConstant
      end where

      ! Force minimun sigma scale in relation to kernel scale
      where ( ( kernelSigmaSupportScale .gt. 3d0*kernelSmoothingScale ) .or. ( &
        densityEstimate .eq. 0d0 )  ) 
        kernelSigmaSupportScale = 3d0*kernelSmoothingScale
      end where 

      !! Limit the maximum/minimum value
      !do nd =1, 3
      !    if ( this%dimensionMask(nd) .eq. 1 ) then 
      !        where ( kernelSigmaSupportScale/this%binSize(nd) .gt. 3d0*this%maxHOverLambda(nd) )
      !            kernelSigmaSupportScale = 3d0*this%binSize(nd)*this%maxHOverLambda(nd) 
      !        end where
      !        !where ( kernelSigmaSupportScale/this%binSize(nd) .lt. this%minHOverLambda(nd) )
      !        !    kernelSigmaSupportScale = this%binSize(nd)*this%minHOverLambda(nd) 
      !        !end where
      !    end if 
      !end do
      !where ( kernelSigmaSupportScale .gt. 3d0*kernelSmoothingScale )
      !    kernelSigmaSupportScale =  3d0*kernelSmoothingScale
      !end where
      !where ( kernelSigmaSupportScale/this%binSize(1) .gt. this%maxHOverLambda(1) )
      !    kernelSigmaSupportScale = this%binSize(1)*this%maxHOverLambda(1) 
      !end where

      ! Done
      return

    end subroutine prComputeSupportScale



    subroutine prComputeCurvatureKernelBandwidth( this, densityEstimate, nEstimate, &
                       kernelSmoothing, kernelSmoothingScale, kernelSmoothingShape, &
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
      doubleprecision, dimension(:),   intent(in)    :: nEstimate
      doubleprecision, dimension(:),   intent(in)    :: densityEstimate
      doubleprecision, dimension(:,:), intent(in)    :: kernelSmoothing
      doubleprecision, dimension(:),   intent(in)    :: kernelSmoothingScale
      doubleprecision, dimension(:,:), intent(in)    :: kernelSmoothingShape
      doubleprecision, dimension(:),   intent(in)    :: kernelSigmaSupportScale
      ! out 
      doubleprecision, dimension(:,:), intent(inout) :: curvatureBandwidth
      ! local 
      doubleprecision, dimension(:),   allocatable   :: nVirtualPowerBeta
      doubleprecision, dimension(:,:), allocatable   :: shapeTerm
      doubleprecision, dimension(:)  , allocatable   :: shapeTermSum
      integer, dimension(3)                          :: shapeTermNums = 1
      integer :: n, nActiveBins, nd
      !----------------------------------------------------------------------------

      ! Allocate local arrays
      nActiveBins = size( nEstimate ) ! Maybe removed
      allocate( shapeTerm( 3, nActiveBins  ) )
      allocate( shapeTermSum(  nActiveBins ) )
      allocate( nVirtualPowerBeta( nActiveBins ) )


      ! Compute virtual particle cloud size
      nVirtualPowerBeta = 0d0
      where ( densityEstimate .gt. 0d0 )  
        nVirtualPowerBeta = ( ( sqrtEightPi*kernelSigmaSupportScale )**nDim*&
            nEstimate**2d0/densityEstimate )**this%betaDimensionConstant
      end where


      ! Compute shape dependent terms
      curvatureBandwidth = 0d0
      shapeTerm = 0d0
      do nd = 1, 3
        if ( this%dimensionMask(nd) .eq. 1 ) then 
          shapeTermNums     = 1
          shapeTermNums(nd) = 5
          ! Compute sum for shape term
          shapeTermSum = 0d0
          do n =1,3
              if ( this%dimensionMask(n) .eq. 1 ) then
                  where( kernelSmoothingShape(n,:) .ne. 0d0 ) 
                      shapeTermSum = shapeTermSum + shapeTermNums(n)/( kernelSmoothingShape(n,:)**2 ) 
                  end where
              end if
          end do 
          where( kernelSmoothingShape(nd,:) .ne. 0d0 ) 
              shapeTerm( nd, : ) = (                                           &
                  ( 1d0/( nDim + 4d0 )/( kernelSmoothingShape( nd, : )**4d0 ) )*   &
                      (                                                        &
                          shapeTermSum                                         &
                      )                                                        &
                  )**( -1d0/( nDim + 6d0 ) )
          end where
          curvatureBandwidth( nd, : ) = &
              this%alphaDimensionConstant*nVirtualPowerBeta*shapeTerm( nd, : )*kernelSmoothingScale

          ! Limit the maximum/minimum value
          !where ( curvatureBandwidth( nd, : )/this%binSize(nd) .gt. 3d0*this%maxHOverLambda(nd) )
          !    curvatureBandwidth( nd, : ) = 3d0*this%binSize(nd)*this%maxHOverLambda(nd) 
          !end where
          !where ( curvatureBandwidth( nd, : )/this%binSize(nd) .lt. this%minHOverLambda(nd) )
          !    curvatureBandwidth( nd, : ) = & 
          !        this%alphaDimensionConstant*nVirtualPowerBeta*shapeTerm( nd, : )*this%binSize(nd)*this%minHOverLambda(nd) 
          !end where

        end if
      end do 

      deallocate( shapeTerm )
      deallocate( nVirtualPowerBeta )
  
      return


    end subroutine prComputeCurvatureKernelBandwidth



    subroutine prComputeOptimalSmoothingAndShape( this, nEstimate, netRoughness, &
                        roughnessXXActive, roughnessYYActive, roughnessZZActive, &
                                            kernelSmoothing, kernelSmoothingOld, & 
                                  kernelSmoothingScale, kernelSmoothingScaleOld, & 
                                                           kernelSmoothingShape  )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType) :: this
        doubleprecision, dimension(:), intent(in)      :: nEstimate 
        doubleprecision, dimension(:), intent(in)      :: netRoughness 
        doubleprecision, dimension(:), intent(in)      :: roughnessXXActive 
        doubleprecision, dimension(:), intent(in)      :: roughnessYYActive 
        doubleprecision, dimension(:), intent(in)      :: roughnessZZActive 
        doubleprecision, dimension(:,:), intent(inout) :: kernelSmoothing
        doubleprecision, dimension(:,:), intent(inout) :: kernelSmoothingOld
        doubleprecision, dimension(:),   intent(inout) :: kernelSmoothingScale
        doubleprecision, dimension(:),   intent(inout) :: kernelSmoothingScaleOld
        doubleprecision, dimension(:,:), intent(inout) :: kernelSmoothingShape
        doubleprecision, dimension(:), allocatable     :: roughnessScale
        integer :: nd
        !------------------------------------------------------------------------------

        allocate( roughnessScale( this%nComputeBins ) )
      
        ! Compute smoothing scale. 
        ! For the cases of really low roghness, use the limit
        ! value to keep bound the smoothing scale
        kernelSmoothingScale = 0d0
        where (abs( netRoughness ) .gt. this%minLimitRoughness )
            kernelSmoothingScale = ( nDim*nEstimate/( ( 4*pi )**( 0.5*nDim )*netRoughness ) )**( 1d0/( nDim + 4d0 ) )
        elsewhere
            ! Estimate a scale based on minLimitRoughness
            kernelSmoothingScale = ( nDim*nEstimate/( ( 4*pi )**( 0.5*nDim )*this%minLimitRoughness ) )**( 1d0/( nDim + 4d0 ) )
        end where

        ! Bounded
        ! FIX ! 
        where( kernelSmoothingScale .gt. this%maxHOverLambda(1)*this%binSize(1) )
            kernelSmoothingScale = this%maxHOverLambda(1)*this%binSize(1)
        end where
         

        ! Compute roughness scale,
        ! Even better: set as netRoughness
        roughnessScale = netRoughness
        ! Limit roughness max roughness scale ?
        where ( roughnessScale .gt. this%maxLimitRoughness ) 
            roughnessScale = this%maxLimitRoughness
        end where 
        !roughnessScale = 1d0
        !do nd=1,3
        !    if ( this%dimensionMask(nd) .eq. 1 ) then
        !        select case (nd) 
        !            case (1)  
        !                roughnessScale = roughnessScale*roughnessXXActive
        !            case (2) 
        !                roughnessScale = roughnessScale*roughnessYYActive
        !            case (3) 
        !                roughnessScale = roughnessScale*roughnessZZActive
        !        end select
        !    end if
        !end do
        !roughnessScale = roughnessScale**( 1d0/nDim )

        ! The fact that roughnessScale can be zero
        ! indicate that one or more kernel dimensions 
        ! where compressed

        ! Compute shape factors and kernelSmoothing
        kernelSmoothing = 0d0
        kernelSmoothingShape = 0d0
        do nd=1,3
            if ( this%dimensionMask(nd) .eq. 1 ) then
                ! For anisotropic kernels
                select case (nd) 
                    case (1)   
                        where (  abs(netRoughness) .gt. this%minLimitRoughness ) 
                            kernelSmoothingShape( nd, : ) = ( roughnessScale/roughnessXXActive )**( 0.25 )
                        end where
                    case (2) 
                        where (  abs(netRoughness) .gt. this%minLimitRoughness ) 
                            kernelSmoothingShape( nd, : ) = ( roughnessScale/roughnessYYActive )**( 0.25 )
                        end where
                    case (3) 
                        where (  abs(netRoughness) .gt. this%minLimitRoughness ) 
                            kernelSmoothingShape( nd, : ) = ( roughnessScale/roughnessZZActive )**( 0.25 )
                        end where
                end select


                ! Minimum roughness sets transition to isotropic        
                !where ( abs(netRoughness) .lt. this%minLimitRoughness ) 
                !    kernelSmoothingShape( nd, : ) = 1d0
                !end where
                where( kernelSmoothingShape(nd,:) .lt. this%minKernelShape )
                    kernelSmoothingShape(nd,:) = this%minKernelShape
                end where
                where( kernelSmoothingShape(nd,:) .gt. this%maxKernelShape )
                    kernelSmoothingShape(nd,:) = this%maxKernelShape
                end where
                kernelSmoothing( nd, : ) = kernelSmoothingShape( nd, : )*kernelSmoothingScale


                ! Limit the growth for maximum/minimum values
                ! What if previous was zero ?
                where ( ( kernelSmoothing(nd,:) .gt. (1d0 + this%maxSmoothingGrowth)*kernelSmoothingOld(nd,:) ) .and. &
                        ( kernelSmoothingOld(nd,:) .ne. 0d0 ) )
                    kernelSmoothing(nd,:) = (1d0 + this%maxSmoothingGrowth)*kernelSmoothingOld(nd,:)
                end where
                !where ( (kernelSmoothing(nd,:) .lt. kernelSmoothingOld(nd,:)/( 1d0 + this%maxSmoothingGrowth ) ) &
                !        .and. ( kernelSmoothing(nd,:) .ne. 0d0 ) )
                !    kernelSmoothing( nd, : ) = kernelSmoothingOld(nd,:)/( 1d0 + this%maxSmoothingGrowth )
                !end where

                ! Limit the maximum/minimum value
                where ( kernelSmoothing( nd, : )/this%binSize(nd) .gt. this%maxHOverLambda(nd) )
                    kernelSmoothing( nd, : ) = this%binSize(nd)*this%maxHOverLambda(nd) 
                end where
                where ( kernelSmoothing( nd, : )/this%binSize(nd) .lt. this%minHOverLambda(nd) )
                    kernelSmoothing( nd, : ) = this%binSize(nd)*this%minHOverLambda(nd) 
                end where
                !where ( kernelSmoothing( nd,:) .lt. this%binSize(nd)*this%minHOverLambda(nd) ) 
                !where ( abs(roughnessScale) .lt. 1d-2 ) 
                !where ( abs(netRoughness) .lt. this%minLimitRoughness ) 
                !    !kernelSmoothing(nd,:) = kernelSmoothingOld(nd,:)
                !    kernelSmoothing(nd,:) = 0d0 
                !end where

                where ( kernelSmoothingScale .gt. 0d0 ) 
                    kernelSmoothingShape(nd,:) = kernelSmoothing(nd,:)/kernelSmoothingScale
                end where  

                !where ( abs(netRoughness) .lt. this%minLimitRoughness ) 
                !    kernelSmoothingShape(nd,:) = 0d0
                !end where


            end if
        end do
      

        ! Should deallocate ?
        deallocate( roughnessScale )

      
        return


    end subroutine prComputeOptimalSmoothingAndShape



    subroutine prComputeKernelDatabaseFlatIndexesLog( this, smoothing, flatDBIndexes, transposeKernel )
      !------------------------------------------------------------------------------
      ! 
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      implicit none
      class( GridProjectedKDEType ) :: this
      doubleprecision, dimension(3), intent(in) :: smoothing
      integer, dimension(2), intent(inout)      :: flatDBIndexes
      logical, intent(inout)                    :: transposeKernel
      ! local 
      integer, dimension(3) :: indexes 
      integer :: nd
      ! A more robust function would be good 
      !------------------------------------------------------------------------------

      ! Initialize indexes 
      indexes(:) = 1
      transposeKernel = .false.

      ! Compute index value if required 
      ! because active dimension
      do nd = 1, 3
        if ( smoothing( nd ) .le. 0d0 ) cycle
        indexes(nd) = min(&
          max(&
            floor(&
              log( smoothing(nd)/this%binSize(nd)/this%minHOverLambda(nd) )/this%deltaHOverLambda(nd)&
            ) + 1, 1 &
          ), &
        this%nDeltaHOverLambda(nd)  )
      end do 
      
      ! 1D
      if ( nDim .eq. 1 ) then 
        flatDBIndexes(1) = maxval(indexes)
        flatDBIndexes(2) = 1 
        ! Done 
        return
      end if 

      ! 2D
      ! Will work properly as long nDeltaHOverLambda
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
      ! Will work properly as long nDeltaHOverLambda
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


      return


    end subroutine prComputeKernelDatabaseFlatIndexesLog


    ! OUTDATED
    ! Kernel Database indexes, Flat
    subroutine prComputeKernelDatabaseFlatIndexesLinear( this, smoothing, flatDBIndexes, transposeKernel )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(3), intent(in) :: smoothing
        integer, dimension(2), intent(inout)      :: flatDBIndexes
        logical, intent(inout)                    :: transposeKernel
        ! local 
        integer, dimension(3) :: indexes 
        integer :: nd 
        !------------------------------------------------------------------------------
        
        indexes = 0

        do nd = 1, nDim
            indexes(nd) = min(&
                max(&
                    floor(&
                        (smoothing(nd)/this%binSize(nd) - this%minHOverLambda(nd))/this%deltaHOverLambda(nd)&
                    ) + 1, 1 &
                ), &
            this%nDeltaHOverLambda(nd)  )
        end do 
       

        ! 1D
        if ( nDim .eq. 1 ) then 
                flatDBIndexes(1) = maxval(indexes)
                flatDBIndexes(2) = 1 
                ! Done 
                return
        end if 

        ! Will work properly as long nDeltaHOverLambda
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


        return


    end subroutine prComputeKernelDatabaseFlatIndexesLinear

    ! DEPRECATION WARNING !
    ! Kernel Database indexes, 3D
    function prComputeKernelDatabaseIndexesLinear( this, smoothing ) result(indexes)
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(3), intent(in) :: smoothing
        integer, dimension(3) :: indexes 
        integer :: nd 
        !------------------------------------------------------------------------------


        indexes = 0

        do nd = 1, nDim
            indexes(nd) = min(&
                max(&
                    floor(&
                        (smoothing(nd)/this%binSize(nd) - this%minHOverLambda(nd))/this%deltaHOverLambda(nd)&
                    ) + 1, 1 &
                ), &
            this%nDeltaHOverLambda(nd)  )
        end do 


        return


    end function prComputeKernelDatabaseIndexesLinear
    ! DEPRECATION WARNING !
    function prComputeKernelDatabaseIndexesLog( this, smoothing ) result(indexes)
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridProjectedKDEType ) :: this
        doubleprecision, dimension(3), intent(in) :: smoothing
        integer, dimension(3) :: indexes 
        integer :: nd 
        !------------------------------------------------------------------------------

        ! Initialize indexes 
        indexes = 1

        ! Compute indexes where smoothing > 0d0
        do nd = 1, 3
            if ( smoothing( nd ) .le. 0d0 ) cycle
            indexes(nd) = min(&
                max(&
                    floor(&
                        log( smoothing(nd)/this%binSize(nd)/this%minHOverLambda(nd) )/this%deltaHOverLambda(nd)&
                    ) + 1, 1 &
                ), &
            this%nDeltaHOverLambda(nd)  )
        end do 


        return


    end function prComputeKernelDatabaseIndexesLog


    subroutine prSetKernelFromDatabase( this, gridCell, kernel,  smoothing )
        implicit none
        class( GridProjectedKDEType ), target       :: this
        type( GridCellType ), intent(inout)         :: gridCell
        class( KernelType ), target, intent(inout)  :: kernel
        doubleprecision, dimension(3), intent(in)   :: smoothing
        !-----------------------------------------------------------

        ! Restart kernel matrices
        call kernel%ResetMatrix()

        ! Compute indexes on kernel database
        call this%ComputeKernelDatabaseFlatIndexes( smoothing,   &
          gridCell%kernelDBFlatIndexes, gridCell%transposeKernel )

        ! Copy kernel from database
        call kernel%CopyFrom( this%kernelDatabaseFlat( gridCell%kernelDBFlatIndexes(1), gridCell%kernelDBFlatIndexes(2) ) )

        if ( gridCell%transposeKernel ) then
            ! Determine spans ( internally computes transpose )
            call kernel%ComputeGridSpansTranspose( gridCell%id, this%nBins,        &
              gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
              gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan, &
                                                               this%dimensionMask  )
        else
            ! Determine spans
            call kernel%ComputeGridSpans( gridCell%id, this%nBins, &
                gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
                gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan, & 
                                                                 this%dimensionMask  )
        end if 


        ! For boundaries 
        gridCell%kernelMatrix = kernel%bmatrix

        ! Done
        return

    end subroutine prSetKernelFromDatabase


    subroutine prSetKernelSigmaFromDatabase( this, gridCell, kernel, smoothing )
        implicit none
        class( GridProjectedKDEType ), target        :: this
        type( GridCellType ),  intent(inout)         :: gridCell
        class( KernelType ), target, intent(inout)   :: kernel
        doubleprecision, dimension(3), intent(in)    :: smoothing
        !-----------------------------------------------------------

        ! Restart kernel matrices
        call kernel%ResetMatrix()

        ! Compute indexes on kernel database
        ! transposeKernelSigma will always be false as this kernel is isotropic
        call this%ComputeKernelDatabaseFlatIndexes( smoothing, &
            gridCell%kernelSigmaDBFlatIndexes, gridCell%transposeKernelSigma ) 

        ! Copy kernel from database
        call kernel%CopyFrom(& 
            this%kernelDatabaseFlat( gridCell%kernelSigmaDBFlatIndexes(1), gridCell%kernelSigmaDBFlatIndexes(2) ) )

        ! Determine spans
        call kernel%ComputeGridSpans( gridCell%id, this%nBins, &
            gridCell%kernelSigmaXGSpan, gridCell%kernelSigmaYGSpan, gridCell%kernelSigmaZGSpan, & 
            gridCell%kernelSigmaXMSpan, gridCell%kernelSigmaYMSpan, gridCell%kernelSigmaZMSpan, & 
                                                                             this%dimensionMask  ) 

        ! Boundary matrix ?
        gridCell%kernelSigmaMatrix = kernel%bmatrix

        ! Done
        return


    end subroutine prSetKernelSigmaFromDatabase


    subroutine prSetKernelSD1DFromDatabase( this, gridCell, kernel, smoothing )
        implicit none
        class( GridProjectedKDEType ), target        :: this
        type( GridCellType ),  intent(inout)         :: gridCell
        class( KernelType ), target, intent(inout)   :: kernel
        doubleprecision, dimension(3), intent(in)    :: smoothing
        !-----------------------------------------------------------

        ! Restart kernel matrices
        call kernel%ResetMatrix()

        ! Compute indexes on kernel database
        gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

        ! Copy kernel from database
        call kernel%CopyFrom( this%kernelSDDatabase1( gridCell%kernelSDDBIndexes(this%idDim1) ) )

        ! Determine spans
        call kernel%ComputeGridSpans( gridCell%id, this%nBins, &
           gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
           gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan, &
                                                                  this%dimensionMask  ) 

        ! For boundaries
        gridCell%kernelSD1Matrix = kernel%bmatrix

        ! Done
        return


    end subroutine prSetKernelSD1DFromDatabase


    subroutine prSetKernelSD2DFromDatabase( this, gridCell, kernel1, kernel2, smoothing )
        implicit none
        class( GridProjectedKDEType ), target      :: this
        type( GridCellType ),  intent(inout)       :: gridCell
        class( KernelType ), target, intent(inout) :: kernel1
        class( KernelType ), target, intent(inout) :: kernel2
        doubleprecision, dimension(3), intent(in)  :: smoothing
        !-----------------------------------------------------------

        ! Restart kernel matrices
        call kernel1%ResetMatrix()
        call kernel2%ResetMatrix()

        ! Compute indexes on kernel database
        gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

        ! Copy kernel from database
        call kernel1%CopyFrom( this%kernelSDDatabase1( gridCell%kernelSDDBIndexes(this%idDim1) ) )

        ! Determine spans
        call kernel1%ComputeGridSpans( gridCell%id, this%nBins, &
           gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
           gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan, &
                                                                     this%dimensionMask  ) 

        ! Copy kernel from database
        call kernel2%CopyFrom( this%kernelSDDatabase2( gridCell%kernelSDDBIndexes(this%idDim2) ) )

        ! Determine spans
        call kernel2%ComputeGridSpans( gridCell%id, this%nBins, &
           gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
           gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan, &
                                                                     this%dimensionMask  ) 

        ! For boundaries
        gridCell%kernelSD1Matrix = kernel1%bmatrix
        gridCell%kernelSD2Matrix = kernel2%bmatrix

        ! Done
        return


    end subroutine prSetKernelSD2DFromDatabase


    subroutine prSetKernelSD3DFromDatabase( this, gridCell, kernel1, kernel2, kernel3, smoothing )
        implicit none
        class( GridProjectedKDEType ), target      :: this
        type( GridCellType ),  intent(inout)       :: gridCell
        class( KernelType ), target, intent(inout) :: kernel1
        class( KernelType ), target, intent(inout) :: kernel2
        class( KernelType ), target, intent(inout) :: kernel3
        doubleprecision, dimension(3), intent(in)  :: smoothing
        !-----------------------------------------------------------

        ! Compute indexes on kernel database
        gridCell%kernelSDDBIndexes = this%ComputeKernelDatabaseIndexes( smoothing )

        ! Copy kernel from database
        call kernel1%CopyFrom( this%kernelSDXDatabase( gridCell%kernelSDDBIndexes(1) ) )

        ! Determine spans
        call kernel1%ComputeGridSpans( gridCell%id, this%nBins, &
           gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
           gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan, &
                                                                     this%dimensionMask  )

        ! Copy kernel from database
        call kernel2%CopyFrom( this%kernelSDYDatabase( gridCell%kernelSDDBIndexes(2) ) )

        ! Determine spans
        call kernel2%ComputeGridSpans( gridCell%id, this%nBins, &
           gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
           gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan, &
                                                                     this%dimensionMask  )

        ! Copy kernel from database
        call kernel3%CopyFrom( this%kernelSDZDatabase( gridCell%kernelSDDBIndexes(3) ) )

        ! Determine spans
        call kernel3%ComputeGridSpans( gridCell%id, this%nBins, &
           gridCell%kernelSD3XGSpan, gridCell%kernelSD3YGSpan, gridCell%kernelSD3ZGSpan, & 
           gridCell%kernelSD3XMSpan, gridCell%kernelSD3YMSpan, gridCell%kernelSD3ZMSpan, &
                                                                     this%dimensionMask  )

        ! For boundaries
        gridCell%kernelSD1Matrix = kernel1%bmatrix
        gridCell%kernelSD2Matrix = kernel2%bmatrix
        gridCell%kernelSD3Matrix = kernel3%bmatrix

        ! Done
        return


    end subroutine prSetKernelSD3DFromDatabase


    subroutine prSetKernelBrute( this, gridCell, kernel,  smoothing )
        implicit none
        class( GridProjectedKDEType ), target      :: this
        type( GridCellType ), intent(inout)        :: gridCell
        class( KernelType ), target, intent(inout) :: kernel
        doubleprecision, dimension(3), intent(in)  :: smoothing
        !-----------------------------------------------------------

        ! Restart kernel matrices
        call kernel%ResetMatrix()

        ! Setup kernel matrix
        call kernel%SetupMatrix( smoothing )

        ! Determine spans
        call kernel%ComputeGridSpans( gridCell%id, this%nBins   , &
            gridCell%kernelXGSpan, gridCell%kernelYGSpan, gridCell%kernelZGSpan, & 
            gridCell%kernelXMSpan, gridCell%kernelYMSpan, gridCell%kernelZMSpan, &
                                                             this%dimensionMask  )

        ! For boundaries 
        gridCell%kernelMatrix = kernel%bmatrix

        ! Done 
        return


    end subroutine prSetKernelBrute


    subroutine prSetKernelSigmaBrute( this, gridCell, kernel, smoothing )
        implicit none
        class( GridProjectedKDEType ), target      :: this
        type( GridCellType ), intent(inout)        :: gridCell
        class( KernelType ), target, intent(inout) :: kernel
        doubleprecision, dimension(3), intent(in)  :: smoothing
        !-----------------------------------------------------------

        ! Restart kernel matrices
        call kernel%ResetMatrix()

        ! Setup kernel matrix
        call kernel%SetupMatrix( smoothing )

        ! Determine spans
        call kernel%ComputeGridSpans( gridCell%id, this%nBins, &
            gridCell%kernelSigmaXGSpan, gridCell%kernelSigmaYGSpan, gridCell%kernelSigmaZGSpan, & 
            gridCell%kernelSigmaXMSpan, gridCell%kernelSigmaYMSpan, gridCell%kernelSigmaZMSpan, &
                                                                            this%dimensionMask  )

        gridCell%kernelSigmaMatrix = kernel%bmatrix


        ! Done
        return


    end subroutine prSetKernelSigmaBrute


    subroutine prSetKernelSD1DBrute( this, gridCell, kernel, smoothing )
        implicit none
        class( GridProjectedKDEType ), target      :: this
        type( GridCellType ),  intent(inout)       :: gridCell
        class( KernelType ), target, intent(inout) :: kernel
        doubleprecision, dimension(3), intent(in)  :: smoothing
        !-----------------------------------------------------------

        ! Restart kernel matrices
        call kernel%ResetMatrix()

        ! Compute matrix
        call kernel%SetupMatrix( smoothing )

        ! Determine spans
        call kernel%ComputeGridSpans( gridCell%id, this%nBins, &
           gridCell%kernelSDXGSpan, gridCell%kernelSDYGSpan, gridCell%kernelSDZGSpan, & 
           gridCell%kernelSDXMSpan, gridCell%kernelSDYMSpan, gridCell%kernelSDZMSpan, & 
                                                                   this%dimensionMask )

        ! For boundaries
        gridCell%kernelSD1Matrix = kernel%bmatrix

        ! Done
        return


    end subroutine prSetKernelSD1DBrute


    subroutine prSetKernelSD2DBrute( this, gridCell, kernel1, kernel2, smoothing )
        implicit none
        class( GridProjectedKDEType ), target      :: this
        type( GridCellType ),  intent(inout)       :: gridCell
        class( KernelType ), target, intent(inout) :: kernel1
        class( KernelType ), target, intent(inout) :: kernel2
        doubleprecision, dimension(3), intent(in)  :: smoothing
        doubleprecision, dimension(3)              :: locSmoothing
        !-----------------------------------------------------------

        ! Restart kernel matrices
        call kernel1%ResetMatrix()
        call kernel2%ResetMatrix()

        ! Fill smoothing, curvature kernel is isotropic
        locSmoothing = 0d0
        locSmoothing(this%idDim1) = smoothing(this%idDim1)
        locSmoothing(this%idDim2) = smoothing(this%idDim1)

        ! Compute matrix
        call kernel1%SetupMatrix( locSmoothing )

        ! Determine spans
        call kernel1%ComputeGridSpans( gridCell%id, this%nBins, &
           gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
           gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan, & 
                                                                      this%dimensionMask )

        ! Fill smoothing, curvature kernel is isotropic
        locSmoothing = 0d0
        locSmoothing(this%idDim1) = smoothing(this%idDim2)
        locSmoothing(this%idDim2) = smoothing(this%idDim2)

        ! Compute matrix
        call kernel2%SetupMatrix( locSmoothing )

        ! Determine spans
        call kernel2%ComputeGridSpans( gridCell%id, this%nBins, &
           gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
           gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan, & 
                                                                      this%dimensionMask )

        ! For boundaries
        gridCell%kernelSD1Matrix = kernel1%bmatrix
        gridCell%kernelSD2Matrix = kernel2%bmatrix

        ! Done
        return


    end subroutine prSetKernelSD2DBrute


    subroutine prSetKernelSD3DBrute( this, gridCell, kernel1, kernel2, kernel3, smoothing )
        implicit none
        class( GridProjectedKDEType ), target      :: this
        type( GridCellType ),  intent(inout)       :: gridCell
        class( KernelType ), target, intent(inout) :: kernel1
        class( KernelType ), target, intent(inout) :: kernel2
        class( KernelType ), target, intent(inout) :: kernel3
        doubleprecision, dimension(3), intent(in)  :: smoothing
        !-----------------------------------------------------------

        ! Restart kernel matrices
        call kernel1%ResetMatrix()
        call kernel2%ResetMatrix()
        call kernel3%ResetMatrix()

        ! Compute matrix
        call kernel1%SetupMatrix( (/smoothing(1),smoothing(1),smoothing(1)/) )

        ! Determine spans
        call kernel1%ComputeGridSpans( gridCell%id, this%nBins, &
           gridCell%kernelSD1XGSpan, gridCell%kernelSD1YGSpan, gridCell%kernelSD1ZGSpan, & 
           gridCell%kernelSD1XMSpan, gridCell%kernelSD1YMSpan, gridCell%kernelSD1ZMSpan, & 
                                                                      this%dimensionMask )

        ! Compute matrix
        call kernel2%SetupMatrix( (/smoothing(2),smoothing(2),smoothing(2)/) )

        ! Determine spans
        call kernel2%ComputeGridSpans( gridCell%id, this%nBins, &
           gridCell%kernelSD2XGSpan, gridCell%kernelSD2YGSpan, gridCell%kernelSD2ZGSpan, & 
           gridCell%kernelSD2XMSpan, gridCell%kernelSD2YMSpan, gridCell%kernelSD2ZMSpan, & 
                                                                      this%dimensionMask )

        ! Compute matrix
        call kernel3%SetupMatrix( (/smoothing(3),smoothing(3),smoothing(3)/) )

        ! Determine spans
        call kernel3%ComputeGridSpans( gridCell%id, this%nBins, &
           gridCell%kernelSD3XGSpan, gridCell%kernelSD3YGSpan, gridCell%kernelSD3ZGSpan, & 
           gridCell%kernelSD3XMSpan, gridCell%kernelSD3YMSpan, gridCell%kernelSD3ZMSpan, & 
                                                                      this%dimensionMask )

        ! For boundaries
        gridCell%kernelSD1Matrix = kernel1%bmatrix
        gridCell%kernelSD2Matrix = kernel2%bmatrix
        gridCell%kernelSD3Matrix = kernel3%bmatrix


        ! Done
        return


    end subroutine prSetKernelSD3DBrute


    ! Utils kernel database
    function prGenerateLogSpaceData( this, initPoint, endPoint, nPoints ) result( output )
      !------------------------------------------------------------------------------
      !
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      implicit none 
      class(GridProjectedKDEType) :: this
      doubleprecision, intent(in) :: initPoint, endPoint
      integer, intent(in)         :: nPoints
      doubleprecision :: deltaExponent, deltaAccumulated
      doubleprecision :: logInit, logEnd, logBase
      doubleprecision, dimension(:), allocatable :: output
      integer :: n
      !-----------------------------------------

      allocate( output( nPoints ) )

      logBase = 10d0
      deltaAccumulated = 0d0
      ! Some sanity to verify init smaller than end

      logInit        = log10( initPoint )
      logEnd         = log10( endPoint  )
      deltaExponent  = ( logEnd - logInit )/( nPoints - 1 )  

      do n = 1, nPoints
        output(n) = logBase**( logInit + deltaAccumulated )
        deltaAccumulated = deltaAccumulated + deltaExponent
      end do
      
      return

    end function prGenerateLogSpaceData



    function prComputeXYTranspose( this, sourceMatrix ) result( transposedMatrix )
      !------------------------------------------------------------------------------
      !
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      implicit none 
      class(GridProjectedKDEType) :: this
      doubleprecision, dimension(:,:,:), intent(in)  :: sourceMatrix
      doubleprecision, dimension(:,:,:), allocatable :: transposedMatrix
      ! local
      integer, dimension(3) :: sourceShape
      integer :: n
      !------------------------------------------------------------------------------
      
      sourceShape = shape( sourceMatrix )

      allocate( transposedMatrix( sourceShape(2), sourceShape(1), sourceShape(3) ) )

      do n = 1, sourceShape(3)
          transposedMatrix(:,:,n) = transpose( sourceMatrix(:,:,n) )
      end do
  
      return 

    end function prComputeXYTranspose


    ! Utils output files
    subroutine prExportDensity( this, outputFileName )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(GridProjectedKDEType) :: this
        character(len=*), intent(in) :: outputFileName
        integer :: ix, iy, iz
        integer :: outputUnit = 555
        !------------------------------------------------------------------------------

        ! Write the output file name
        ! Add some default
        open( outputUnit, file=outputFileName, status='replace' )

        ! Following column-major nesting
        do iz = 1, this%nBins(3)
            do iy = 1, this%nBins(2)
                do ix = 1, this%nBins(1)
                    if ( this%densityEstimateGrid( ix, iy, iz ) .le. 0d0 ) cycle
                    ! THIS FORMAT MAY BE DYNAMIC ACCORDING TO THE TOTAL NUMBER OF PARTICLES
                    !write(outputUnit,"(I8,I8,I8,es18.9e3,I8)") ix, iy, iz, & 
                    write(outputUnit,"(I8,I8,I8,2es18.9e3)") ix, iy, iz, & 
                      this%densityEstimateGrid( ix, iy, iz ), this%histogram%counts( ix, iy, iz ) 
                end do
            end do
        end do

        ! Finished
        close(outputUnit)


    end subroutine prExportDensity


    subroutine prExportDensityUnit( this, outputUnit, outputDataId, particleGroupId )
      !------------------------------------------------------------------------------
      ! 
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      implicit none 
      class(GridProjectedKDEType) :: this
      integer, intent(in) :: outputUnit
      integer, intent(in) :: outputDataId
      integer, intent(in) :: particleGroupId
      integer :: ix, iy, iz
      integer :: countNonZero, counter
      !------------------------------------------------------------------------------

      counter      = 0
      countNonZero = count(this%densityEstimateGrid /=0d0)
      if( allocated( this%outputBinIds ) ) deallocate( this%outputBinIds )
      allocate( this%outputBinIds(countNonZero,3) )
       
      ! Following column-major nesting
      do iz = 1, this%nBins(3)
        do iy = 1, this%nBins(2)
          do ix = 1, this%nBins(1)
            if ( this%densityEstimateGrid( ix, iy, iz ) .le. 0d0 ) cycle
            ! THIS FORMAT MAY BE DYNAMIC ACCORDING TO THE TOTAL NUMBER OF PARTICLES/COLUMNS
            ! write(outputUnit,"(I8,I8,I8,I8,I6,es18.9e3,I8)") outputDataId, particleGroupId, &
            write(outputUnit,"(I8,I8,I8,I8,I6,2es18.9e3)") outputDataId, particleGroupId, &
                ix, iy, iz, this%densityEstimateGrid( ix, iy, iz ), &
                               this%histogram%counts( ix, iy, iz ) 
            counter = counter + 1 
            this%outputBinIds(counter,:) = (/ix,iy,iz/)
          end do
        end do
      end do


    end subroutine prExportDensityUnit


    subroutine prExportOptimizationVariables( this, outputFileName, &
        densityEstimateArray, kernelSmoothing, kernelSigmaSupportScale, &
        curvatureBandwidth, nEstimate, netRoughness )
        !------------------------------------------------------------------------------
        implicit none 
        class(GridProjectedKDEType) :: this
        character(len=500), intent(in) :: outputFileName
        doubleprecision, dimension(:)  ,intent(in) :: densityEstimateArray
        doubleprecision, dimension(:,:),intent(in) :: kernelSmoothing 
        doubleprecision, dimension(:)  ,intent(in) :: kernelSigmaSupportScale
        doubleprecision, dimension(:,:),intent(in) :: curvatureBandwidth
        doubleprecision, dimension(:)  ,intent(in) :: nEstimate
        doubleprecision, dimension(:)  ,intent(in) :: netRoughness
        integer :: ix, iy, iz, n
        integer :: outputUnit = 555
        !------------------------------------------------------------------------------

        ! Write the output file name
        ! Add some default
        open( outputUnit, file=outputFileName, status='replace' )


        do n = 1, this%nComputeBins
            ix = this%computeBinIds( 1, n )
            iy = this%computeBinIds( 2, n )
            iz = this%computeBinIds( 3, n )
            ! THIS FORMAT MAY BE DYNAMIC ACCORDING TO THE TOTAL NUMBER OF PARTICLES
            write(outputUnit,&
                "(I6,I6,I6,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8,F16.8)") &
                ix, iy, iz,& 
                densityEstimateArray( n ),& 
                kernelSmoothing(1,n), kernelSmoothing(2,n), kernelSmoothing(3,n),& 
                kernelSigmaSupportScale(n), &
                curvatureBandwidth(1,n), curvatureBandwidth(2,n), curvatureBandwidth(3,n), &
                nEstimate(n), netRoughness(n)
        end do

        ! Finished
        close(outputUnit)


    end subroutine prExportOptimizationVariables


    subroutine prExportOptimizationVariablesExtended( this, outputFileName, &
      densityEstimateArray, kernelSmoothing, kernelSmoothingScale,        & 
                           kernelSmoothingShape, kernelSigmaSupportScale, &
      curvatureBandwidth, nEstimate, roughnessXXArray, roughnessYYArray, &
                                           roughnessZZArray, netRoughness )
      !------------------------------------------------------------------------------
      implicit none 
      class(GridProjectedKDEType) :: this
      character(len=500), intent(in) :: outputFileName
      doubleprecision, dimension(:)  ,intent(in) :: densityEstimateArray
      doubleprecision, dimension(:,:),intent(in) :: kernelSmoothing 
      doubleprecision, dimension(:)  ,intent(in) :: kernelSmoothingScale
      doubleprecision, dimension(:,:),intent(in) :: kernelSmoothingShape
      doubleprecision, dimension(:)  ,intent(in) :: kernelSigmaSupportScale
      doubleprecision, dimension(:,:),intent(in) :: curvatureBandwidth
      doubleprecision, dimension(:)  ,intent(in) :: nEstimate
      doubleprecision, dimension(:)  ,intent(in) :: netRoughness
      doubleprecision, dimension(:)  ,intent(in) :: roughnessXXArray   
      doubleprecision, dimension(:)  ,intent(in) :: roughnessYYArray
      doubleprecision, dimension(:)  ,intent(in) :: roughnessZZArray
      integer :: ix, iy, iz, n
      integer :: outputUnit = 555
      !------------------------------------------------------------------------------

      ! Write the output file name
      ! Add some default
      open( outputUnit, file=outputFileName, status='replace' )

      do n = 1, this%nComputeBins
        ix = this%computeBinIds( 1, n )
        iy = this%computeBinIds( 2, n )
        iz = this%computeBinIds( 3, n )
        !"(I6,I6,I6,es18.9e3,3es18.9e3,3es18.9e3,es18.9e3,es18.9e3,3es18.9e3,es18.9e3, & 
        !  3es18.9e3,es18.9e3)")  &
        write(outputUnit,"(3I6,17es18.9e3)") ix, iy, iz,& 
          densityEstimateArray( n ),& 
          kernelSmoothing(1,n), kernelSmoothing(2,n), kernelSmoothing(3,n),& 
          kernelSmoothingShape(1,n), kernelSmoothingShape(2,n), kernelSmoothingShape(3,n),& 
          kernelSmoothingScale(n), kernelSigmaSupportScale(n), &
          curvatureBandwidth(1,n), curvatureBandwidth(2,n), curvatureBandwidth(3,n), &
          nEstimate(n), roughnessXXArray(n), roughnessYYArray(n), &
          roughnessZZArray(n), netRoughness(n)
      end do

      ! Finished
      close(outputUnit)

    end subroutine prExportOptimizationVariablesExtended



    subroutine prExportOptimizationVariablesExtendedError( this, outputFileName, &
        densityEstimateArray, kernelSmoothing, kernelSmoothingScale,        & 
                             kernelSmoothingShape, kernelSigmaSupportScale, &
        curvatureBandwidth, nEstimate, roughnessXXArray, roughnessYYArray, &
                                             roughnessZZArray, netRoughness, &
                                relativeDensityChange, relativeRoughnessChange )
        !------------------------------------------------------------------------------
        implicit none 
        class(GridProjectedKDEType) :: this
        character(len=500), intent(in) :: outputFileName
        doubleprecision, dimension(:)  ,intent(in) :: densityEstimateArray
        doubleprecision, dimension(:,:),intent(in) :: kernelSmoothing 
        doubleprecision, dimension(:)  ,intent(in) :: kernelSmoothingScale
        doubleprecision, dimension(:,:),intent(in) :: kernelSmoothingShape
        doubleprecision, dimension(:)  ,intent(in) :: kernelSigmaSupportScale
        doubleprecision, dimension(:,:),intent(in) :: curvatureBandwidth
        doubleprecision, dimension(:)  ,intent(in) :: nEstimate
        doubleprecision, dimension(:)  ,intent(in) :: netRoughness
        doubleprecision, dimension(:)  ,intent(in) :: roughnessXXArray   
        doubleprecision, dimension(:)  ,intent(in) :: roughnessYYArray
        doubleprecision, dimension(:)  ,intent(in) :: roughnessZZArray
        doubleprecision, dimension(:)  ,intent(in) :: relativeDensityChange
        doubleprecision, dimension(:)  ,intent(in) :: relativeRoughnessChange
        integer :: ix, iy, iz, n
        integer :: outputUnit = 555
        !------------------------------------------------------------------------------

        ! Write the output file name
        ! Add some default
        open( outputUnit, file=outputFileName, status='replace' )


        do n = 1, this%nComputeBins
            ix = this%computeBinIds( 1, n )
            iy = this%computeBinIds( 2, n )
            iz = this%computeBinIds( 3, n )
            !    "(I6,I6,I6,es18.9e3,3es18.9e3,3es18.9e3,& 
            !      es18.9e3,es18.9e3,3es18.9e3,es18.9e3, & 
            !      3es18.9e3,es18.9e3,es18.9e3,es18.9e3)")  &
            write(outputUnit,"(3I6,19es18.9e3)") ix, iy, iz, & 
                densityEstimateArray( n ),& 
                kernelSmoothing(1,n), kernelSmoothing(2,n), kernelSmoothing(3,n),& 
                kernelSmoothingShape(1,n), kernelSmoothingShape(2,n), kernelSmoothingShape(3,n),& 
                kernelSmoothingScale(n), kernelSigmaSupportScale(n), &
                curvatureBandwidth(1,n), curvatureBandwidth(2,n), curvatureBandwidth(3,n), &
                nEstimate(n), roughnessXXArray(n), roughnessYYArray(n), &
                roughnessZZArray(n), netRoughness(n), relativeDensityChange(n),&
                relativeRoughnessChange(n)
        end do

        ! Finished
        close(outputUnit)


    end subroutine prExportOptimizationVariablesExtendedError


    subroutine prWriteErrorMetricsRecord(this, outputUnit, loopId,   &  
              relativeDensity, relativeRoughness, relativeSmoothing, &
           nFractionDensity, nFractionRoughness, nFractionSmoothing, &
                                              error1, error2, error3 )
        !------------------------------------------------------------------------------
        !------------------------------------------------------------------------------
        implicit none 
        class(GridProjectedKDEType) :: this
        integer,intent(in) :: outputUnit, loopId
        doubleprecision, intent(in) :: relativeDensity
        doubleprecision, intent(in) :: relativeRoughness
        doubleprecision, intent(in) :: relativeSmoothing
        doubleprecision, intent(in) :: nFractionDensity
        doubleprecision, intent(in) :: nFractionRoughness
        doubleprecision, intent(in) :: nFractionSmoothing
        doubleprecision, intent(in) :: error1, error2, error3
        !------------------------------------------------------------------------------
    
        write(outputUnit, '(I8,9es18.9e3)') loopId, &
          relativeDensity, nFractionDensity,     & 
          relativeRoughness, nFractionRoughness, & 
          relativeSmoothing, nFractionSmoothing, & 
          error1, error2, error3 


    end subroutine prWriteErrorMetricsRecord


end module GridProjectedKDEModule
