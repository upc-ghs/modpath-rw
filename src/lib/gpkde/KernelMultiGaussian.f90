module KernelMultiGaussianModule
  !--------------------------------------------------------------------------------
  ! Module that provides grid-projected MultiGaussian kernels 
  !--------------------------------------------------------------------------------
  use PrecisionModule, only : fp 
  use ConstantsModule, only : fZERO, fONE, fTWO, fTHREE, fFIVE, pi, sqrtPi, sqrtTwo
  implicit none
  !--------------------------------------------------------------------------------

  ! Parameters
  integer, parameter                 :: defaultKernelRange   = 3
  integer, parameter                 :: defaultKernelSDRange = 4
  ! Common vars, will be updated according to init function
  integer                            :: nDim          = 3
  real(fp)                           :: fNDim         = 3.0_fp
  integer, dimension(3)              :: dimensionMask = (/1,1,1/)
  integer                            :: idDim1, idDim2
  integer, dimension(:), allocatable :: dimensions

  ! Set default access status to private
  private

  type, public, abstract :: KernelType
      
    ! Properties
    real(fp), dimension(3) :: binSize     = fZERO
    real(fp), dimension(3) :: bandwidth   = fZERO
    integer                :: matrixRange 
    integer , dimension(3) :: matrixPositiveShape = 0 
    logical                :: shouldIntegrateOne 
    real(fp), dimension(:,:,:), allocatable :: matrix

  contains 
      
    ! Procedures
    procedure, non_overridable :: Initialize                   => prInitialize 
    procedure, non_overridable :: Reset                        => prReset      
    procedure, non_overridable :: ResetMatrix                  => prResetMatrix
    procedure, non_overridable :: ResetMatrixOld               => prResetMatrixOld
    procedure, non_overridable :: ComputeSpansBounded          => prComputeSpansBounded
    procedure, non_overridable :: ComputeSpansBoundedTranspose => prComputeSpansBoundedTranspose
    procedure, non_overridable :: GenerateZeroPositiveGrid     => prGenerateZeroPositiveGrid
    procedure, non_overridable :: UnfoldZeroPositiveMatrix     => prUnfoldZeroPositiveMatrix
    procedure, non_overridable :: SetupMatrix                  => prSetupMatrix
    procedure, non_overridable :: CopyFrom                     => prCopyFrom
    procedure( ComputeKernelMatrix ), deferred  :: ComputeMatrix 

  end type

  ! Management of module dimensions 
  public :: InitializeKernelDimensions
  public :: ResetKernelDimensions
      
  ! MultiGaussian W
  type, public, extends( KernelType ) :: KernelMultiGaussianType
    ! Properties
  contains
    ! Procedures
    procedure :: ComputeMatrix => prComputeKernelWMatrix
  end type
  
  ! Second Derivative X
  type, public, extends( KernelType ) :: KernelSecondDerivativeXType
    ! Properties
  contains
    ! Procedures
    procedure :: ComputeMatrix => prComputeKernelVXMatrix
  end type

  ! Second Derivative Y
  type, public, extends( KernelType ) :: KernelSecondDerivativeYType
    ! Properties
  contains
    ! Procedures
    procedure :: ComputeMatrix => prComputeKernelVYMatrix
  end type

  ! Second Derivative Z
  type, public, extends( KernelType ) :: KernelSecondDerivativeZType
    ! Properties
  contains
    ! Procedures
    procedure :: ComputeMatrix => prComputeKernelVZMatrix
  end type

  ! Interfaces
  abstract interface
    subroutine ComputeKernelMatrix( this, zPXGrid, zPYgrid, zPZGrid )
      !------------------------------------------------------------------------------
      !
      !------------------------------------------------------------------------------
      ! Specifications 
      !------------------------------------------------------------------------------
      import :: KernelType
      import :: fp
      class( KernelType )                     :: this 
      integer , dimension(:,:,:), intent(in)  :: zPXGrid, zPYgrid, zPZGrid
      real(fp), dimension(:)    , allocatable :: hDelta
      real(fp), dimension(:,:,:), allocatable :: zeroPositiveMatrix
      integer :: nx, ny, nz
      !------------------------------------------------------------------------------
    end subroutine ComputeKernelMatrix
  end interface

contains

  ! Common object methods
  subroutine prInitialize( this, binSize, matrixRange, dimensionMask )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( KernelType )                           :: this 
    real(fp), dimension(:)                        :: binSize
    integer , intent(in), optional                :: matrixRange
    integer , dimension(:), intent(in), optional  :: dimensionMask
    !------------------------------------------------------------------------------

    ! Assign binSize 
    this%binSize = binSize 

    ! Range for computing kernel extent
    if ( present( matrixRange ) ) then 
      this%matrixRange = matrixRange
    else
      this%matrixRange = defaultKernelRange 
    end if


  end subroutine prInitialize


  subroutine prReset( this )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( KernelType ) :: this 
    !------------------------------------------------------------------------------

    this%bandwidth = fZERO
    this%binSize   = fZERO
    this%matrixPositiveShape = 0
    this%matrixRange = defaultKernelRange ! Set to zero ?

    if ( allocated( this%matrix ) ) deallocate( this%matrix )

  end subroutine prReset


  subroutine prResetMatrix( this )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( KernelType ) :: this 
    !------------------------------------------------------------------------------

    if ( allocated( this%matrix ) ) deallocate( this%matrix )
    this%bandwidth = fZERO
    this%matrixPositiveShape = 0

  end subroutine prResetMatrix


  subroutine prResetMatrixOld( this )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( KernelType ) :: this 
    !------------------------------------------------------------------------------

    if ( allocated( this%matrix ) ) deallocate( this%matrix )

    this%bandwidth = fZERO
    this%matrixPositiveShape = 0


  end subroutine prResetMatrixOld


  subroutine prCopyFrom( this, source )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( KernelType ), intent(inout) :: this 
    class( KernelType ), intent(in)   :: source 
    !------------------------------------------------------------------------------

    this%matrixPositiveShape = source%matrixPositiveShape
    this%matrix              = source%matrix

  end subroutine prCopyFrom


  subroutine prComputeSpansBounded( this, gridIndexes, gridShape, &
                                 xGridSpan, yGridSpan, zGridSpan, &
                           xKernelSpan, yKernelSpan, zKernelSpan  )
    !------------------------------------------------------------------------------
    ! Calculate both grid and kernel spans and apply 
    ! reflection boundary correction if necessary
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( KernelType )                  :: this
    integer, dimension(3), intent(in)    :: gridShape
    integer, dimension(3), intent(in)    :: gridIndexes
    integer, dimension(2), intent(inout) :: xGridSpan, yGridSpan, zGridSpan
    integer, dimension(2), intent(inout) :: xKernelSpan, yKernelSpan, zKernelSpan
    logical, dimension(3)                :: isBoundary
    integer, dimension(3)                :: boundDir
    integer, dimension(3)                :: boundLoc
    integer                              :: nd, did
    !------------------------------------------------------------------------------

    ! Spans in grid (ORIGINAL) 
    xGridSpan(1) = max( gridIndexes(1) - this%matrixPositiveShape(1), 1)
    xGridSpan(2) = min( gridIndexes(1) + this%matrixPositiveShape(1), gridShape(1) )
    yGridSpan(1) = max( gridIndexes(2) - this%matrixPositiveShape(2), 1)
    yGridSpan(2) = min( gridIndexes(2) + this%matrixPositiveShape(2), gridShape(2) )
    zGridSpan(1) = max( gridIndexes(3) - this%matrixPositiveShape(3), 1)
    zGridSpan(2) = min( gridIndexes(3) + this%matrixPositiveShape(3), gridShape(3) )

    ! Spans in kernel matrix
    xKernelSpan = xGridSpan + this%matrixPositiveShape(1) - gridIndexes(1) + 1
    yKernelSpan = yGridSpan + this%matrixPositiveShape(2) - gridIndexes(2) + 1
    zKernelSpan = zGridSpan + this%matrixPositiveShape(3) - gridIndexes(3) + 1

    isBoundary = .false.
    boundLoc   = 0
    boundDir   = 0
    do nd=1,nDim
      did = dimensions(nd)
      if ( ( gridIndexes(did) - this%matrixPositiveShape(did) ) .lt. 1 ) then 
        isBoundary(did) = .true.
        boundLoc(did)   = 1 - gridIndexes(did) + this%matrixPositiveShape(did)
        boundDir(did)   = 1
      else if ( ( gridIndexes(did) + this%matrixPositiveShape(did) ) .gt. gridShape(did) ) then 
        isBoundary(did) = .true.
        boundLoc(did)   = 2*this%matrixPositiveShape(did) + 1 &
          - ( gridIndexes(did) + this%matrixPositiveShape(did) - gridShape(did) ) + 1 
        boundDir(did)   = 2
      end if 
    end do 

    if ( any( isBoundary ) ) then 
      call prKernelReflection( this, gridIndexes, gridShape, &
                            xGridSpan, yGridSpan, zGridSpan, &
                      xKernelSpan, yKernelSpan, zKernelSpan, &
                             isBoundary, boundLoc, boundDir  )
    end if 

  end subroutine prComputeSpansBounded


  subroutine prKernelReflection( this, gridIndexes, gridShape,  &
                              xGridSpan, yGridSpan, zGridSpan,  &
                         xKernelSpan, yKernelSpan, zKernelSpan, & 
                                 isBoundary, boundLoc, boundDir )
    !------------------------------------------------------------------------------
    ! Apply kernel reflection directly on this%matrix 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( KernelType ) :: this
    integer, dimension(3), intent(in)    :: gridShape
    integer, dimension(3), intent(in)    :: gridIndexes
    integer, dimension(2), intent(inout) :: xGridSpan, yGridSpan, zGridSpan
    integer, dimension(2), intent(inout) :: xKernelSpan, yKernelSpan, zKernelSpan
    logical, dimension(3), intent(in)    :: isBoundary 
    integer, dimension(3), intent(in)    :: boundLoc
    integer, dimension(3), intent(in)    :: boundDir
    integer, dimension(:), allocatable   :: kernelShape
    integer                              :: lenb, nd, did
    !------------------------------------------------------------------------------

    ! If boundary
    kernelShape  = shape(this%matrix)

    ! Kernel reflection ! 
    do nd=1,nDim
      did = dimensions(nd)
      if ( .not. isBoundary(did) ) cycle
      select case(did)
      case(1)
        ! Reflection X
        select case( boundDir(did) ) 
        case(1)
          ! WEST
          lenb = boundLoc(did)
          this%matrix( boundLoc(did) + 1: boundLoc(did) + lenb, :, :) = &
          this%matrix( boundLoc(did) + 1: boundLoc(did) + lenb, :, :) + &
          this%matrix( boundLoc(did):1:-1, :, :)
          this%matrix( :boundLoc(did), :, :) = 0
        case(2)
          ! EAST 
          lenb = kernelShape(did) - boundLoc(did) + 1 
          this%matrix( boundLoc(did) - lenb: boundLoc(did) - 1, :, :) = &
          this%matrix( boundLoc(did) - lenb: boundLoc(did) - 1, :, :) + &
          this%matrix( kernelShape(did): boundLoc(did) :-1, :, :)
          this%matrix( boundLoc(did):, :, :) = 0
        end select    
      case(2)
        ! Reflection Y
        select case( boundDir(did) ) 
        case(1)
          ! SOUTH
          lenb = boundLoc(did)
          this%matrix( :, boundLoc(did) + 1: boundLoc(did) + lenb, :) = &
          this%matrix( :, boundLoc(did) + 1: boundLoc(did) + lenb, :) + &
          this%matrix( :, boundLoc(did):1:-1, :)
          this%matrix( :, :boundLoc(did), :) = 0
        case(2)
          ! NORTH
          lenb = kernelShape(did) - boundLoc(did) + 1 
          this%matrix( :, boundLoc(did) - lenb: boundLoc(did) - 1, :) = &
          this%matrix( :, boundLoc(did) - lenb: boundLoc(did) - 1, :) + &
          this%matrix( :, kernelShape(did): boundLoc(did) :-1, :)
          this%matrix( :, boundLoc(did):, :) = 0
        end select    
      case(3)
        ! Reflection Z
        select case( boundDir(did) ) 
        case(1)
          ! BOTTOM
          lenb = boundLoc(did)
          this%matrix( :, :, boundLoc(did) + 1: boundLoc(did) + lenb) = &
          this%matrix( :, :, boundLoc(did) + 1: boundLoc(did) + lenb) + &
          this%matrix( :, :, boundLoc(did):1:-1)
          this%matrix( :, :, :boundLoc(did)) = 0
        case(2)
          ! TOP
          lenb = kernelShape(did) - boundLoc(did) + 1 
          this%matrix( :, :, boundLoc(did) - lenb: boundLoc(did) - 1) = &
          this%matrix( :, :, boundLoc(did) - lenb: boundLoc(did) - 1) + &
          this%matrix( :, :, kernelShape(did): boundLoc(did) :-1)
          this%matrix( :, :, boundLoc(did):) = 0
        end select
      end select
    end do 


  end subroutine prKernelReflection 


  subroutine prComputeSpansBoundedTranspose( this, gridIndexes, gridShape, &
                                          xGridSpan, yGridSpan, zGridSpan, &
                                    xKernelSpan, yKernelSpan, zKernelSpan  ) 
    !------------------------------------------------------------------------------
    !  
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( KernelType ) :: this
    integer, dimension(3), intent(in) :: gridShape
    integer, dimension(3), intent(in) :: gridIndexes
    integer, dimension(2), intent(inout) :: xGridSpan, yGridSpan, zGridSpan
    integer, dimension(2), intent(inout) :: xKernelSpan, yKernelSpan, zKernelSpan
    logical, dimension(3)                :: isBoundary
    integer, dimension(3)                :: boundDir
    integer, dimension(3)                :: boundLoc
    !integer                              :: nd, did
    !------------------------------------------------------------------------------

    ! Spans in grid
    !( notice the switch in matrixPositiveShape indexes,
    !  compare against ComputeGridSpans )
    xGridSpan(1) = max( gridIndexes(1) - this%matrixPositiveShape(2), 1)
    xGridSpan(2) = min( gridIndexes(1) + this%matrixPositiveShape(2), gridShape(1) )
    yGridSpan(1) = max( gridIndexes(2) - this%matrixPositiveShape(1), 1)
    yGridSpan(2) = min( gridIndexes(2) + this%matrixPositiveShape(1), gridShape(2) )
    zGridSpan(1) = max( gridIndexes(3) - this%matrixPositiveShape(3), 1)
    zGridSpan(2) = min( gridIndexes(3) + this%matrixPositiveShape(3), gridShape(3) )

    ! Spans in transposed kernel matrix
    xKernelSpan = xGridSpan + this%matrixPositiveShape(2) - gridIndexes(1) + 1
    yKernelSpan = yGridSpan + this%matrixPositiveShape(1) - gridIndexes(2) + 1
    zKernelSpan = zGridSpan + this%matrixPositiveShape(3) - gridIndexes(3) + 1

    isBoundary = .false.
    boundLoc   = 0
    boundDir   = 0

    ! X
    if ( dimensionMask(1) .eq. 1 ) then 
      if ( ( gridIndexes(1) - this%matrixPositiveShape(2) ) .lt. 1 ) then 
        isBoundary(1)  = .true.
        boundLoc(1)    = 1 - gridIndexes(1) + this%matrixPositiveShape(2)
        boundDir(1)    = 1
      else if ( ( gridIndexes(1) + this%matrixPositiveShape(2) ) .gt. gridShape(1) ) then 
        isBoundary(1)  = .true.
        boundLoc(1)    = 2*this%matrixPositiveShape(2) + 1 &
            - ( gridIndexes(1) + this%matrixPositiveShape(2) - gridShape(1) ) + 1
        boundDir(1)    = 2
      end if 
    end if 

    ! Y
    if ( dimensionMask(2) .eq. 1 ) then 
      if ( ( gridIndexes(2) - this%matrixPositiveShape(1) ) .lt. 1 ) then 
        isBoundary(2)  = .true.
        boundLoc(2)    = 1 - gridIndexes(2) + this%matrixPositiveShape(1)
        boundDir(2)    = 1
      else if ( ( gridIndexes(2) + this%matrixPositiveShape(1) ) .gt. gridShape(2) ) then 
        isBoundary(2)  = .true.
        boundLoc(2)    = 2*this%matrixPositiveShape(1) + 1 &
            - ( gridIndexes(2) + this%matrixPositiveShape(1) - gridShape(2) ) + 1
        boundDir(2)    = 2
      end if 
    end if 

    ! Z
    if ( dimensionMask(3) .eq. 1 ) then 
      if ( ( gridIndexes(3) - this%matrixPositiveShape(3) ) .lt. 1 ) then 
        isBoundary(3)  = .true.
        boundLoc(3)    = 1 - gridIndexes(3) + this%matrixPositiveShape(3)
        boundDir(3)    = 1
      else if ( ( gridIndexes(3) + this%matrixPositiveShape(3) ) .gt. gridShape(3) ) then 
        isBoundary(3)  = .true.
        boundLoc(3)    = 2*this%matrixPositiveShape(3) + 1 &
            - ( gridIndexes(3) + this%matrixPositiveShape(3) - gridShape(3) ) + 1
        boundDir(3)    = 2
      end if 
    end if 

    call prTransposeXYMatrix( this ) 
    if ( any( isBoundary ) ) then 
      call prKernelReflection( this, gridIndexes, gridShape, &
                            xGridSpan, yGridSpan, zGridSpan, &
                      xKernelSpan, yKernelSpan, zKernelSpan, &
                             isBoundary, boundLoc, boundDir  )
    end if 


  end subroutine prComputeSpansBoundedTranspose



  subroutine prGenerateZeroPositiveGrid( this, zPXGrid, zPYGrid, zPZGrid  )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( KernelType ) :: this 
    integer, dimension(:,:,:), intent(inout) :: zPXGrid, zPYGrid, zPZGrid
    ! local
    integer :: nx, ny, nz
    integer :: i
    !------------------------------------------------------------------------------

    nx = this%matrixPositiveShape(1) 
    ny = this%matrixPositiveShape(2) 
    nz = this%matrixPositiveShape(3) 

    ! Positive octant
    zPXGrid = spread(  spread( [(i, i=0, nx)], 2, ny + 1 ), 3, nz + 1 )
    zPYGrid = spread(  spread( [(i, i=0, ny)], 1, nx + 1 ), 3, nz + 1 )
    zPZGrid = reshape( spread( [(i, i=0, nz)], 1, (nx + 1)*( ny + 1 ) ), &
                                              [ nx + 1, ny + 1, nz + 1 ] )

    return

  end subroutine prGenerateZeroPositiveGrid


  subroutine prUnfoldZeroPositiveMatrix( this, sourceZeroPositive, targetMatrix )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( KernelType ) :: this
    real(fp), dimension(:,:,:), intent(in)    :: sourceZeroPositive
    real(fp), dimension(:,:,:), intent(inout) :: targetMatrix
    ! local
    integer :: nx, ny, nz
    !------------------------------------------------------------------------------
    ! VERIFY WHAT HAPPENS WITH OCTANTS IN 1D/2D

    nx = this%matrixPositiveShape(1)
    ny = this%matrixPositiveShape(2)
    nz = this%matrixPositiveShape(3)

    targetMatrix( nx+1:2*nx+1 , ny+1:2*ny+1 , nz+1:2*nz+1 ) = sourceZeroPositive                                   ! Octant III
    targetMatrix( 1:nx        , ny+1:2*ny+1 , nz+1:2*nz+1 ) = sourceZeroPositive(nx+1:2:-1, :         , :        ) ! Octant OII
    targetMatrix( nx+1:2*nx+1 , 1:ny        , nz+1:2*nz+1 ) = sourceZeroPositive(:        , ny+1:2:-1 , :        ) ! Octant IOI
    targetMatrix( 1:nx        , 1:ny        , nz+1:2*nz+1 ) = sourceZeroPositive(nx+1:2:-1, ny+1:2:-1 , :        ) ! Octant OOI
    targetMatrix( nx+1:2*nx+1 , ny+1:2*ny+1 , 1:nz        ) = sourceZeroPositive(:        , :         , nz+1:2:-1) ! Octant IIO 
    targetMatrix( 1:nx        , ny+1:2*ny+1 , 1:nz        ) = sourceZeroPositive(nx+1:2:-1, :         , nz+1:2:-1) ! Octant OIO
    targetMatrix( nx+1:2*nx+1 , 1:ny        , 1:nz        ) = sourceZeroPositive(:        , ny+1:2:-1 , nz+1:2:-1) ! Octant IOO
    targetMatrix( 1:nx        , 1:ny        , 1:nz        ) = sourceZeroPositive(nx+1:2:-1, ny+1:2:-1 , nz+1:2:-1) ! Octant OOO


    return

  end subroutine prUnfoldZeroPositiveMatrix


  subroutine prSetupMatrix( this, bandwidth )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( KernelType ) :: this 
    real(fp), dimension(3), intent(in) :: bandwidth
    ! local
    integer, dimension(:,:,:), allocatable    :: zPXGrid, zPYGrid, zPZGrid
    !------------------------------------------------------------------------------

    ! Assign kernel properties
    this%bandwidth = bandwidth

    ! It is expected that for less than 3D, bandwidth of compressed 
    ! dimension is zero

    ! This means that matrixPositiveShape will have a zero in the compressed dimension
    where ( this%binSize .ne. fZERO ) 
      this%matrixPositiveShape = ceiling( this%matrixRange*this%bandwidth/this%binSize )
    elsewhere
      this%matrixPositiveShape = 0
    end where
    
    allocate( zPXGrid( this%matrixPositiveShape(1) + 1, this%matrixPositiveShape(2) + 1, this%matrixPositiveShape(3) + 1 ) )
    allocate( zPYGrid( this%matrixPositiveShape(1) + 1, this%matrixPositiveShape(2) + 1, this%matrixPositiveShape(3) + 1 ) )
    allocate( zPZGrid( this%matrixPositiveShape(1) + 1, this%matrixPositiveShape(2) + 1, this%matrixPositiveShape(3) + 1 ) )

    call this%GenerateZeroPositiveGrid( zPXGrid, zPYGrid, zPZGrid )
    call this%ComputeMatrix( zPXGrid, zPYGrid, zPZGrid )

    ! Necessary ?
    deallocate( zPXGrid )
    deallocate( zPYGrid )
    deallocate( zPZGrid )


    return

  end subroutine prSetupMatrix


  ! ComputeXYTranspose
  function prComputeXYTranspose( sourceMatrix ) result( transposedMatrix )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    real(fp), dimension(:,:,:), intent(in)  :: sourceMatrix
    real(fp), dimension(:,:,:), allocatable :: transposedMatrix
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


  subroutine prTransposeXYMatrix( this )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class( KernelType ), target :: this
    real(fp), dimension(:,:,:), pointer :: sourceMatrix
    real(fp), dimension(:,:,:), allocatable :: transposedMatrix
    ! local
    integer, dimension(3) :: sourceShape
    integer :: n
    !------------------------------------------------------------------------------
    
    sourceMatrix => this%matrix
    sourceShape = shape( sourceMatrix )
    allocate( transposedMatrix( sourceShape(2), sourceShape(1), sourceShape(3) ) )
    do n = 1, sourceShape(3)
      transposedMatrix(:,:,n) = transpose( sourceMatrix(:,:,n) )
    end do
 
    call move_alloc( transposedMatrix, this%matrix ) 

    return 

  end subroutine prTransposeXYMatrix



  ! Module methods
  subroutine InitializeKernelDimensions( inDimensionMask )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    integer, dimension(:), intent(in) :: inDimensionMask
    integer :: nd, dcount
    !------------------------------------------------------------------------------

    dimensionMask = inDimensionMask
    nDim          = sum( dimensionMask )
    fNDim         = real(nDim,fp)

    ! Save dim mask into dimensions 
    if ( allocated( dimensions ) ) deallocate( dimensions )
    allocate( dimensions( nDim  ) )
    dcount= 0
    do nd = 1, 3
      if ( dimensionMask(nd) .eq. 0 ) cycle
      dcount = dcount + 1
      dimensions(dcount) = nd
    end do

    ! The old way of detecting dims
    if ( nDim .eq. 1 ) then 
      do nd=1,3
        if ( dimensionMask(nd) .eq. 1 ) then 
          idDim1 = nd 
          exit
        end if 
      end do
    end if 
    if ( nDim .eq. 2 ) then
      idDim1 = 0
      idDim2 = 0 
      do nd=1,3
        if ( (dimensionMask(nd) .eq. 1) .and. (idDim1.eq.0) ) then 
          idDim1 = nd
        elseif ( dimensionMask(nd) .eq. 1 ) then 
          idDim2 = nd
          exit
        end if 
      end do
    end if 


  end subroutine InitializeKernelDimensions


  subroutine ResetKernelDimensions( )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    !------------------------------------------------------------------------------

    dimensionMask = (/1,1,1/)
    nDim          = 3 
    fNDim         = 3.0_fp 
    if( allocated( dimensions ) ) deallocate( dimensions ) 

  end subroutine ResetKernelDimensions


  ! KernelMultiGaussianW
  subroutine prComputeKernelWMatrix( this, zPXGrid, zPYGrid, zPZGrid )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none
    class( KernelMultiGaussianType ) :: this 
    integer, dimension(:,:,:), intent(in)   :: zPXGrid, zPYgrid, zPZGrid
    ! local 
    real(fp), dimension(3)                  :: hDelta
    real(fp), dimension(:,:,:), allocatable :: zeroPositiveMatrix
    integer  :: nx, ny, nz
    integer  :: nd
    real(fp) :: summatrix
    !------------------------------------------------------------------------------

    this%shouldIntegrateOne = .true.

    ! Initialize zeroPositiveMatrix
    nx = this%matrixPositiveShape(1)
    ny = this%matrixPositiveShape(2)
    nz = this%matrixPositiveShape(3)
    allocate( zeroPositiveMatrix( nx+1, ny+1, nz+1 ) ) 


    ! Kerne%matrix allocation ( consider doing this only if grid size changed )
    if ( allocated( this%matrix ) ) deallocate( this%matrix )
    allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )

    if ( all( this%bandwidth .eq. fZERO ) ) then 
      this%matrix = fZERO
      return
    end if

    ! Compute normalized smoothing bandwidth/Delta
    hDelta = fZERO
    where( this%binSize .ne. fZERO ) 
      hDelta = this%bandwidth/this%binSize
    end where

    ! Compute kernel
    zeroPositiveMatrix(:,:,:) = (0.5**fNDim)
    do nd=1,3
      if ( hDelta(nd) .le. fZERO ) cycle
      select case(nd)
      case(1)
          zeroPositiveMatrix = zeroPositiveMatrix*(         &
            erf( ( zPXGrid + 0.5 )/( hDelta(1)*sqrtTwo ) ) & 
          - erf( ( zPXGrid - 0.5 )/( hDelta(1)*sqrtTwo ) ) )
      case(2)
          zeroPositiveMatrix = zeroPositiveMatrix*(         &
            erf( ( zPYGrid + 0.5 )/( hDelta(2)*sqrtTwo ) ) & 
          - erf( ( zPYGrid - 0.5 )/( hDelta(2)*sqrtTwo ) ) )
      case(3)
          zeroPositiveMatrix = zeroPositiveMatrix*(         &
            erf( ( zPZGrid + 0.5 )/( hDelta(3)*sqrtTwo ) ) & 
          - erf( ( zPZGrid - 0.5 )/( hDelta(3)*sqrtTwo ) ) )
      end select
    end do 

    ! Unfold
    call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrix )

    ! Normalization correction
    summatrix = sum( this%matrix )
    if ( summatrix .ne. fZERO ) this%matrix = this%matrix/summatrix

    return

  end subroutine prComputeKernelWMatrix


  ! KernelSecondDerivativeX
  subroutine prComputeKernelVXMatrix( this, zPXGrid, zPYGrid, zPZGrid )
    !------------------------------------------------------------------------------
    ! Integrated X second derivative
    !
    ! Applies kernel corrections A.3 and A.4 in Sole-Mari et al. 2019
    !
    ! Returns the product Delta_i**2*V^i
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class( KernelSecondDerivativeXType )    :: this 
    integer, dimension(:,:,:), intent(in)   :: zPXGrid, zPYgrid, zPZGrid
    ! local
    real(fp), dimension(3)                  :: hDelta
    real(fp), dimension(:,:,:), allocatable :: zeroPositiveMatrix
    integer  :: nx, ny, nz
    integer  :: nd
    real(fp) :: aDenom, aNum, aCoeff, summatrixsq, normprod
    !------------------------------------------------------------------------------

    this%shouldIntegrateOne = .false.

    ! Initialize zeroPositiveMatrix
    nx = this%matrixPositiveShape(1)
    ny = this%matrixPositiveShape(2)
    nz = this%matrixPositiveShape(3)
    allocate( zeroPositiveMatrix( nx+1, ny+1, nz+1 ) ) 


    ! Kernel%matrix allocation ( consider doing this only if grid size changed )
    if ( allocated( this%matrix ) ) deallocate( this%matrix )
    allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )

    if ( all( this%bandwidth .eq. fZERO ) ) then 
      this%matrix = fZERO
      return
    end if

    ! Compute normalized smoothing bandwidth/Delta
    hDelta = fZERO
    where( this%binSize .ne. fZERO ) 
      hDelta = this%bandwidth/this%binSize
    end where


    ! Compute kernel
    zeroPositiveMatrix(:,:,:) = fONE
    do nd=1,3
      if ( hDelta(nd) .le. fZERO ) cycle
      select case(nd)
      case(1)
        zeroPositiveMatrix = zeroPositiveMatrix*( -fONE/( (fTWO**(fNDim-0.5))*sqrtPi*( hDelta(1)**fTHREE ) ) )*(&
                ( zPXGrid + 0.5 )*exp( -fONE*( ( zPXGrid + 0.5 )**fTWO )/( fTWO*( hDelta(1)**fTWO ) ) ) - &
                ( zPXGrid - 0.5 )*exp( -fONE*( ( zPXGrid - 0.5 )**fTWO )/( fTWO*( hDelta(1)**fTWO ) ) ) )
      case(2)
        zeroPositiveMatrix = zeroPositiveMatrix*( erf( ( zPYGrid + 0.5 )/( hDelta(2)*sqrtTwo ) ) - & 
                                                  erf( ( zPYGrid - 0.5 )/( hDelta(2)*sqrtTwo ) ) )
      case(3)
        zeroPositiveMatrix = zeroPositiveMatrix*( erf( ( zPZGrid + 0.5 )/( hDelta(3)*sqrtTwo ) ) - &
                                                  erf( ( zPZGrid - 0.5 )/( hDelta(3)*sqrtTwo ) ) )
      end select
    end do 

    ! Unfold
    call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrix )

    ! Kernel corrections
    aNum   = sum( this%matrix, mask=( this%matrix.gt.fZERO ) )
    aDenom = sum( this%matrix, mask=( this%matrix.lt.fZERO ) )
    aCoeff = fONE
    if ( aDenom .gt. fZERO ) aCoeff = -fONE*aNum/aDenom
    where ( this%matrix.gt.fZERO )
      this%matrix = aCoeff*this%matrix
    end where

    ! Health
    summatrixsq = sum( this%matrix**fTWO )
    if ( summatrixsq.le.fZERO ) then 
      this%matrix = fZERO
      return
    end if

    ! Kernel correction
    ! Note: the cell volume is implicit in the non-dimensional bandwidths,
    ! further considering that this kernel is Delta_i**2*V^i
    do nd=1,3
      if ( hDelta(nd) .le. fZERO ) cycle
      select case(nd)
      case(1)
        ! In order to handle limit cases, precompute the normalization product
        normprod  = ( fTWO**(fNDim + fTWO)*pi**(0.5*fNDim)*(hDelta(1)**fFIVE)*summatrixsq )
        if ( normprod.le.fZERO ) this%matrix = fZERO; return; 
        this%matrix = &
          this%matrix*sqrt(fTHREE/normprod)
      case(2)
        this%matrix = this%matrix/sqrt( hDelta(2) )
      case(3)
        this%matrix = this%matrix/sqrt( hDelta(3) )
      end select
    end do 

    return

  end subroutine prComputeKernelVXMatrix


  ! KernelSecondDerivativeY
  subroutine prComputeKernelVYMatrix( this, zPXGrid, zPYGrid, zPZGrid )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class( KernelSecondDerivativeYType )    :: this 
    integer, dimension(:,:,:), intent(in)   :: zPXGrid, zPYgrid, zPZGrid
    ! local
    real(fp), dimension(3)                  :: hDelta
    real(fp), dimension(:,:,:), allocatable :: zeroPositiveMatrix
    integer  :: nx, ny, nz
    integer  :: nd
    real(fp) :: aDenom, aNum, aCoeff, summatrixsq, normprod
    !------------------------------------------------------------------------------

    this%shouldIntegrateOne = .false.

    ! Initialize zeroPositiveMatrix
    nx = this%matrixPositiveShape(1)
    ny = this%matrixPositiveShape(2)
    nz = this%matrixPositiveShape(3)
    allocate( zeroPositiveMatrix( nx+1, ny+1, nz+1 ) ) 


    ! Kernel%matrix allocation ( consider doing this only if grid size changed )
    if ( allocated( this%matrix ) ) deallocate( this%matrix )
    allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )

    if ( all( this%bandwidth .eq. fZERO ) ) then 
      this%matrix = fZERO
      return
    end if

    ! Compute normalized smoothing bandwidth/Delta
    hDelta = fZERO
    where( this%binSize .ne. fZERO ) 
      hDelta = this%bandwidth/this%binSize
    end where


    ! Compute kernel
    zeroPositiveMatrix(:,:,:) = fONE
    do nd=1,3
      if ( hDelta(nd) .le. fZERO ) cycle
      select case(nd)
      case(1)
        zeroPositiveMatrix = zeroPositiveMatrix*( erf( ( zPXGrid + 0.5 )/( hDelta(1)*sqrtTwo ) ) - & 
                                                  erf( ( zPXGrid - 0.5 )/( hDelta(1)*sqrtTwo ) ) )
      case(2)
        zeroPositiveMatrix = zeroPositiveMatrix*( -fONE/( (fTWO**(fNDim-0.5))*sqrtPi*( hDelta(2)**fTHREE ) ) )*(&
                ( zPYGrid + 0.5 )*exp( -fONE*( ( zPYGrid + 0.5 )**fTWO )/( fTWO*( hDelta(2)**fTWO ) ) ) - &
                ( zPYGrid - 0.5 )*exp( -fONE*( ( zPYGrid - 0.5 )**fTWO )/( fTWO*( hDelta(2)**fTWO ) ) ) )
      case(3)
        zeroPositiveMatrix = zeroPositiveMatrix*( erf( ( zPZGrid + 0.5 )/( hDelta(3)*sqrtTwo ) ) - &
                                                  erf( ( zPZGrid - 0.5 )/( hDelta(3)*sqrtTwo ) ) )
      end select
    end do 

    ! Unfold
    call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrix )

    ! Kernel corrections
    aNum   = sum( this%matrix, mask=( this%matrix.gt.fZERO ) )
    aDenom = sum( this%matrix, mask=( this%matrix.lt.fZERO ) )
    aCoeff = fONE
    if ( aDenom .gt. fZERO ) aCoeff = -fONE*aNum/aDenom
    where ( this%matrix.gt.fZERO )
      this%matrix = aCoeff*this%matrix
    end where

    ! Health
    summatrixsq = sum( this%matrix**fTWO )
    if ( summatrixsq.le.fZERO ) then 
      this%matrix = fZERO
      return
    end if

    ! Kernel correction
    ! Note: the cell volume is implicit in the non-dimensional bandwidths,
    ! further considering that this kernel is Delta_i**2*V^i
    do nd=1,3
      if ( hDelta(nd) .le. fZERO ) cycle
      select case(nd)
      case(1)
        this%matrix = this%matrix/sqrt( hDelta(1) )
      case(2)
        ! In order to handle limit cases, precompute the normalization product
        normprod  = ( fTWO**(fNDim + fTWO)*pi**(0.5*fNDim)*(hDelta(2)**fFIVE)*summatrixsq )
        if ( normprod.le.fZERO ) this%matrix = fZERO; return; 
        this%matrix =& 
          this%matrix*sqrt(fTHREE/normprod)
      case(3)
        this%matrix = this%matrix/sqrt( hDelta(3) )
      end select
    end do 

    return

  end subroutine prComputeKernelVYMatrix


  ! KernelSecondDerivativeZ
  subroutine prComputeKernelVZMatrix( this, zPXGrid, zPYGrid, zPZGrid )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class( KernelSecondDerivativeZType )    :: this 
    integer, dimension(:,:,:), intent(in)   :: zPXGrid, zPYgrid, zPZGrid
    ! local
    real(fp), dimension(3)                  :: hDelta
    real(fp), dimension(:,:,:), allocatable :: zeroPositiveMatrix
    integer  :: nx, ny, nz
    integer  :: nd
    real(fp) :: aDenom, aNum, aCoeff, summatrixsq, normprod
    !------------------------------------------------------------------------------

    this%shouldIntegrateOne = .false.

    ! Initialize zeroPositiveMatrix
    nx = this%matrixPositiveShape(1)
    ny = this%matrixPositiveShape(2)
    nz = this%matrixPositiveShape(3)
    allocate( zeroPositiveMatrix( nx+1, ny+1, nz+1 ) ) 

    ! Kernel%matrix allocation ( consider doing this only if grid size changed )
    if ( allocated( this%matrix ) ) deallocate( this%matrix )
    allocate( this%matrix( 2*nx + 1, 2*ny + 1, 2*nz + 1 ) )

    if ( all( this%bandwidth .eq. fZERO ) ) then 
      this%matrix = fZERO
      return
    end if

    ! Compute normalized smoothing bandwidth/Delta
    hDelta = fZERO
    where( this%binSize .ne. fZERO ) 
      hDelta = this%bandwidth/this%binSize
    end where

    ! Compute kernel
    zeroPositiveMatrix(:,:,:) = fONE
    do nd=1,3
      if ( hDelta(nd) .le. fZERO ) cycle
      select case(nd)
      case(1)
        zeroPositiveMatrix = zeroPositiveMatrix*( erf( ( zPXGrid + 0.5 )/( hDelta(1)*sqrtTwo ) ) - & 
                                                  erf( ( zPXGrid - 0.5 )/( hDelta(1)*sqrtTwo ) ) )
      case(2)
        zeroPositiveMatrix = zeroPositiveMatrix*( erf( ( zPYGrid + 0.5 )/( hDelta(2)*sqrtTwo ) ) - &
                                                  erf( ( zPYGrid - 0.5 )/( hDelta(2)*sqrtTwo ) ) )
      case(3)
        zeroPositiveMatrix = zeroPositiveMatrix*( -fONE/( ( fTWO**(fNDim-0.5) )*sqrtPi*( hDelta(3)**fTHREE ) ) )*(&
                ( zPZGrid + 0.5 )*exp( -fONE*( ( zPZGrid + 0.5 )**fTWO )/( fTWO*( hDelta(3)**fTWO ) ) ) - &
                ( zPZGrid - 0.5 )*exp( -fONE*( ( zPZGrid - 0.5 )**fTWO )/( fTWO*( hDelta(3)**fTWO ) ) ) )
      end select
    end do 

    ! Unfold
    call this%UnfoldZeroPositiveMatrix( zeroPositiveMatrix, this%matrix )

    ! Kernel corrections
    aNum   = sum( this%matrix, mask=( this%matrix.gt.fZERO ) )
    aDenom = sum( this%matrix, mask=( this%matrix.lt.fZERO ) )
    aCoeff = fONE
    if ( aDenom .gt. fZERO ) aCoeff = -fONE*aNum/aDenom
    where ( this%matrix.gt.fZERO )
      this%matrix = aCoeff*this%matrix
    end where

    ! Health
    summatrixsq = sum( this%matrix**fTWO )
    if ( summatrixsq.le.fZERO ) then 
      this%matrix = fZERO
      return
    end if

    ! Kernel correction
    ! Note: the cell volume is implicit in the non-dimensional bandwidths,
    ! further considering that this kernel is Delta_i**2*V^i
    do nd=1,3
      if ( hDelta(nd) .le. fZERO ) cycle
      select case(nd)
      case(1)
        this%matrix = this%matrix/sqrt( hDelta(1) )
      case(2)
        this%matrix = this%matrix/sqrt( hDelta(2) )
      case(3)
        ! In order to handle limit cases, precompute the normalization product
        normprod  = ( fTWO**(fNDim + fTWO)*pi**(0.5*fNDim)*(hDelta(3)**fFIVE)*summatrixsq )
        if ( normprod.le.fZERO ) this%matrix = fZERO; return; 
        this%matrix =& 
          this%matrix*sqrt(fTHREE/normprod)
      end select
    end do 


    return


  end subroutine prComputeKernelVZMatrix


end module KernelMultiGaussianModule
