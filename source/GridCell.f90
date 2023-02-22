module GridCellModule

    use KernelMultiGaussianModule, only : KernelMultiGaussianType, &
                                      KernelSecondDerivativeXType, &
                                      KernelSecondDerivativeYType, &
                                      KernelSecondDerivativeZType, &
                                      KernelType
    implicit none

    type GridCellType
    
        integer, dimension(3) :: id
        logical               :: convergence = .false.

        ! Kernel matrices
        doubleprecision, dimension(:,:,:), allocatable :: kernelSigmaMatrix
        doubleprecision, dimension(:,:,:), allocatable :: kernelMatrix
        doubleprecision, dimension(:,:,:), allocatable :: kernelSD1Matrix
        doubleprecision, dimension(:,:,:), allocatable :: kernelSD2Matrix
        doubleprecision, dimension(:,:,:), allocatable :: kernelSD3Matrix

        ! Kernel indexes
        integer, dimension(3) :: kernelDBIndexes      = 0
        integer, dimension(3) :: kernelSigmaDBIndexes = 0
        integer, dimension(3) :: kernelSDDBIndexes    = 0

        ! Flat db stuff
        integer, dimension(2) :: kernelDBFlatIndexes      = 0
        integer, dimension(2) :: kernelSigmaDBFlatIndexes = 0
        logical               :: transposeKernel = .false.
        logical               :: transposeKernelSigma = .false.


        ! For the cases with zero support
        logical               :: skipKernelSigma = .false.


        ! Spans
        integer, dimension(2) :: kernelXGSpan = 0
        integer, dimension(2) :: kernelYGSpan = 0
        integer, dimension(2) :: kernelZGSpan = 0
        integer, dimension(2) :: kernelXMSpan = 0
        integer, dimension(2) :: kernelYMSpan = 0
        integer, dimension(2) :: kernelZMSpan = 0

        integer, dimension(2) :: kernelSigmaXGSpan = 0
        integer, dimension(2) :: kernelSigmaYGSpan = 0
        integer, dimension(2) :: kernelSigmaZGSpan = 0
        integer, dimension(2) :: kernelSigmaXMSpan = 0
        integer, dimension(2) :: kernelSigmaYMSpan = 0
        integer, dimension(2) :: kernelSigmaZMSpan = 0

        integer, dimension(2) :: kernelSDXGSpan = 0
        integer, dimension(2) :: kernelSDYGSpan = 0
        integer, dimension(2) :: kernelSDZGSpan = 0
        integer, dimension(2) :: kernelSDXMSpan = 0
        integer, dimension(2) :: kernelSDYMSpan = 0
        integer, dimension(2) :: kernelSDZMSpan = 0

        integer, dimension(2) :: kernelSD1XGSpan = 0
        integer, dimension(2) :: kernelSD1YGSpan = 0
        integer, dimension(2) :: kernelSD1ZGSpan = 0
        integer, dimension(2) :: kernelSD1XMSpan = 0
        integer, dimension(2) :: kernelSD1YMSpan = 0
        integer, dimension(2) :: kernelSD1ZMSpan = 0

        integer, dimension(2) :: kernelSD2XGSpan = 0
        integer, dimension(2) :: kernelSD2YGSpan = 0
        integer, dimension(2) :: kernelSD2ZGSpan = 0
        integer, dimension(2) :: kernelSD2XMSpan = 0
        integer, dimension(2) :: kernelSD2YMSpan = 0
        integer, dimension(2) :: kernelSD2ZMSpan = 0

        integer, dimension(2) :: kernelSD3XGSpan = 0
        integer, dimension(2) :: kernelSD3YGSpan = 0
        integer, dimension(2) :: kernelSD3ZGSpan = 0
        integer, dimension(2) :: kernelSD3XMSpan = 0
        integer, dimension(2) :: kernelSD3YMSpan = 0
        integer, dimension(2) :: kernelSD3ZMSpan = 0

    contains

        procedure :: Initialize => prInitialize
        procedure :: Reset      => prReset

    end type GridCellType


contains


    subroutine prInitialize( this, id )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridCellType ) :: this
        integer, dimension(3), intent(in) :: id
        !------------------------------------------------------------------------------

        this%id = id

    end subroutine prInitialize


    subroutine prReset( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none
        class( GridCellType ) :: this
        !------------------------------------------------------------------------------

        ! Kernel indexes
        this%kernelDBIndexes          = 0
        this%kernelSigmaDBIndexes     = 0
        this%kernelSDDBIndexes        = 0
        this%kernelDBFlatIndexes      = 0
        this%kernelSigmaDBFlatIndexes = 0
        this%transposeKernel          = .false.
        this%transposeKernelSigma     = .false.
        this%skipKernelSigma          = .false.
        this%kernelXGSpan = 0
        this%kernelYGSpan = 0
        this%kernelZGSpan = 0
        this%kernelXMSpan = 0
        this%kernelYMSpan = 0
        this%kernelZMSpan = 0
        this%kernelSigmaXGSpan = 0
        this%kernelSigmaYGSpan = 0
        this%kernelSigmaZGSpan = 0
        this%kernelSigmaXMSpan = 0
        this%kernelSigmaYMSpan = 0
        this%kernelSigmaZMSpan = 0
        this%kernelSDXGSpan = 0
        this%kernelSDYGSpan = 0
        this%kernelSDZGSpan = 0
        this%kernelSDXMSpan = 0
        this%kernelSDYMSpan = 0
        this%kernelSDZMSpan = 0

    end subroutine





end module GridCellModule
