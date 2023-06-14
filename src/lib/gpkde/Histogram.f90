module HistogramModule
  !------------------------------------------------------------------------------
  ! Histogram on a regular grid
  !  - Count elements as individual points
  !  - Count weighted elements 
  !------------------------------------------------------------------------------
  use PrecisionModule, only : fp
  use ConstantsModule, only : fZERO, fONE, fTWO, fTHREE
  use iso_c_binding
  !------------------------------------------------------------------------------
  implicit none

  ! Set default access status to private
  private

  type, public :: HistogramType
    ! Properties
    real(fp), dimension(:,:,:), allocatable :: counts 
    real(fp), dimension(:,:,:), allocatable :: wcounts 
    integer , dimension(3)                  :: domainGridSize
    integer , dimension(3)                  :: gridSize
    integer , dimension(:), pointer         :: nBins
    real(fp), dimension(3)                  :: binSize
    real(fp)                                :: binVolume
    real(fp)                                :: binDistance
    integer , dimension(:,:), allocatable   :: activeBinIds
    integer                                 :: nActiveBins
    integer , dimension(:,:), allocatable   :: boundingBoxBinIds
    integer                                 :: nBBoxBins 
    real(fp), dimension(3)                  :: domainOrigin ! of the reconstruction grid 
    real(fp), dimension(3)                  :: gridOrigin   ! while allocating from dataset
    real(fp), dimension(:), pointer         :: origin       ! point to the proper origin for indexes 
    logical                                 :: adaptGridToCoords
    integer, dimension(:), allocatable      :: dimensions
    integer                                 :: nDim
    integer                                 :: nPoints 
    real(fp)                                :: nEffective 
    real(fp)                                :: totalMass, effectiveMass 
    logical                                 :: isWeighted
    real(fp)                                :: maxCount
    real(fp)                                :: maxRawDensity
    integer                                 :: effectiveWeightFormat = 0
  ! HistogramType contains
  contains
    ! Procedures
    procedure :: Initialize                       => prInitialize
    procedure :: Reset                            => prReset
    procedure :: ComputeCounts                    => prComputeCounts
    procedure :: ComputeCountsWeighted            => prComputeCountsWeighted
    procedure :: ComputeCountsAndWeights          => prComputeCountsAndWeights
    procedure :: ComputeEffectiveCountsAndWeights => prComputeEffectiveCountsAndWeights
    procedure :: ComputeActiveBinIds              => prComputeActiveBinIds
    procedure :: ComputeActiveBinIdsSliced        => prComputeActiveBinIdsSliced
    procedure :: EstimateBinSizeFD                => prEstimateBinSizeFD
    procedure :: EstimateBinSizeScott             => prEstimateBinSizeScott
  end type 


  interface
    subroutine prSortDblArray1D( array, elem_count, elem_size, compare) bind(C,name="qsort")
      !--------------------------------------------------------------------------------------------- 
      ! Interface sort to standard C library qsort. 
      ! Adapted from:
      ! https://stackoverflow.com/questions/20941575/sorting-in-fortran-undefined-reference-to-qsort
      !--------------------------------------------------------------------------------------------- 
      import
      real(c_double),intent(inout) :: array(*)
      integer(c_size_t),value      :: elem_count
      integer(c_size_t),value      :: elem_size
      type(c_funptr),value         :: compare !int(*compare)(const void *, const void *)
    end subroutine prSortDblArray1D 
  end interface


! HistogramModule contains
contains

  subroutine prInitialize( this, & 
        domainGridSize, binSize, & 
                   domainOrigin, & 
              adaptGridToCoords  )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class(HistogramType), target                 :: this
    integer , dimension(3), intent(in)           :: domainGridSize
    real(fp), dimension(3), intent(in)           :: binSize
    real(fp), dimension(3), intent(in), optional :: domainOrigin
    logical ,               intent(in), optional :: adaptGridToCoords
    integer :: nd, dcount
    !------------------------------------------------------------------------------

    ! Stop if all bin sizes are wrong
    if ( all( binSize .lt. fZERO ) ) then 
      write(*,*)'Error: while initializing Histogram, all binSizes are .lt. 0. Stop.'
      stop 
    end if 

    ! Stop if any nBins .lt. 1
    if ( any( domainGridSize .lt. 1 ) ) then 
      write(*,*)'Error: while initializing Histogram, some domainGridSize .lt. 1. Stop.'
      stop
    end if

    ! Allocate grid with nBins ? 
    if( present(adaptGridToCoords) ) then 
      this%adaptGridToCoords = adaptGridToCoords
    else
      this%adaptGridToCoords = .false.
    end if

    ! Assign dim properties
    this%nDim = sum((/1,1,1/), mask=(binSize.gt.fZERO))
    if ( this%nDim .le. 0 ) then 
      write(*,*)'Error: while initializing Histogram, nDim .le. 0. Stop.'
      stop
    end if 

    ! Save dim mask into dimensions
    ! Notice that for histogram purposes, a dimension 
    ! will be marked as inactive IF the binSize in said
    ! dimension is 0. This addresses those cases
    ! where, for example, a 2D reconstruction is made over 
    ! a slice of a 3D distribution of particles. 
    if ( allocated( this%dimensions ) ) deallocate( this%dimensions ) 
    allocate( this%dimensions( this%nDim  ) )
    dcount= 0
    do nd = 1, 3
      if ( binSize(nd) .le. fZERO ) cycle
      dcount = dcount + 1
      this%dimensions(dcount) = nd
    end do 

    ! domainOrigin
    if ( present( domainOrigin ) ) then 
      this%domainOrigin = domainOrigin
    else 
      this%domainOrigin = fZERO 
    end if

    ! Initialize variables
    this%domainGridSize = domainGridSize
    this%binSize        = binSize
    this%binVolume      = product( binSize, mask=(binSize.gt.fZERO) ) 
    this%binDistance    = ( this%binVolume )**(fONE/real(this%nDim,fp))

    ! Allocate and initialize histogram counts, 
    ! if not adapting to the particle distribution follows the domain grid. 
    if ( .not. this%adaptGridToCoords ) then 
      if ( allocated( this%counts ) ) deallocate( this%counts ) 
      allocate( this%counts( domainGridSize(1), domainGridSize(2), domainGridSize(3) ) )
      this%counts   = fZERO
      this%origin   => this%domainOrigin
      this%nBins    => this%domainGridSize
      this%gridSize = this%domainGridSize
    end if 


  end subroutine prInitialize


  subroutine prReset( this )
    !------------------------------------------------------------------------------
    ! 
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class(HistogramType) :: this
    !------------------------------------------------------------------------------

    if( allocated(this%counts )) deallocate( this%counts  )
    if( allocated(this%wcounts)) deallocate( this%wcounts )
    this%domainGridSize = 0
    this%gridSize       = 0
    this%nBins          => null()
    this%binSize        = fZERO
    this%binVolume      = fZERO
    this%binDistance    = fZERO
    if( allocated(this%activeBinIds) ) deallocate( this%activeBinIds ) 
    this%nActiveBins = 0 
    if( allocated(this%boundingBoxBinIds) ) deallocate( this%boundingBoxBinIds )
    this%nBBoxBins = 0
    this%domainOrigin = fZERO
    this%gridOrigin = fZERO
    this%origin => null()
    this%adaptGridToCoords = .false.
    if( allocated( this%dimensions) ) deallocate( this%dimensions) 
    this%nDim = 3
    this%nPoints = 0
    this%nEffective = fZERO
    this%totalMass = fZERO
    this%effectiveMass = fZERO
    this%isWeighted = .false.
    this%maxCount = 0
    this%maxRawDensity = fZERO
    this%effectiveWeightFormat = 0

  end subroutine prReset


  subroutine prComputeCounts( this, dataPoints, exact )
    !------------------------------------------------------------------------------
    ! 
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class(HistogramType) :: this
    real(fp), dimension(:,:), intent(in) :: dataPoints
    logical , intent(in), optional       :: exact 
    integer , dimension(2)               :: nPointsShape
    integer , dimension(3)               :: gridIndexes
    integer :: np, nd, did
    logical :: inside
    integer :: exactIndex 
    !------------------------------------------------------------------------------

    exactIndex = 1
    if( present(exact) ) then 
      if ( exact ) exactIndex = 0
    end if

    ! Reset counts
    this%counts  = fZERO
    nPointsShape = shape(dataPoints) ! (nParticles,3)

    do np = 1, nPointsShape(1)
      ! Initialize point
      inside      = .true.
      gridIndexes = 1
      ! Compute gridIndex and verify  
      ! only for active dimensions if inside or not.
      ! Histogram considers active dimensions those where
      ! binSize is non zero.
      do nd = 1,this%nDim
        did = this%dimensions(nd)
        gridIndexes(did) = floor(( dataPoints(np,did) - this%origin(did) )/this%binSize(did)) + exactIndex
        if( (gridIndexes(did) .gt. this%nBins(did)) .or.&
            (gridIndexes(did) .le. 0) ) then
          inside = .false. 
          exit
        end if 
      end do
      if ( .not. inside ) cycle

      ! Increase counter
      this%counts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) = &
          this%counts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) + fONE
    end do

    this%totalMass     = sum(this%counts)
    this%nPoints       = int(this%totalMass)
    this%nEffective    = this%nPoints
    this%effectiveMass = fONE
    this%isWeighted    = .false.
    this%maxCount      = maxval(this%counts)
    this%maxRawDensity = this%maxCount/this%binVolume


  end subroutine prComputeCounts


  subroutine prComputeCountsWeighted( this, dataPoints, weights, exact )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class(HistogramType) :: this
    real(fp), dimension(:,:), intent(in) :: dataPoints
    real(fp), dimension(:), intent(in)   :: weights
    logical , intent(in), optional       :: exact 
    integer , dimension(2)               :: nPointsShape
    integer , dimension(3)               :: gridIndexes
    integer :: np, nd, did
    logical :: inside
    integer :: exactIndex 
    !------------------------------------------------------------------------------

    exactIndex = 1
    if( present(exact) ) then 
      if ( exact ) exactIndex = 0
    end if

    ! Reset counts
    this%counts  = fZERO
    nPointsShape = shape(dataPoints) ! (nParticles,3)

    do np = 1, nPointsShape(1)
      ! Initialize point
      inside      = .true.
      gridIndexes = 1
      ! Compute gridIndex and verify  
      ! only for active dimensions if inside or not.
      ! Histogram considers active dimensions those where
      ! binSize is non zero.
      do nd = 1,this%nDim
        did = this%dimensions(nd)
        gridIndexes(did) = floor(( dataPoints(np,did) - this%origin(did) )/this%binSize(did)) + exactIndex
        if( (gridIndexes(did) .gt. this%nBins(did)) .or.&
            (gridIndexes(did) .le. 0) ) then
          inside = .false. 
          exit
        end if 
      end do
      if ( .not. inside ) cycle

      ! Increase counter
      this%counts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) = &
          this%counts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) + weights(np)

    end do 
 
    ! Get total mass 
    this%totalMass = sum(weights)
    ! Compute the effective number of particles 
    ! and mass according to effectiveWeightFormat
    select case(this%effectiveWeightFormat)
    case(1)
      ! 1: using the simple average mass
      this%nEffective = nPointsShape(1)
      this%effectiveMass = this%totalMass/this%nEffective
    case default
      ! Defaults to Kish (1965,1992)
      this%nEffective    = this%totalMass**2/sum(weights**2)
      this%effectiveMass = this%totalMass/this%nEffective
    end select
    ! Transform mass counts into n counts
    this%counts        = this%counts/this%effectiveMass
    this%nPoints       = size(weights)
    this%isWeighted    = .true.
    this%maxCount      = maxval(this%counts)
    this%maxRawDensity = this%maxCount/this%binVolume


  end subroutine prComputeCountsWeighted


  subroutine prComputeCountsAndWeights( this, dataPoints, weights, exact )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class(HistogramType) :: this
    real(fp), dimension(:,:), intent(in) :: dataPoints
    real(fp), dimension(:), intent(in)   :: weights
    logical , intent(in), optional       :: exact 
    integer , dimension(2)               :: nPointsShape
    integer , dimension(3)               :: gridIndexes
    integer :: np, nd, did
    logical :: inside
    integer :: exactIndex 
    !------------------------------------------------------------------------------

    exactIndex = 1
    if( present(exact) ) then 
      if ( exact ) exactIndex = 0
    end if

    ! Reset counts
    this%counts  = fZERO
    this%wcounts = fZERO
    nPointsShape = shape(dataPoints) ! (nParticles,3)

    do np = 1, nPointsShape(1)
      ! Initialize point
      inside      = .true.
      gridIndexes = 1
      ! Compute gridIndex and verify  
      ! only for active dimensions if inside or not.
      ! Histogram considers active dimensions those where
      ! binSize is non zero.
      do nd = 1,this%nDim
        did = this%dimensions(nd)
        gridIndexes(did) = floor(( dataPoints(np,did) - this%origin(did) )/this%binSize(did)) + exactIndex
        if( (gridIndexes(did) .gt. this%nBins(did)) .or.&
            (gridIndexes(did) .le. 0) ) then
          inside = .false. 
          exit
        end if 
      end do
      if ( .not. inside ) cycle

      ! Increase counter
      this%counts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) = &
          this%counts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) + fONE
      this%wcounts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) = &
          this%wcounts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) + weights(np)

    end do 
 
    ! Histogram metric, no transformation of counts 
    this%totalMass     = sum(weights)
    this%nEffective    = nPointsShape(1)
    this%nPoints       = size(weights)
    this%effectiveMass = fONE
    this%isWeighted    = .true.
    this%maxCount      = maxval(this%counts)
    this%maxRawDensity = this%maxCount/this%binVolume


  end subroutine prComputeCountsAndWeights


  subroutine prComputeEffectiveCountsAndWeights( this, dataPoints, weights, exact )
    !------------------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class(HistogramType) :: this
    real(fp), dimension(:,:), intent(in) :: dataPoints
    real(fp), dimension(:), intent(in)   :: weights
    logical , intent(in), optional       :: exact 
    integer , dimension(2)               :: nPointsShape
    integer , dimension(3)               :: gridIndexes
    integer :: np, nd, did
    logical :: inside
    integer :: exactIndex 
    !------------------------------------------------------------------------------

    exactIndex = 1
    if( present(exact) ) then 
      if ( exact ) exactIndex = 0
    end if

    ! Reset counts
    this%counts  = fZERO
    this%wcounts = fZERO
    nPointsShape = shape(dataPoints) ! (nParticles,3)

    do np = 1, nPointsShape(1)
      ! Initialize point
      inside      = .true.
      gridIndexes = 1
      ! Compute gridIndex and verify  
      ! only for active dimensions if inside or not.
      ! Histogram considers active dimensions those where
      ! binSize is non zero.
      do nd = 1,this%nDim
        did = this%dimensions(nd)
        gridIndexes(did) = floor(( dataPoints(np,did) - this%origin(did) )/this%binSize(did)) + exactIndex
        if( (gridIndexes(did) .gt. this%nBins(did)) .or.&
            (gridIndexes(did) .le. 0) ) then
          inside = .false. 
          exit
        end if 
      end do
      if ( .not. inside ) cycle

      ! Increase counter
      this%counts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) = &
          this%counts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) + weights(np)**2
      this%wcounts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) = &
          this%wcounts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) + weights(np)

    end do 
 
    ! Histogram metric, transform counts into the 
    ! nEffective for each histogram cell M**2/sum mp**2
    where( this%counts.ne.fZERO )
      this%counts=this%wcounts**2/this%counts
    end where
    this%totalMass     = sum(weights)
    this%nEffective    = sum(this%counts)
    this%nPoints       = size(weights)
    this%effectiveMass = fONE
    this%isWeighted    = .true.
    this%maxCount      = maxval(this%counts)
    this%maxRawDensity = this%maxCount/this%binVolume


  end subroutine prComputeEffectiveCountsAndWeights


  subroutine prComputeActiveBinIds( this )
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! Specifications 
    !------------------------------------------------------------------------------
    implicit none 
    class(HistogramType) :: this
    integer              :: ix, iy, iz
    integer              :: icount = 1 
    !------------------------------------------------------------------------------

    ! Reset icount
    icount = 1

    this%nActiveBins = count( this%counts/=fZERO )

    ! Reallocate activeBinIds to new size 
    if ( allocated( this%activeBinIds ) )  deallocate( this%activeBinIds )
    allocate( this%activeBinIds( 3, this%nActiveBins ) )
    
    ! Following column-major nesting
    ! This could be in parallel with OpenMP (?)
    do iz = 1, this%nBins(3)
      do iy = 1, this%nBins(2)
        do ix = 1, this%nBins(1)
          if ( this%counts( ix, iy, iz ) .eq. fZERO ) cycle
          this%activeBinIds( :, icount ) = [ ix, iy, iz ]
          icount = icount + 1
        end do
      end do
    end do


  end subroutine prComputeActiveBinIds


  subroutine prComputeActiveBinIdsSliced( this, slicedDimension, sliceIndex )
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  ! Specifications 
  !------------------------------------------------------------------------------
  implicit none 
  class(HistogramType) :: this
  integer, intent(in)  :: slicedDimension
  integer, intent(in)  :: sliceIndex
  integer              :: ix, iy, iz
  integer              :: icount = 1 
  !------------------------------------------------------------------------------

    ! Reset icount
    icount = 1

    if ( ( slicedDimension.lt.1 ).or.( slicedDimension.gt.3 ) ) then 
      write(*,*)'Error: while computing sliced active bins, invalid sliced dimension.'
      stop 
    end if
    if ( ( sliceIndex.le.0 ).or.( sliceIndex.gt.this%nBins(slicedDimension) ) ) then 
      write(*,*)'Error: while computing sliced active bins, invalid slice index.'
      stop 
    end if

    ! Compute active bins on a slice
    select case(slicedDimension)
    case(3)
      this%nActiveBins = count( this%counts(:,:,sliceIndex)/=fZERO )
      ! Reallocate activeBinIds to new size 
      if ( allocated( this%activeBinIds ) )  deallocate( this%activeBinIds )
      allocate( this%activeBinIds( 3, this%nActiveBins ) )
      ! Following column-major nesting
      do iy = 1, this%nBins(2)
        do ix = 1, this%nBins(1)
          if ( this%counts( ix, iy, sliceIndex ) .eq. fZERO ) cycle
          this%activeBinIds( :, icount ) = [ ix, iy, 1 ]
          icount = icount + 1
        end do
      end do
    case(2)
      this%nActiveBins = count( this%counts(:,sliceIndex,:)/=fZERO )
      ! Reallocate activeBinIds to new size 
      if ( allocated( this%activeBinIds ) )  deallocate( this%activeBinIds )
      allocate( this%activeBinIds( 3, this%nActiveBins ) )
      ! Following column-major nesting
      do iz = 1, this%nBins(3)
        do ix = 1, this%nBins(1)
          if ( this%counts( ix, sliceIndex, iz ) .eq. fZERO ) cycle
          this%activeBinIds( :, icount ) = [ ix, 1, iz ]
          icount = icount + 1
        end do
      end do
    case(1)
      this%nActiveBins = count( this%counts(sliceIndex,:,:)/=fZERO )
      ! Reallocate activeBinIds to new size 
      if ( allocated( this%activeBinIds ) )  deallocate( this%activeBinIds )
      allocate( this%activeBinIds( 3, this%nActiveBins ) )
      ! Following column-major nesting
      do iz = 1, this%nBins(3)
        do iy = 1, this%nBins(2)
          if ( this%counts( sliceIndex, iy, iz ) .eq. fZERO ) cycle
          this%activeBinIds( :, icount ) = [ 1, iy, iz ]
          icount = icount + 1
        end do
      end do
    end select


  end subroutine prComputeActiveBinIdsSliced


  subroutine prEstimateBinSizeScott( this, dataPoints, binSize ) 
  !------------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------------
  ! Specifications 
  !------------------------------------------------------------------------------
  implicit none 
  class(HistogramType) :: this
  real(fp), dimension(:,:), intent(in)  :: dataPoints ! (npoints, 3)
  real(fp), dimension(3), intent(inout) :: binSize    ! for 3D but now only use first
  real(fp) :: avg, var
  integer  :: nPoints
  integer, dimension(2) :: nPointsShape
  integer  :: n 
  !------------------------------------------------------------------------------
   
    nPointsShape = shape(dataPoints) 
    nPoints = nPointsShape(1)

    avg = sum(dataPoints(:,1))/real(nPoints,fp)
    var = fZERO
 
    do n =1, nPoints
      var = var + (dataPoints(n,1) - avg)**2d0
    end do 
    var = var/real(nPoints,fp)

    binSize(1) = 3.49*sqrt(var)/( (real(nPoints,fp))**(fONE/fTHREE) )


  end subroutine prEstimateBinSizeScott


  subroutine prEstimateBinSizeFD( this, dataPoints, binSize ) 
  !------------------------------------------------------------------------------
  ! Freedman Diaconis
  !------------------------------------------------------------------------------
  ! Specifications 
  !------------------------------------------------------------------------------
  use iso_c_binding
  implicit none 
  class(HistogramType) :: this
  real(fp), dimension(:,:), intent(in) :: dataPoints ! (npoints, 3)
  real(fp), dimension(3), intent(inout) :: binSize    ! for 3D but now only use first
  !real(fp) :: avg, var
  integer  :: nPoints
  integer, dimension(2) :: nPointsShape
  !integer  :: n 
  real(fp), dimension(:), allocatable :: array
  integer  :: residue
  real(fp) :: residuefp
  real(fp) :: Q1, Q3, IQR
  !------------------------------------------------------------------------------
    
    nPointsShape = shape(dataPoints) 
    nPoints = nPointsShape(1)
    if ( allocated(array) ) deallocate( array ) 
    allocate( array(nPoints) )
    array(:) = dataPoints(:,1)

    ! Preprocessor flags for type selection !
    call prSortDblArray1D( array, int(nPoints,c_size_t), int(fp,c_size_t), c_funloc( comparedbl ) ) 

    residue = modulo( nPoints, 2 )
    if ( residue .eq. 0 ) then
      ! check the number of elements in half the list is even
      ! if it is, the median if Q1 is the average between n/4 and n/4+1 
      residuefp = modulo( real(nPoints,fp)/fTWO, fTWO)
      if ( residuefp.eq.0 ) then
        Q1 = ( array(int(nPoints/4)) + array(int(nPoints/4)+1) )/fTWO
        Q3 = ( array(int(nPoints/2) + int(nPoints/4)) + array(int(nPoints/2) +int(nPoints/4)+1) )/fTWO
      else
      ! if odd, the median is the center value      
        Q1 = array(int(nPoints/4)+1)
        Q3 = array(nPoints/2+int(nPoints/4)+1)
      end if 
    else
      ! if the length of the original list is odd..
      ! need to verify the same as before but for nPoints - 1. 
      residuefp = modulo( real((nPoints-1),fp)/fTWO, fTWO)
      if ( residuefp.eq.0 ) then
        Q1 = ( array(int((nPoints-1)/4)) + array(int((nPoints-1)/4)+1) )/fTWO
        Q3 = ( array(int(nPoints/2) + 1 + int((nPoints-1)/4)     ) + &
               array(int(nPoints/2) + 1 + int((nPoints-1)/4) + 1 ) )/fTWO
      else
      ! if odd, the median is the center value      
        Q1 = array( int((nPoints-1)/4) + 1 )
        Q3 = array( int(nPoints/2) + 1 + int((nPoints-1)/4) + 1 )
      end if 
    end if
    
    ! The interquartile range
    IQR = Q3-Q1 

    ! The estimated bin size
    binSize(1) = fTWO*IQR/(real(nPoints,fp)**(fONE/fTHREE))


  end subroutine prEstimateBinSizeFD


  integer(c_int) function comparedbl( a, b ) bind(C)
    use iso_c_binding
    real(c_double) a, b
    if ( a .lt. b ) comparedbl = -1
    if ( a .eq. b ) comparedbl = 0
    if ( a .gt. b ) comparedbl = 1
  end function comparedbl


end module HistogramModule

