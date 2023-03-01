module HistogramModule
    !------------------------------------------------------------------------------
    ! Histogram related operations
    !  - Count elements on a regular grid
    !  - Count weighted elements on a regular grid 
    !------------------------------------------------------------------------------
    implicit none

    ! Set default access status to private
    private

    type, public :: HistogramType

        ! Properties
        doubleprecision, dimension(:,:,:), allocatable :: counts 
        integer, dimension(:,:,:), allocatable     :: ncounts 
        doubleprecision, dimension(:,:,:), allocatable     :: avgmbin
        integer, dimension(3)                      :: nBins
        doubleprecision, dimension(3)              :: binSize
        doubleprecision                            :: binVolume
        doubleprecision                            :: binDistance
        integer, dimension(:,:), allocatable       :: activeBinIds
        integer                                    :: nActiveBins
        integer, dimension(:,:), allocatable       :: boundingBoxBinIds
        integer                                    :: nBBoxBins 
        doubleprecision, dimension(3)              :: domainOrigin ! of the reconstruction grid 
        integer, dimension(:), allocatable         :: dimensions
        integer                                    :: nDim

    contains

        ! Procedures
        procedure :: Initialize    => prInitialize
        procedure :: Reset         => prReset
        procedure :: ComputeCounts => prComputeCounts
        procedure :: ComputeCountsWeighted => prComputeCountsWeighted
        procedure :: Export        => prExport
        procedure :: ComputeActiveBinIds => prComputeActiveBinIds
        procedure :: ComputeBoundingBox  => prComputeBoundingBox

    end type 



contains


    subroutine prInitialize( this, nBins, binSize, dimensionMask, domainOrigin )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(HistogramType)                      :: this
        integer, dimension(3), intent(in)         :: nBins
        doubleprecision, dimension(3), intent(in) :: binSize
        integer, dimension(3), optional           :: dimensionMask
        integer, dimension(3)                     :: locDimensionMask
        doubleprecision, dimension(3), optional   :: domainOrigin
        integer :: nd, dcount
        !------------------------------------------------------------------------------

        ! Stop if all bin sizes are wrong
        if ( all( binSize .lt. 0d0 ) ) then 
          write(*,*)'Error while initializing Histogram, all binSizes are .lt. 0d0. Stop.'
          stop 
        end if 

        ! Stop if any nBins .lt. 1
        if ( any( nBins .lt. 1 ) ) then 
          write(*,*)'Error while initializing Histogram, some nBins .lt. 1. Stop.'
          stop
        end if

        ! dimensionMask
        if( present(dimensionMask) ) then 
          locDimensionMask = dimensionMask
        else
          locDimensionMask = (/1,1,1/)
        end if


        ! Assign dim properties
        !this%nDim = sum(dimensionMask, mask=(locDimensionMask.eq.1))
        this%nDim = sum((/1,1,1/), mask=(binSize.gt.0d0))
        if ( this%nDim .le. 0 ) then 
          write(*,*)'Error while initializing Histogram, nDim .le. 0. Stop.'
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
          if ( binSize(nd) .le. 0d0 ) cycle
          !if ( locDimensionMask(nd) .eq. 0 ) cycle
          dcount = dcount + 1
          this%dimensions(dcount) = nd
        end do 


        ! domainOrigin
        if ( present( domainOrigin ) ) then 
          this%domainOrigin = domainOrigin
        else 
          this%domainOrigin = 0 
        end if


        ! Initialize variables
        this%nBins     = nBins
        this%binSize   = binSize
        !this%binVolume = product( binSize, mask=(locDimensionMask.eq.1) ) 
        !this%binDistance = ( this%binVolume )**(1d0/sum(locDimensionMask))
        this%binVolume = product( binSize, mask=(binSize.gt.0d0) ) 
        this%binDistance = ( this%binVolume )**(1d0/this%nDim)


        ! Allocate and initialize histogram counts
        allocate( this%counts( nBins(1), nBins(2), nBins(3) ) )
        allocate( this%ncounts( nBins(1), nBins(2), nBins(3) ) )
        allocate( this%avgmbin( nBins(1), nBins(2), nBins(3) ) )
        this%counts = 0


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

        this%nBins       = 0
        this%nActiveBins = 0 
        this%binSize     = 0d0
        this%binVolume   = 0d0

        if( allocated(this%counts) ) deallocate( this%counts )
        if( allocated(this%ncounts) ) deallocate( this%ncounts )
        if( allocated(this%avgmbin) ) deallocate( this%avgmbin )
        if( allocated(this%activeBinIds) ) deallocate( this%activeBinIds ) 
        if( allocated(this%boundingBoxBinIds) ) deallocate( this%boundingBoxBinIds ) 
        if( allocated( this%dimensions) ) deallocate( this%dimensions) 

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
        doubleprecision, dimension(:,:), intent(in) :: dataPoints
        logical, intent(in), optional      :: exact 
        integer, dimension(2)              :: nPointsShape
        integer                            :: np, nd, did
        integer, dimension(3)              :: gridIndexes
        logical :: inside
        integer :: exactIndex 
        !------------------------------------------------------------------------------

        exactIndex = 1
        if( present(exact) ) then 
          if ( exact ) exactIndex = 0
        end if

        ! Reset counts
        this%counts  = 0
        nPointsShape = shape(dataPoints)

        ! This could be done with OpenMP (?) 
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
              gridIndexes(did) = floor(( dataPoints(np,did) - this%domainOrigin(did))/this%binSize(did)) + exactIndex
              if( (gridIndexes(did) .gt. this%nBins(did)) .or.&
                  (gridIndexes(did) .le. 0) ) then
                inside = .false. 
                exit
              end if 
            end do
            if ( .not. inside ) cycle

            ! Increase counter
            this%counts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) = &
                this%counts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) + 1d0

        end do 
       
    end subroutine prComputeCounts


    subroutine prComputeCountsWeighted( this, dataPoints, weights, exact )
        !------------------------------------------------------------------------------
        ! 
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(HistogramType) :: this
        doubleprecision, dimension(:,:), intent(in) :: dataPoints
        doubleprecision, dimension(:), intent(in)   :: weights
        logical, intent(in), optional      :: exact 
        integer, dimension(2)              :: nPointsShape
        integer                            :: np, nd, did
        integer, dimension(3)              :: gridIndexes
        logical :: inside
        integer :: exactIndex 
        !------------------------------------------------------------------------------

        exactIndex = 1
        if( present(exact) ) then 
          if ( exact ) exactIndex = 0
        end if

        ! Reset counts
        this%counts  = 0
        this%ncounts  = 0
        this%avgmbin  = 0
        nPointsShape = shape(dataPoints)

        ! This could be done with OpenMP (?) 
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
              gridIndexes(did) = floor(( dataPoints(np,did) - this%domainOrigin(did))/this%binSize(did)) + exactIndex
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
            this%ncounts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) = &
                this%ncounts( gridIndexes(1), gridIndexes(2), gridIndexes(3) ) + 1

        end do 
      
        where (this%ncounts.ne.0)  
          this%avgmbin = this%counts/this%ncounts
        end where
        where (this%avgmbin.ne.0)  
          this%counts = this%counts/this%avgmbin
        end where

    end subroutine prComputeCountsWeighted



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

        this%nActiveBins = count( this%counts/=0d0 )
        !this%nActiveBins = count( this%counts/=0 )

        ! Reallocate activeBinIds to new size 
        if ( allocated( this%activeBinIds ) )  deallocate( this%activeBinIds )
        allocate( this%activeBinIds( 3, this%nActiveBins ) )
        
        ! Following column-major nesting
        ! This could be in parallel with OpenMP (?)
        do iz = 1, this%nBins(3)
            do iy = 1, this%nBins(2)
                do ix = 1, this%nBins(1)
                    if ( this%counts( ix, iy, iz ) .eq. 0d0 ) cycle
                    !if ( this%counts( ix, iy, iz ) .eq. 0 ) cycle
                    this%activeBinIds( :, icount ) = [ ix, iy, iz ]
                    icount = icount + 1
                end do
            end do
        end do


    end subroutine prComputeActiveBinIds


    
    subroutine prComputeBoundingBox( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class( HistogramType ) :: this
        integer, dimension(2)  :: xBounds
        integer, dimension(2)  :: yBounds
        integer, dimension(2)  :: zBounds
        integer                :: ix, iy, iz
        integer                :: icount = 1 
        !------------------------------------------------------------------------------

        ! If no active bins, compute them 
        if ( .not. allocated( this%activeBinIds ) ) then 
            call this%ComputeActiveBinIds()
        end if 

        ! Get bounding box boundaries
        xBounds(1) = minval( this%activeBinIds( 1, : ) )
        xBounds(2) = maxval( this%activeBinIds( 1, : ) )
        yBounds(1) = minval( this%activeBinIds( 2, : ) )
        yBounds(2) = maxval( this%activeBinIds( 2, : ) )
        zBounds(1) = minval( this%activeBinIds( 3, : ) )
        zBounds(2) = maxval( this%activeBinIds( 3, : ) )

        ! Compute number of bins
        this%nBBoxBins = ( xBounds(2) - xBounds(1) + 1 )*&
                         ( yBounds(2) - yBounds(1) + 1 )*&
                         ( zBounds(2) - zBounds(1) + 1 )


        ! Maybe something that verifies the number of bins 
        if ( allocated( this%boundingBoxBinIds ) ) deallocate( this%boundingBoxBinIds )
        allocate( this%boundingBoxBinIds( 3, this%nBBoxBins ) )


        do iz = 1, this%nBins(3)
            do iy = 1, this%nBins(2)
                do ix = 1, this%nBins(1)
                    if ( &
                        ! OUTSIDE
                        ( ix .lt. xBounds(1) ) .or. ( ix .gt. xBounds(2) ) .or. &
                        ( iy .lt. yBounds(1) ) .or. ( iy .gt. yBounds(2) ) .or. &
                        ( iz .lt. zBounds(1) ) .or. ( iz .gt. zBounds(2) )      &
                    ) cycle
                    this%boundingBoxBinIds( :, icount ) = [ ix, iy, iz ]
                    icount = icount + 1
                end do
            end do
        end do


        return
        

    end subroutine 




    subroutine prExport( this )
        !------------------------------------------------------------------------------
        ! 
        !
        !------------------------------------------------------------------------------
        ! Specifications 
        !------------------------------------------------------------------------------
        implicit none 
        class(HistogramType) :: this
        integer              :: ix, iy, iz
        integer              :: histogramUnit = 555
        character(len=200)   :: outputFileName
        !------------------------------------------------------------------------------

        ! Write the output file name

        write( unit=outputFileName, fmt='(a)' )'histogram_output_.hist'
        !write( unit=outputFileName, fmt='(a)' )'histogram_output_'//trim(adjustl(tempTimeId))//'.hist'
        open(  histogramUnit, file=outputFileName, status='replace' )

        ! Write data, only non-zero values
        do ix = 1, this%nBins(1)
            do iy = 1, this%nBins(2)
                do iz = 1, this%nBins(3)
                    if ( this%counts( ix, iy, iz ) .eq. 0 ) cycle
                    write(histogramUnit,"(I6,I6,I6,es18.9e3)") ix, iy, iz, this%counts( ix, iy, iz )
                    !write(histogramUnit,"(I6,I6,I6,I6)") ix, iy, iz, this%counts( ix, iy, iz )
                end do
            end do
        end do

        ! Finished
        close(histogramUnit)


    end subroutine prExport




end module HistogramModule

