module TransportModelDataModule
  use ModflowRectangularGridModule,only : ModflowRectangularGridType
  use ParticleTrackingOptionsModule,only : ParticleTrackingOptionsType
  implicit none
  !---------------------------------------------------------------------------------------------------------------

  ! Set default access status to private
  private

    type,public :: TransportModelDataType

      logical :: Initialized = .false.
      doubleprecision,dimension(:),allocatable :: AlphaLong
      doubleprecision,dimension(:),allocatable :: AlphaTrans
      doubleprecision :: DMol
      class(ModflowRectangularGridType),pointer :: Grid => null()

    contains

      procedure :: Initialize=>pr_Initialize
      procedure :: Reset=>pr_Reset
      procedure :: ReadData=>pr_ReadData

    end type


contains


    subroutine pr_Initialize( this, grid )
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(TransportModelDataType) :: this
    class(ModflowRectangularGridType),intent(inout),pointer :: grid
    integer :: cellCount, gridType
    !---------------------------------------------------------------------------------------------------------------

        this%Initialized = .false.
        
        ! Call Reset to make sure that all arrays are initially unallocated
        call this%Reset()
        
        ! Return if the grid cell count equals 0
        cellCount = grid%CellCount
        if(cellCount .le. 0) return
        
        ! Check budget reader and grid data for compatibility and allocate appropriate cell-by-cell flow arrays
        gridType = grid%GridType
        this%Grid => grid

        this%Initialized = .true.


    end subroutine pr_Initialize


    subroutine pr_ReadData(this, inUnit, inFile, outUnit, trackingOptions )
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    use utl7module,only : urdcom, upcase
    use UTL8MODULE,only : urword, ustop, u1dint, u1drel, u1ddbl, u8rdcom,         &
      u3dintmp, u3dintmpusg, u3ddblmp, u3ddblmpusg
    implicit none
    class(TransportModelDataType) :: this
    integer,intent(in) :: inUnit, outUnit
    character(len=200), intent(in)    :: inFile 
    type(ParticleTrackingOptionsType),intent(inout) :: trackingOptions
    integer,dimension(:),allocatable :: cellsPerLayer
    integer,dimension(:),allocatable :: layerTypes
    !class(ModflowRectangularGridType),intent(inout) :: grid
    class(ModflowRectangularGridType), pointer :: grid => null()
    character(len=200) :: line
    character(len=16) :: txt
    integer :: n, m, nn, length, iface, errorCode, layer
    integer :: icol, istart, istop
    integer :: iodispersion = 0
    doubleprecision :: r
    character(len=24),dimension(2) :: aname
    data aname(1) /'          BOUNDARY ARRAY'/
    data aname(2) /'            DISPERSIVITY'/
    integer :: tempAlphaUnit = 666
    character(len=200) :: tempAlphaFile
    !---------------------------------------------------------------------------------------------------------------

        
        ! Open dispersion unit 
        open( inUnit, file=inFile, status='old', access='sequential')

        ! Initialize pointer
        grid => null() 
        grid => this%Grid

        if(allocated(this%AlphaLong)) deallocate(this%AlphaLong)
        if(allocated(this%AlphaTrans)) deallocate(this%AlphaTrans)
        allocate(this%AlphaTrans(grid%CellCount))
        allocate(this%AlphaLong(grid%CellCount))
        allocate(cellsPerLayer(grid%LayerCount))
        do n = 1, grid%LayerCount
            cellsPerLayer(n) = grid%GetLayerCellCount(n)
        end do
    
        ! Write header to the listing file
        write(outUnit, *)
        write(outUnit, '(1x,a)') 'MODPATH RWPT dispersion data file'
        write(outUnit, '(1x,a)') '----------------------------'


        ! Read dispersivity variable, eventually file
        ! These methods follor OPEN/CLOSE, CONSTANT input format
        ! and variables are expected to be defined for each layer

        ! ALPHALONG
        if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, this%AlphaLong, ANAME(2))                      
        else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              this%AlphaLong, aname(2), cellsPerLayer)
        else
            write(outUnit,*) 'Invalid grid type specified when reading ALPHALONG array data.'
            write(outUnit,*) 'Stopping.'
            call ustop(' ')          
        end if
        

        ! ALPHATRANS
        if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, this%AlphaTrans, ANAME(2))                      
        else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              this%AlphaTrans, aname(2), cellsPerLayer)
        else
              write(outUnit,*) 'Invalid grid type specified when reading ALPHATRANS array data.'
              write(outUnit,*) 'Stopping.'
              call ustop(' ')          
        end if


        ! Dmol
        read( inUnit, * ) line
        icol = 1
        call urword(line,icol,istart,istop,3,n,r,0,0)
        this%DMol = r
        trackingOptions%Dmol = r

        ! Time Step kind 
        read(inUnit, '(a)') line
        icol = 1
        call urword(line,icol,istart,istop,1,n,r,0,0)
        line = line(istart:istop)
        ! Give a Courant
        if ( line .eq. 'CONSTANT_CU' ) then
            trackingOptions%timeStepKind = 1
            read( inUnit, * ) line
            icol = 1
            call urword(line,icol,istart,istop,3,n,r,0,0)
            trackingOptions%timeStepParameters(1) = r
        ! Give a Peclet
        else if ( line .eq. 'CONSTANT_PE' ) then
            trackingOptions%timeStepKind = 2
            read( inUnit, * ) line
            icol = 1
            call urword(line,icol,istart,istop,3,n,r,0,0)
            trackingOptions%timeStepParameters(2) = r
        ! Minimum between Courant and Peclet 
        else if ( line .eq. 'MIN_CU_PE' ) then
            trackingOptions%timeStepKind = 3
            read( inUnit, * ) line
            icol = 1
            call urword(line,icol,istart,istop,3,n,r,0,0)
            trackingOptions%timeStepParameters(1) = r
            read( inUnit, * ) line
            icol = 1
            call urword(line,icol,istart,istop,3,n,r,0,0)
            trackingOptions%timeStepParameters(2) = r
        else
            call ustop('RWPT: Invalid options for time step selection. Stop.')
        end if


        ! Advection Integration Kind
        read(inUnit, '(a)', iostat=iodispersion) line
        if ( iodispersion .lt. 0 ) then 
            ! end of file
            trackingOptions%advectionKind = 1
            write(outUnit,'(A)') 'RWPT: Advection integration not specified. Default to EXPONENTIAL.'
        else
            icol = 1
            call urword(line,icol,istart,istop,1,n,r,0,0)
            line = line(istart:istop)
            select case(line)
                case('EXPONENTIAL')
                    trackingOptions%advectionKind = 1
                    write(outUnit,'(A)') 'RWPT: Advection integration is EXPONENTIAL.'
                case('EULERIAN')
                    trackingOptions%advectionKind = 2
                    write(outUnit,'(A)') 'RWPT: Advection integration is EULERIAN.'
                case default
                    trackingOptions%advectionKind = 1
                    write(outUnit,'(A)') 'RWPT: Advection integration not specified. Default to EXPONENTIAL.'
            end select
        end if


        ! Domain dimensions option
        read(inUnit, '(a)', iostat=iodispersion) line
        if ( iodispersion .lt. 0 ) then 
            ! end of file
            trackingOptions%twoDimensions = .false.
            write(outUnit,'(A)') 'RWPT: Number of dimensions not specified. Defaults to 3D.'
        else
            icol = 1
            call urword(line,icol,istart,istop,1,n,r,0,0)
            line = line(istart:istop)
            select case(line)
                case('2D')
                    trackingOptions%twoDimensions = .true.
                    write(outUnit,'(A)') 'RWPT: Selected 2D domain solver.'
                case('3D')
                    trackingOptions%twoDimensions = .false.
                    write(outUnit,'(A)') 'RWPT: Selected 3D domain solver.'
                case default
                    trackingOptions%twoDimensions = .false.
                    write(outUnit,'(A)') 'RWPT: Invalid option for domain solver, defaults to 3D.'
            end select
        end if 


        ! Close dispersion data file
        close( inUnit )


    end subroutine pr_ReadData


    subroutine pr_Reset(this)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(TransportModelDataType) :: this
    !---------------------------------------------------------------------------------------------------------------
       
        this%Grid => null()
        this%DMol = 0
        if(allocated(this%AlphaLong)) deallocate( this%AlphaLong )
        if(allocated(this%AlphaTrans)) deallocate( this%AlphaTrans )

    end subroutine pr_Reset



end module TransportModelDataModule
