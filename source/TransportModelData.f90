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

      !logical :: SteadyState = .true.
      !doubleprecision,dimension(:),pointer :: AlphaLong
      !doubleprecision,dimension(:),pointer :: AlphaTrans
      !doubleprecision,dimension(:),allocatable :: Porosity

      class(ModflowRectangularGridType),pointer :: Grid => null()

    contains

      procedure :: Initialize=>pr_Initialize
      procedure :: Reset=>pr_Reset
      procedure :: ReadData=>pr_ReadData
      !procedure :: SetDispersivities=>pr_SetDispersivities

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
        !select case (gridType)
        !    case (1)
        !        if((budgetReader%GetBudgetType() .ne. 1)) return
        !        if((headReader%GridStyle .ne. 1) .or. (headReader%CellCount .ne. cellCount)) return
        !        flowArraySize = budgetReader%GetTransportArraySize()
        !        if(flowArraySize .ne. cellCount) return
        !        allocate(this%TransportsRightFace(flowArraySize))
        !        allocate(this%TransportsFrontFace(flowArraySize))
        !        allocate(this%TransportsLowerFace(flowArraySize))
        !        allocate(this%TransportsJA(0))
        !    case (2)
        !        if((budgetReader%GetBudgetType() .ne. 2)) return
        !        if((headReader%GridStyle .ne. 2) .or. (headReader%CellCount .ne. cellCount)) return
        !        flowArraySize = budgetReader%GetTransportArraySize()
        !        if(flowArraySize .ne. grid%JaCount) return
        !        allocate(this%TransportsJA(flowArraySize))
        !        allocate(this%TransportsRightFace(0))
        !        allocate(this%TransportsFrontFace(0))
        !        allocate(this%TransportsLowerFace(0))
        !    case (3, 4)
        !        if((budgetReader%GetBudgetType() .ne. 2)) return
        !        if((headReader%GridStyle .ne. 1) .or. (headReader%CellCount .ne. cellCount)) return
        !        flowArraySize = budgetReader%GetTransportArraySize()
        !        if(flowArraySize .ne. grid%JaCount) return
        !        allocate(this%TransportsJA(flowArraySize))
        !        allocate(this%TransportsRightFace(0))
        !        allocate(this%TransportsFrontFace(0))
        !        allocate(this%TransportsLowerFace(0))          
        !    !case (4)
        !    !    ! Not implemented
        !    !    return
        !    case default
        !        return
        !end select
        !
        !! Set pointers to budgetReader and grid. Assign tracking options.
        !this%HeadReader => headReader
        !this%BudgetReader => budgetReader
        this%Grid => grid
        !this%HNoTransport = hNoTransport
        !this%HDry = hDry
        !
        !! Allocate the rest of the arrays
        !allocate(this%IBoundTS(cellCount))
        !allocate(this%Heads(cellCount))
        !allocate(this%SourceTransports(cellCount))
        !allocate(this%SinkTransports(cellCount))
        !allocate(this%StorageTransports(cellCount))
        !allocate(this%SubFaceTransportsComputed(cellCount))
        !allocate(this%BoundaryTransports(cellCount * 6))
        !allocate(this%SubFaceTransports(cellCount * 4))
        !
        !! Allocate buffers for reading array and list data
        !allocate(this%ListItemBuffer(this%BudgetReader%GetMaximumListItemCount()))
        !allocate(this%ArrayBufferDbl(this%BudgetReader%GetMaximumArrayItemCount()))  
        !allocate(this%ArrayBufferInt(this%BudgetReader%GetMaximumArrayItemCount()))
        
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


        ! Read longitudinal dispersivity string file
        read( inUnit, '(a)' ) line
        icol = 1
        call urword(line,icol,istart,istop,0,n,r,0,0)
        tempAlphaFile = line(istart:istop)
        open( tempAlphaUnit, file=tempAlphaFile, status='old', access='sequential' )
        write(outUnit, '(1x,a)') 'AlphaLong file is ', tempAlphaFile
        write(outUnit, *)

        ! ALPHALONG
        if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(tempAlphaUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, this%AlphaLong, ANAME(2))                      
        else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(tempAlphaUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              this%AlphaLong, aname(2), cellsPerLayer)
        else
              write(outUnit,*) 'Invalid grid type specified when reading ALPHALONG array data.'
              write(outUnit,*) 'Stopping.'
              call ustop(' ')          
        end if
        
        close( tempAlphaUnit )

        ! Read transverse dispersivity string file
        read( inUnit, '(a)' ) line
        icol = 1
        call urword(line,icol,istart,istop,0,n,r,0,0)
        tempAlphaFile = line(istart:istop)
        open( tempAlphaUnit, file=tempAlphaFile, status='old', access='sequential' )
        write(outUnit, '(1x,a)') 'AlphaTrans file is ', tempAlphaFile
        write(outUnit, *)

        ! ALPHATRANS
        if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3ddblmp(tempAlphaUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, this%AlphaTrans, ANAME(2))                      
        else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3ddblmpusg(tempAlphaUnit, outUnit, grid%CellCount, grid%LayerCount,  &
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
       
        !this%ReferenceTime = 0.0d0
        !this%StoppingTime = 0.0d0
        !this%CurrentStressPeriod = 0
        !this%CurrentTimeStep = 0
        !this%HeadReader => null()
        !this%BudgetReader => null()
        this%Grid => null()
        !this%AlphaLong => null()
        !this%AlphaTrans => null()
        this%DMol = 0
        if(allocated(this%AlphaLong)) deallocate( this%AlphaLong )
        if(allocated(this%AlphaTrans)) deallocate( this%AlphaTrans )

        !if(allocated(this%IBoundTS)) deallocate(this%IBoundTS)
        !if(allocated(this%ArrayBufferInt)) deallocate(this%ArrayBufferInt)
        !if(allocated(this%Heads)) deallocate(this%Heads)
        !if(allocated(this%TransportsJA)) deallocate(this%TransportsJA)
        !if(allocated(this%TransportsRightFace)) deallocate(this%TransportsRightFace)
        !if(allocated(this%TransportsFrontFace)) deallocate(this%TransportsFrontFace)
        !if(allocated(this%TransportsLowerFace)) deallocate(this%TransportsLowerFace)
        !if(allocated(this%SourceTransports)) deallocate(this%SourceTransports)
        !if(allocated(this%SinkTransports)) deallocate(this%SinkTransports)
        !if(allocated(this%StorageTransports)) deallocate(this%StorageTransports)
        !if(allocated(this%BoundaryTransports)) deallocate(this%BoundaryTransports)
        !if(allocated(this%SubFaceTransports)) deallocate(this%SubFaceTransports)
        !if(allocated(this%ArrayBufferDbl)) deallocate(this%ArrayBufferDbl)
        !if(allocated(this%ListItemBuffer)) deallocate(this%ListItemBuffer)
        !if(allocated(this%SubFaceTransportsComputed)) deallocate(this%SubFaceTransportsComputed)
        !this%IBound => null()
        !this%Porosity => null()
        !this%Retardation => null()
        !this%Zones => null()


    end subroutine pr_Reset



    !subroutine pr_SetDispersivities(this, alphaLong, alphaTrans, arraySize)
    !!***************************************************************************************************************
    !!
    !!***************************************************************************************************************
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !implicit none
    !class(TransportModelDataType) :: this
    !integer,intent(in) :: arraySize
    !doubleprecision,dimension(arraySize),intent(in),target :: alphaLong
    !doubleprecision,dimension(arraySize),intent(in),target :: alphaTrans
    !!---------------------------------------------------------------------------------------------------------------


    !    if(arraySize .ne. this%Grid%CellCount) then
    !        write(*,*) "TransportModelDataType: Dispersivities array size does not match the cell count for the grid. stop"
    !        stop
    !    end if
    !    
    !    this%AlphaLong => alphaLong
    !    this%AlphaTrans => alphaTrans


    !end subroutine pr_SetDispersivities


    !! COULD BE REUSED FOR TRANSPORT BOUNDARY CONDITIONS
    !subroutine pr_SetIbound(this, ibound, arraySize)
    !!***************************************************************************************************************
    !!
    !!***************************************************************************************************************
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !implicit none
    !class(TransportModelDataType) :: this
    !integer,intent(in) :: arraySize
    !integer :: n
    !integer,dimension(arraySize),intent(in),target :: ibound
    !!---------------------------------------------------------------------------------------------------------------
    !  

    !    if(arraySize .ne. this%Grid%CellCount) then
    !        write(*,*) "TransportModelDataType: The IBound array size does not match the cell count for the grid. stop"
    !        stop
    !    end if
    !    
    !    this%IBound => ibound
    !    ! Initialize the IBoundTS array to the same values as IBound whenever the IBound array is set.
    !    ! The values of IBoundTS will be updated for dry cells every time that data for a time step is loaded.
    !    do n = 1, arraySize
    !        this%IBoundTS(n) = this%IBound(n)
    !    end do

    !end subroutine pr_SetIbound


end module TransportModelDataModule


! THRASH 





    !! READ MPBAS DATA
    !! READ AND PRINT COMMENTS.
    !call u8rdcom(inUnit,outUnit,line,errorCode)
 
    !! No flow and dry cell head flags
    !if(grid%GridType .gt. 2) then
    !    this%HNoFlow = 1.0E+30
    !    this%HDry = -1.0E+30
    !    read(line, *) this%DefaultIfaceCount
    !    write(outUnit,'(1X,A,1PG12.5,A)') 'Aquifer head is set to ', this%HNoFlow,    &
    !      ' at cells with IDOMAIN = 0.'
    !    write(outUnit,'(1X,A,1PG12.5,A)') 'Aquifer head is set to ', this%Hdry,       &
    !      ' at all dry cells when Newton-Raphson rewetting is not active.'
    !    write(outUnit, *)
    !else
    !    read(line,*) this%HNoFlow, this%HDry
    !    read(inUnit, *) this%DefaultIfaceCount
    !    write(outUnit,'(1X,A,1PG12.5,A)') 'Aquifer head is set to ', this%HNoFlow,    &
    !      ' at all no-flow cells (IBOUND=0).'
    !    write(outUnit,'(1X,A,1PG12.5,A)') 'Aquifer head is set to ', this%Hdry,       &
    !      ' at all dry cells.'
    !    write(outUnit, *)
    !end if
    !
    !    
    !    ! READ NUMBER OF STRESS PACKAGES THAT DO NOT HAVE AUXILIARY IFACE SUPPORT
    !    allocate (this%DefaultIfaceValues(this%DefaultIfaceCount))
    !    allocate (this%DefaultIfaceLabels(this%DefaultIfaceCount))
    !    
    !    if(this%DefaultIfaceCount .gt. 0) then
    !      ! READ BUDGET LABELS AND DEFAULT IFACE
    !      do n = 1, this%DefaultIfaceCount
    !        read(inUnit,'(A)') line
    !        call utrimall(line)
    !        length = len_trim(line)
    !        if(length .eq. 0) THEN
    !          call ustop('Stress package labels cannot be blank strings.')
    !        else
    !          if(length .gt. 16) length = 16
    !          this%DefaultIfaceLabels(n) = line
    !          call upcase(this%DefaultIfaceLabels(n))
    !        end if

    !        ! READ DEFAULT IFACE          
    !        read(inUnit,*) this%DefaultIfaceValues(n)
    !        
    !      end do

    !      ! Check for duplicate labels        
    !      do n = 1, this%DefaultIfaceCount
    !        txt = this%DefaultIfaceLabels(n)
    !        do m = 1, this%DefaultIfaceCount
    !          if(m .NE. n) then
    !            if(txt .eq. this%DefaultIfaceLabels(m)) THEN
    !             call ustop('Duplicate stress package labels are not allowed.')
    !            end if
    !          end if
    !        end do
    !      end do
    !      
    !      write(outUnit, *)
    !      write(outUnit,'(1X,A,A)') 'Default stress package boundary face options (IFACE):'
    !      if(this%DefaultIfaceCount .eq. 0) THEN
    !        write(outUnit,'(3X,A)') 'None were specified'
    !      else
    !        do n = 1, this%DefaultIfaceCount
    !          iface = this%DefaultIfaceValues(n)
    !          if(iface .eq. 0) THEN
    !            write(outUnit,'(3X,A,A)') trim(this%DefaultIfaceLabels(n)),       &
    !              ' will be treated as internal stresses (IFACE = 0)'
    !          else if((iface .lt. 0) .or. (iface .gt. 6)) THEN
    !            call ustop(' IFACE must have a value between 0 and 6.')
    !          else
    !            write(outUnit,'(3X,A,A,I2)') trim(this%DefaultIfaceLabels(n)),    &
    !              ' will be assigned to face ',IFACE
    !          end if
    !        end do
    !      end if
    !      
    !    end if
    !
    !    ! LAYTYP
    !    ! Read LAYTYP only for MODFLOW-2005 and MODFLOW-USG datasets.
    !    ! Loop through layers and set the cell type for all cells in a 
    !    ! layer equal to the LAYTYP value for the layer.
    !    ! LAYTYP is not read for MODFLOW-6 datasets. Instead, cell types 
    !    ! for all cells are read directly for from the binary grid file.
    !    if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 2)) then
    !        write(outUnit, *)
    !        write(outUnit,'(1X,A)') 'Layer type (LAYTYP)'
    !        allocate(layerTypes(grid%LayerCount))
    !        read(inUnit,*) (layerTypes(layer), layer = 1, grid%LayerCount)
    !        write(outUnit,'(1X,80I2)') (layerTypes(layer), layer = 1, grid%LayerCount)
    !        n = 0
    !        do layer = 1, grid%LayerCount
    !          if(layerTypes(layer) .lt. 0) layerTypes(layer) = 1
    !          do m = 1, cellsPerLayer(layer)
    !              n = n + 1
    !              grid%CellType(n) = layerTypes(layer)
    !          end do
    !        end do
    !    end if
    !    
    !    ! IBOUND
    !    if(grid%GridType .eq. 1) then
    !          call u3dintmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
    !            grid%ColumnCount, grid%CellCount, this%IBound, ANAME(1))                      
    !    else if(grid%GridType .eq. 2) then
    !        call u3dintmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount, this%IBound, &
    !          aname(1), cellsPerLayer)
    !    else if(grid%GridType .ge. 3) then
    !        do n = 1, grid%CellCount
    !            if(grid%GetIDomain(n) .lt. 1) then
    !                this%IBound(n) = 0
    !            else
    !                this%IBound(n) = 1
    !            end if
    !        end do
    !    else
    !          write(outUnit,*) 'Invalid grid type specified when reading IBOUND array data.'
    !          write(outUnit,*) 'Stopping.'
    !          call ustop(' ')          
    !    end if









      !doubleprecision :: ReferenceTime = 0d0
      !doubleprecision :: StoppingTime = 0d0
      !doubleprecision :: HDry = 0d0
      !doubleprecision :: HNoTransport = 0d0
      !doubleprecision,dimension(:),pointer :: Porosity
      !doubleprecision,dimension(:),pointer :: Retardation
      !type(HeadReadertype)  , pointer :: HeadReader   => null()
      !type(BudgetReaderType), pointer :: BudgetReader => null()
      !class(ModflowRectangularGridType),pointer :: Grid => null()
      !type(BudgetListItemType),allocatable,dimension(:) :: ListItemBuffer
      !logical,allocatable,dimension(:), private :: SubFaceTransportsComputed
      !! Private variables
      !integer,private :: CurrentStressPeriod = 0
      !integer,private :: CurrentTimeStep = 0
      !integer :: DefaultIfaceCount
      !character(len=16),dimension(20) :: DefaultIfaceLabels
      !integer,dimension(20) :: DefaultIfaceValues
      !integer,allocatable,dimension(:) :: IBoundTS
      !integer,allocatable,dimension(:) :: ArrayBufferInt
      !doubleprecision,allocatable,dimension(:) :: Heads
      !doubleprecision,allocatable,dimension(:) :: TransportsJA
      !doubleprecision,allocatable,dimension(:) :: TransportsRightFace
      !doubleprecision,allocatable,dimension(:) :: TransportsFrontFace
      !doubleprecision,allocatable,dimension(:) :: TransportsLowerFace
      !doubleprecision,allocatable,dimension(:) :: SourceTransports
      !doubleprecision,allocatable,dimension(:) :: SinkTransports
      !doubleprecision,allocatable,dimension(:) :: StorageTransports
      !doubleprecision,allocatable,dimension(:) :: BoundaryTransports
      !doubleprecision,allocatable,dimension(:) :: SubFaceTransports
      !doubleprecision,allocatable,dimension(:) :: ArrayBufferDbl



    !subroutine pr_SetZones(this, zones, arraySize)
    !!***************************************************************************************************************
    !!
    !!***************************************************************************************************************
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !implicit none
    !class(TransportModelDataType) :: this
    !integer,intent(in) :: arraySize
    !integer,dimension(arraySize),intent(in),target :: zones
    !!---------------------------------------------------------------------------------------------------------------
    !  

    !    if(arraySize .ne. this%Grid%CellCount) then
    !        write(*,*) "TransportModelDataType: The Zones array size does not match the cell count for the grid. stop"
    !        stop
    !    end if
    !    
    !    this%Zones => zones
    !

    !end subroutine pr_SetZones



    !subroutine pr_SetLayerTypes(this, layerTypes, arraySize)
    !!***************************************************************************************************************
    !!
    !!***************************************************************************************************************
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !  implicit none
    !  class(TransportModelDataType) :: this
    !  integer,intent(in) :: arraySize
    !  integer,dimension(arraySize),intent(in),target :: layerTypes
    !!---------------------------------------------------------------------------------------------------------------
    !
    !  
    !  if(arraySize .ne. this%Grid%LayerCount) then
    !      write(*,*) "TransportModelDataType: The LayerTypes array size does not match the layer count for the grid. stop"
    !      stop
    !  end if
    !  
    !  this%LayerTypes => layerTypes
    !
    !
    !end subroutine pr_SetLayerTypes


    !subroutine pr_LoadTimeStep(this, stressPeriod, timeStep)
    !!***************************************************************************************************************
    !!
    !!***************************************************************************************************************
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !implicit none
    !class(TransportModelDataType) :: this
    !integer,intent(in) :: stressPeriod, timeStep
    !integer :: firstRecord, lastRecord, n, m, firstNonBlank, lastNonBlank, &
    !  trimmedLength
    !integer :: spaceAssigned, status,cellCount, iface, index,              &
    !  boundaryTransportsOffset, listItemBufferSize, cellNumber, layer
    !type(BudgetRecordHeaderType) :: header
    !character(len=16) :: textLabel
    !doubleprecision :: top 
    !real :: HDryTol, HDryDiff
    !!---------------------------------------------------------------------------------------------------------------
    !  
    !    call this%ClearTimeStepBudgetData()
    !    call this%BudgetReader%GetRecordHeaderRange(stressPeriod, timeStep, firstRecord, lastRecord)
    !    if(firstRecord .eq. 0) return
    !
    !    cellCount = this%Grid%CellCount
    !    listItemBufferSize = size(this%ListItemBuffer)
    !    
    !    ! Set steady state = true, then change it if the budget file contains storage
    !    this%SteadyState = .true.
    !    
    !    ! Load heads for this time step
    !    call this%HeadReader%FillTimeStepHeadBuffer(stressPeriod, timeStep, &
    !      this%Heads, cellCount, spaceAssigned)
    !    
    !    ! Fill IBoundTS array and set the SaturatedTop array for the Grid.
    !    ! The saturated top is set equal to the top for confined cells and water table cells 
    !    ! where the head is above the top or below the bottom.
    !    HDryTol = abs(epsilon(HDryTol)*sngl(this%HDry))
    !    if(this%Grid%GridType .gt. 2) then
    !        do n = 1, cellCount
    !            this%Grid%SaturatedTop(n) = this%Grid%Top(n)
    !            this%StorageTransports(n) = 0.0
    !            this%IBoundTS(n) = this%IBound(n)
    !            layer = this%Grid%GetLayer(n)
    !            if(this%Grid%CellType(n) .eq. 1) then
    !                HDryDiff = sngl(this%Heads(n)) - sngl(this%HDry)
    !                if(abs(HDryDiff) .lt. HDryTol) then
    !                    this%IBoundTS(n) = 0
    !                    if(this%Heads(n) .lt. this%Grid%Bottom(n)) then
    !                        this%IBoundTS(n) = 0
    !                        this%Grid%SaturatedTop(n) = this%Grid%Bottom(n)
    !                    end if
    !                end if
    !                if(this%IBoundTS(n) .ne. 0) then
    !                    if((this%Heads(n) .le. this%Grid%Top(n)) .and. &
    !                      (this%Heads(n) .ge. this%Grid%Bottom(n))) then
    !                        this%Grid%SaturatedTop(n) = this%Heads(n)
    !                    end if
    !                end if
    !            end if
    !        end do
    !        
    !    else
    !        do n = 1, cellCount
    !            this%Grid%SaturatedTop(n) = this%Grid%Top(n)
    !            this%StorageTransports(n) = 0.0
    !            this%IBoundTS(n) = this%IBound(n)
    !            layer = this%Grid%GetLayer(n)
    !            if(this%Grid%CellType(n) .eq. 1) then
    !                HDryDiff = sngl(this%Heads(n)) - sngl(this%HDry)
    !                if((abs(HDryDiff) .lt. HDryTol) .or. (this%Heads(n) .gt. 1.0d+6)) then
    !                    this%IBoundTS(n) = 0
    !                end if
    !                if(this%IBoundTS(n) .ne. 0) then
    !                    if((this%Heads(n) .le. this%Grid%Top(n)) .and. &
    !                      (this%Heads(n) .ge. this%Grid%Bottom(n))) then
    !                        this%Grid%SaturatedTop(n) = this%Heads(n)
    !                    end if
    !                end if
    !            end if
    !        end do
    !    end if
    !    
    !    ! Loop through record headers
    !    do n = firstRecord, lastRecord
    !         header = this%BudgetReader%GetRecordHeader(n)
    !         textLabel = header%TextLabel
    !         call TrimAll(textLabel, firstNonBlank, lastNonBlank, trimmedLength)
    !         
    !         select case(textLabel(firstNonBlank:lastNonBlank))
    !         case('CONSTANT HEAD', 'CHD')
    !              ! Read constant head flows into the sinkTransports and sourceTransports arrays.
    !              ! For a standard budget file, Method = 0. For a compact budget file,
    !              ! Method = 2.
    !              if(header%Method .eq. 0) then
    !                  call this%BudgetReader%FillRecordDataBuffer(header,             &
    !                    this%ArrayBufferDbl, cellCount, spaceAssigned, status)
    !                  if(cellCount .eq. spaceAssigned) then
    !                      do m = 1, spaceAssigned
    !                          if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
    !                              this%SourceTransports(m) = this%SourceTransports(m) +         &
    !                                this%ArrayBufferDbl(m)
    !                          else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
    !                              this%SinkTransports(m) = this%SinkTransports(m) +             &
    !                                this%ArrayBufferDbl(m)
    !                          end if
    !                      end do
    !                  end if
    !              else if(header%Method .eq. 2) then
    !                  call this%BudgetReader%FillRecordDataBuffer(header,             &
    !                    this%ListItemBuffer, listItemBufferSize, spaceAssigned, status)
    !                  if(spaceAssigned .gt. 0) then
    !                      do m = 1, spaceAssigned
    !                          cellNumber = this%ListItemBuffer(m)%CellNumber
    !                          if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
    !                              this%SourceTransports(cellNumber) =                      &
    !                                this%SourceTransports(cellNumber) + this%ListItemBuffer(m)%BudgetValue
    !                          else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
    !                              this%SinkTransports(cellNumber) =                        &
    !                                this%SinkTransports(cellNumber) + this%ListItemBuffer(m)%BudgetValue
    !                          end if
    !                      end do
    !                  end if
    !              else if((header%Method .eq. 5) .or. (header%Method .eq. 6)) then
    !                  call this%BudgetReader%FillRecordDataBuffer(header,             &
    !                    this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
    !                    status)
    !                  if(spaceAssigned .gt. 0) then
    !                      do m = 1, spaceAssigned
    !                          call this%CheckForDefaultIface(header%TextLabel, iface)
    !                          index = header%FindAuxiliaryNameIndex('IFACE')
    !                          if(index .gt. 0) then
    !                              iface = int(this%ListItemBuffer(m)%AuxiliaryValues(index))
    !                          end if
    !                          
    !                          cellNumber = this%ListItemBuffer(m)%CellNumber
    !                          if(iface .gt. 0) then
    !                              boundaryTransportsOffset = 6 * (cellNumber - 1)
    !                              this%BoundaryTransports(boundaryTransportsOffset + iface) =   &
    !                                this%BoundaryTransports(boundaryTransportsOffset + iface) + &
    !                                this%ListItemBuffer(m)%BudgetValue
    !                          else
    !                              if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
    !                                  this%SourceTransports(cellNumber) =                  &
    !                                    this%SourceTransports(cellNumber) +                &
    !                                    this%ListItemBuffer(m)%BudgetValue
    !                              else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
    !                                  this%SinkTransports(cellNumber) =                    &
    !                                    this%SinkTransports(cellNumber) +                  &
    !                                    this%ListItemBuffer(m)%BudgetValue
    !                              end if
    !                          end if
    !                      end do
    !                  end if
    !              end if
    !             
    !         case('STORAGE', 'STO-SS', 'STO-SY')
    !              ! Read storage for all cells into the StorageTransports array.
    !              ! Method should always be 0 or 1, but check anyway to be sure.
    !              if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
    !                  if(header%ArrayItemCount .eq. cellCount) then
    !                      call this%BudgetReader%FillRecordDataBuffer(header,         &
    !                        this%ArrayBufferDbl, cellCount, spaceAssigned, status)
    !                      if(cellCount .eq. spaceAssigned) then
    !                          do m = 1, spaceAssigned
    !                              this%StorageTransports(m) = this%StorageTransports(m) + this%ArrayBufferDbl(m)
    !                              if(this%StorageTransports(m) .ne. 0.0) this%SteadyState = .false.
    !                          end do
    !                      end if
    !                  end if
    !              end if
    !             
    !         case('FLOW JA FACE', 'FLOW-JA-FACE')
    !              ! Read connected face flows into the TransportsJA array for unstructured grids.
    !              if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
    !                  ! Method should always be 0 or 1 for flow between grid cells. 
    !                  if(header%ArrayItemCount .eq. this%BudgetReader%GetTransportArraySize()) then
    !                      call this%BudgetReader%FillRecordDataBuffer(header,         &
    !                        this%TransportsJA, header%ArrayItemCount, spaceAssigned,       &
    !                        status)
    !                  end if
    !              else if(header%Method .eq. 6) then
    !                  ! Method code 6 indicates flow to or from cells in the current model grid
    !                  ! and another connected model grid in a multi-model MODFLOW-6 simulation. 
    !                  ! Treat flows to or from connected model grids as distributed source/sink flows 
    !                  ! for the current grid.
    !                  call this%BudgetReader%FillRecordDataBuffer(header,             &
    !                    this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
    !                    status)
    !                  if(spaceAssigned .gt. 0) then
    !                      do m = 1, spaceAssigned
    !                          cellNumber = this%ListItemBuffer(m)%CellNumber
    !                          if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
    !                              this%SourceTransports(cellNumber) =                  &
    !                                  this%SourceTransports(cellNumber) +                &
    !                                  this%ListItemBuffer(m)%BudgetValue
    !                          else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
    !                              this%SinkTransports(cellNumber) =                    &
    !                                  this%SinkTransports(cellNumber) +                  &
    !                                  this%ListItemBuffer(m)%BudgetValue
    !                          end if
    !                      end do
    !                  end if
    !              end if
    !             
    !         case('FLOW RIGHT FACE')
    !              ! Read flows across the right face for structured grids.
    !              ! Method should always be 0 or 1, but check anyway to be sure.
    !              if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
    !                  if(header%ArrayItemCount .eq. this%BudgetReader%GetTransportArraySize()) then
    !                      call this%BudgetReader%FillRecordDataBuffer(header,         &
    !                        this%TransportsRightFace, header%ArrayItemCount, spaceAssigned,&
    !                        status)
    !                  end if
    !              end if
    !             
    !         case('FLOW FRONT FACE')
    !              ! Read flows across the front face for structured grids.
    !              ! Method should always be 0 or 1, but check anyway to be sure.
    !              if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
    !                  if(header%ArrayItemCount .eq. this%BudgetReader%GetTransportArraySize()) then
    !                      call this%BudgetReader%FillRecordDataBuffer(header,         &
    !                        this%TransportsFrontFace, header%ArrayItemCount, spaceAssigned,&
    !                        status)
    !                  end if
    !              end if
    !             
    !         case('FLOW LOWER FACE')
    !              ! Read flows across the lower face for structured grids.
    !              ! Method should always be 0 or 1, but check anyway to be sure.
    !              if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
    !                  if(header%ArrayItemCount .eq. this%BudgetReader%GetTransportArraySize()) then
    !                      call this%BudgetReader%FillRecordDataBuffer(header,         &
    !                        this%TransportsLowerFace, header%ArrayItemCount, spaceAssigned,&
    !                        status)
    !                  end if
    !              end if
    !         
    !          case default
    !              ! Now handle any other records in the budget file.
    !               if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
    !                  if(header%ArrayItemCount .eq. cellCount) then
    !                      call this%BudgetReader%FillRecordDataBuffer(header,         &
    !                        this%ArrayBufferDbl, cellCount, spaceAssigned, status)
    !                      if(cellCount .eq. spaceAssigned) then
    !                          call this%CheckForDefaultIface(header%TextLabel, iface)
    !                          if(iface .gt. 0) then
    !                              do m = 1, spaceAssigned
    !                                  boundaryTransportsOffset = 6 * (m - 1)
    !                                  this%BoundaryTransports(boundaryTransportsOffset + iface) =   &
    !                                    this%BoundaryTransports(boundaryTransportsOffset + iface) + &
    !                                    this%ArrayBufferDbl(m)
    !                              end do
    !                          else
    !                              do m = 1, spaceAssigned
    !                                  if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
    !                                      this%SourceTransports(m) = this%SourceTransports(m) +     &
    !                                        this%ArrayBufferDbl(m)
    !                                  else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
    !                                      this%SinkTransports(m) = this%SinkTransports(m) +         &
    !                                        this%ArrayBufferDbl(m)
    !                                  end if
    !                              end do
    !                          end if
    !                      end if
    !                  end if
    !               else if(header%Method .eq. 3) then
    !                  call this%BudgetReader%FillRecordDataBuffer(header,             &
    !                    this%ArrayBufferDbl, this%ArrayBufferInt,                     &
    !                    header%ArrayItemCount, spaceAssigned, status)
    !                  if(header%ArrayItemCount .eq. spaceAssigned) then
    !                      call this%CheckForDefaultIface(header%TextLabel, iface)
    !                      if(iface .gt. 0) then
    !                          do m = 1, spaceAssigned
    !                              cellNumber = this%ArrayBufferInt(m)
    !                              boundaryTransportsOffset = 6 * (cellNumber - 1)
    !                              this%BoundaryTransports(boundaryTransportsOffset + iface) =   &
    !                                this%BoundaryTransports(boundaryTransportsOffset + iface) + &
    !                                this%ArrayBufferDbl(m)
    !                          end do
    !                      else            
    !                          do m = 1, spaceAssigned
    !                              cellNumber = this%ArrayBufferInt(m)
    !                              if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
    !                                  this%SourceTransports(cellNumber) =                  &
    !                                    this%SourceTransports(cellNumber) +                &
    !                                    this%ArrayBufferDbl(m)
    !                              else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
    !                                  this%SinkTransports(cellNumber) =                    &
    !                                    this%SinkTransports(cellNumber) +                  &
    !                                    this%ArrayBufferDbl(m)
    !                              end if
    !                          end do
    !                      end if
    !                  end if
    !               else if(header%Method .eq. 4) then
    !                  call this%BudgetReader%FillRecordDataBuffer(header,             &
    !                    this%ArrayBufferDbl, header%ArrayItemCount, spaceAssigned,    &
    !                    status)
    !                  if(header%ArrayItemCount .eq. spaceAssigned) then
    !                      call this%CheckForDefaultIface(header%TextLabel, iface)
    !                      if(iface .gt. 0) then
    !                          do m = 1, spaceAssigned
    !                              boundaryTransportsOffset = 6 * (m - 1)
    !                              this%BoundaryTransports(boundaryTransportsOffset + iface) =   &
    !                                this%BoundaryTransports(boundaryTransportsOffset + iface) + &
    !                                this%ArrayBufferDbl(m)
    !                          end do
    !                      else            
    !                          do m = 1, spaceAssigned
    !                              if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
    !                                  this%SourceTransports(m) = this%SourceTransports(m) +     &
    !                                    this%ArrayBufferDbl(m)
    !                              else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
    !                                  this%SinkTransports(m) = this%SinkTransports(m) +         &
    !                                    this%ArrayBufferDbl(m)
    !                              end if
    !                          end do
    !                      end if
    !                  end if
    !              else if(header%Method .eq. 2) then
    !                  call this%BudgetReader%FillRecordDataBuffer(header,             &
    !                    this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
    !                    status)
    !                  if(spaceAssigned .gt. 0) then
    !                      call this%CheckForDefaultIface(header%TextLabel, iface)
    !                      if(iface .gt. 0) then
    !                          do m = 1, spaceAssigned
    !                              cellNumber = this%ListItemBuffer(m)%CellNumber
    !                              boundaryTransportsOffset = 6 * (cellNumber - 1)
    !                              this%BoundaryTransports(boundaryTransportsOffset + iface) =   &
    !                                this%BoundaryTransports(boundaryTransportsOffset + iface) + &
    !                                this%ListItemBuffer(m)%BudgetValue
    !                          end do
    !                      else            
    !                          do m = 1, spaceAssigned
    !                              cellNumber = this%ListItemBuffer(m)%CellNumber
    !                              if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
    !                                  this%SourceTransports(cellNumber) =                  &
    !                                    this%SourceTransports(cellNumber) +                &
    !                                    this%ListItemBuffer(m)%BudgetValue
    !                              else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
    !                                  this%SinkTransports(cellNumber) =                    &
    !                                    this%SinkTransports(cellNumber) +                  &
    !                                    this%ListItemBuffer(m)%BudgetValue
    !                              end if
    !                          end do
    !                      end if
    !                  end if
    !              else if((header%Method .eq. 5) .or. (header%Method .eq. 6)) then
    !                  call this%BudgetReader%FillRecordDataBuffer(header,             &
    !                    this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
    !                    status)
    !                  if(spaceAssigned .gt. 0) then
    !                      do m = 1, spaceAssigned
    !                          call this%CheckForDefaultIface(header%TextLabel, iface)
    !                          index = header%FindAuxiliaryNameIndex('IFACE')
    !                          if(index .gt. 0) then
    !                              iface = int(this%ListItemBuffer(m)%AuxiliaryValues(index))
    !                          end if
    !                          
    !                          cellNumber = this%ListItemBuffer(m)%CellNumber
    !                          if(iface .gt. 0) then
    !                              boundaryTransportsOffset = 6 * (cellNumber - 1)
    !                              this%BoundaryTransports(boundaryTransportsOffset + iface) =   &
    !                                this%BoundaryTransports(boundaryTransportsOffset + iface) + &
    !                                this%ListItemBuffer(m)%BudgetValue
    !                          else
    !                              if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
    !                                  this%SourceTransports(cellNumber) =                  &
    !                                    this%SourceTransports(cellNumber) +                &
    !                                    this%ListItemBuffer(m)%BudgetValue
    !                              else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
    !                                  this%SinkTransports(cellNumber) =                    &
    !                                    this%SinkTransports(cellNumber) +                  &
    !                                    this%ListItemBuffer(m)%BudgetValue
    !                              end if
    !                          end if
    !                      end do
    !                  end if
    !              end if
    !         
    !         end select
    !         
    !    end do
    !
    !    this%CurrentStressPeriod = stressPeriod
    !    this%CurrentTimeStep = timeStep
    !
    !end subroutine pr_LoadTimeStep


    !subroutine pr_ClearTimeStepBudgetData(this)
    !!***************************************************************************************************************
    !!
    !!***************************************************************************************************************
    !!
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !implicit none
    !class(TransportModelDataType) :: this
    !integer :: cellCount, n, arraySize
    !!---------------------------------------------------------------------------------------------------------------
    !  
    !    this%CurrentStressPeriod = 0
    !    this%CurrentTimeStep = 0
    !    
    !    if(allocated(this%SinkTransports)) then
    !        cellCount = this%Grid%CellCount
    !        do n = 1, cellCount
    !            this%IBoundTS(n) = this%IBound(n)
    !            this%Heads(n) = 0.0d0
    !            this%SourceTransports(n) = 0.0d0
    !            this%SinkTransports(n) = 0.0d0
    !            this%StorageTransports(n) = 0.0d0
    !            this%SubFaceTransportsComputed(n) = .false.
    !        end do
    !        
    !        arraySize = cellCount * 6
    !        do n = 1, arraySize
    !            this%BoundaryTransports(n) = 0.0d0
    !        end do
    !        
    !        arraySize = cellCount * 4
    !        do n = 1, arraySize
    !            this%SubFaceTransports(n) = 0.0d0
    !        end do
    !    
    !        arraySize = this%BudgetReader%GetTransportArraySize()
    !        if(this%Grid%GridType .eq. 1) then
    !            do n = 1, arraySize
    !            this%TransportsRightFace(n) = 0.0d0
    !            this%TransportsFrontFace(n) = 0.0d0
    !            this%TransportsLowerFace(n) = 0.0d0
    !            end do
    !        else if(this%Grid%GridType .eq. 2) then
    !            do n = 1, arraySize
    !                this%TransportsJA(n) = 0.0d0
    !            end do
    !        end if
    !       
    !    end if
 

    !end subroutine pr_ClearTimeStepBudgetData


    !function pr_GetCurrentStressPeriod(this) result(stressPeriod)
    !!***************************************************************************************************************
    !!
    !!***************************************************************************************************************
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !implicit none
    !class(TransportModelDataType) :: this
    !integer :: stressPeriod
    !!---------------------------------------------------------------------------------------------------------------
    ! 

    !    stressPeriod = this%CurrentStressPeriod
    ! 

    !end function pr_GetCurrentStressPeriod


    !function pr_GetCurrentTimeStep(this) result(timeStep)
    !!***************************************************************************************************************
    !!
    !!***************************************************************************************************************
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !implicit none
    !class(TransportModelDataType) :: this
    !integer :: timeStep
    !!---------------------------------------------------------------------------------------------------------------
    !

    !    timeStep = this%CurrentTimeStep
    !  

    !end function pr_GetCurrentTimeStep



    !subroutine pr_SetPorosity(this, porosity, arraySize)
    !!***************************************************************************************************************
    !!
    !!***************************************************************************************************************
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !implicit none
    !class(TransportModelDataType) :: this
    !integer,intent(in) :: arraySize
    !doubleprecision,dimension(arraySize),intent(in),target :: porosity
    !!---------------------------------------------------------------------------------------------------------------


    !    if(arraySize .ne. this%Grid%CellCount) then
    !        write(*,*) "TransportModelDataType: The Porosity array size does not match the cell count for the grid. stop"
    !        stop
    !    end if
    !    
    !    this%Porosity => porosity


    !end subroutine pr_SetPorosity


    !subroutine pr_SetRetardation(this, retardation, arraySize)
    !!***************************************************************************************************************
    !!
    !!***************************************************************************************************************
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !implicit none
    !class(TransportModelDataType) :: this
    !integer,intent(in) :: arraySize
    !doubleprecision,dimension(arraySize),intent(in),target :: retardation
    !!---------------------------------------------------------------------------------------------------------------
  

    !    if(arraySize .ne. this%Grid%CellCount) then
    !        write(*,*) "TransportModelDataType: The Retardation array size does not match the cell count for the grid. stop"
    !        stop
    !    end if
    !    
    !    this%Retardation => retardation
 

    !end subroutine pr_SetRetardation

  
    !subroutine pr_SetDefaultIface(this, defaultIfaceLabels, defaultIfaceValues, arraySize)
    !!***************************************************************************************************************
    !! Description goes here
    !!***************************************************************************************************************
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !use UtilMiscModule,only : utrimall
    !implicit none
    !class(TransportModelDataType) :: this
    !integer,intent(in) :: arraySize
    !integer,dimension(arraySize),intent(in) :: defaultIfaceValues
    !character(len=16),dimension(arraySize),intent(in) :: defaultIfaceLabels
    !integer :: n, firstNonBlank, lastNonBlank, trimmedLength
    !character(len=16) :: label
    !!---------------------------------------------------------------------------------------------------------------
    !  

    !    this%DefaultIfaceCount = 0
    !    do n = 1, 20
    !        this%DefaultIfaceValues(n) = 0
    !        this%DefaultIfaceLabels(n) = '                '
    !    end do
    !    
    !    do n = 1, arraySize
    !        this%DefaultIfaceValues(n) = defaultIfaceValues(n)
    !        label = defaultIfaceLabels(n)
    !        call utrimall(label)
    !        this%DefaultIfaceLabels(n) = label
    !    end do
    !    this%DefaultIfaceCount = arraySize
   

    !end subroutine pr_SetDefaultIface


    !subroutine pr_CheckForDefaultIface(this, textLabel, iface)
    !!***************************************************************************************************************
    !!
    !!***************************************************************************************************************
    !!
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !use UtilMiscModule,only : utrimall
    !implicit none
    !class(TransportModelDataType) :: this
    !character*(*), intent(in) :: textLabel
    !integer,intent(inout) :: iface
    !integer :: n
    !character(len=16) :: label
    !!---------------------------------------------------------------------------------------------------------------
    !  
    !    iface = 0
    !    label = textLabel
    !    call utrimall(label)
    !    do n = 1, this%DefaultIfaceCount
    !        if(label .eq. this%DefaultIfaceLabels(n)) then
    !            iface = this%DefaultIfaceValues(n)
    !            return
    !        end if
    !    end do
    !  
    !end subroutine pr_CheckForDefaultIface
