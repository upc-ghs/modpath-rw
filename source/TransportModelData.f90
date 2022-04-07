module TransportModelDataModule
  use ModflowRectangularGridModule,only : ModflowRectangularGridType
  use ParticleTrackingOptionsModule,only : ParticleTrackingOptionsType
  use ParticleGroupModule,only : ParticleGroupType
  use StartingLocationReaderModule,only : ReadAndPrepareLocations
  implicit none
  !---------------------------------------------------------------------------------------------------------------

  ! Set default access status to private
  private

    type,public :: TransportModelDataType

      logical :: Initialized = .false.
      doubleprecision,dimension(:),allocatable :: AlphaLong
      doubleprecision,dimension(:),allocatable :: AlphaTrans
      integer,allocatable,dimension(:)         :: ICBound
      integer,allocatable,dimension(:)         :: ICBoundTS
      doubleprecision :: DMol

      ! grid
      class(ModflowRectangularGridType),pointer :: Grid => null()

      ! Private variables
      integer,private :: CurrentStressPeriod = 0
      integer,private :: CurrentTimeStep = 0

    contains

      procedure :: Initialize=>pr_Initialize
      procedure :: Reset=>pr_Reset
      procedure :: ReadData=>pr_ReadData
      procedure :: LoadTimeStep=>pr_LoadTimeStep

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

        if(allocated(this%AlphaLong)) deallocate(this%AlphaLong)
        if(allocated(this%AlphaTrans)) deallocate(this%AlphaTrans)
        if(allocated(this%ICBoundTS)) deallocate(this%ICBoundTS)
        if(allocated(this%ICBound)) deallocate(this%ICBound)
        allocate(this%AlphaTrans(grid%CellCount))
        allocate(this%AlphaLong(grid%CellCount))
        allocate(this%ICBoundTS(grid%CellCount))
        allocate(this%ICBound(grid%CellCount))


        this%Initialized = .true.


    end subroutine pr_Initialize


    subroutine pr_ReadData(this, inUnit, inFile, outUnit, ibound, grid, trackingOptions )
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
    class(ModflowRectangularGridType),intent(in) :: grid
    integer,dimension(grid%CellCount),intent(in) :: ibound
    character(len=200) :: line
    character(len=16) :: txt
    integer :: n, m, nn, length, iface, errorCode, layer
    integer :: icol, istart, istop
    integer :: iodispersion = 0
    doubleprecision :: r
    character(len=24),dimension(4) :: aname
    data aname(1) /'          BOUNDARY ARRAY'/
    data aname(2) /'            DISPERSIVITY'/
    data aname(3) /'            ICBOUND'/
    data aname(4) /'            IC'/
    integer :: tempAlphaUnit = 666
    character(len=200) :: tempAlphaFile

    ! Initial conditions
    integer :: nInitialConditions, particleCount, seqNumber, slocUnit
    integer :: releaseOption, releaseTimeCount
    doubleprecision :: initialReleaseTime, releaseInterval
    doubleprecision,dimension(:),allocatable :: releaseTimes
    doubleprecision :: frac, tinc
    type(ParticleGroupType),dimension(:),allocatable :: particleGroups

    !---------------------------------------------------------------------------------------------------------------

        
        ! Open dispersion unit 
        open( inUnit, file=inFile, status='old', access='sequential')


        ! Required for u3d
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


        ! ICBOUND
        if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
            call u3dintmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
              grid%ColumnCount, grid%CellCount, this%ICBound, ANAME(3))
        else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
            call u3dintmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,  &
              this%ICBound, aname(3), cellsPerLayer)
        else
              write(outUnit,*) 'Invalid grid type specified when reading ICBOUND array data.'
              write(outUnit,*) 'Stopping.'
              call ustop(' ')          
        end if


        ! INITIAL CONDITIONS
        ! MORE LIKE SPECIES AND EACH WITH AN INITIAL DISTRIBUTION
        ! THESE ARE TRANSFORMED INTO MASS PARTICLES 
        read(inUnit, *) nInitialConditions
        write(outUnit,'(/A,I5)') 'Number of initial conditions = ', nInitialConditions
        particleCount = 0
        slocUnit = 0
        seqNumber = 0
        ! ADD THEM TO EXISTING PARTICLE GROUPS ?
        if(nInitialConditions .gt. 0) then

            ! Would be like a realloc to extend 
            allocate(particleGroups(nInitialConditions))

            ! Loop over initial conditions
            do n = 1, nInitialConditions

                particleGroups(n)%Group = n


                read(inUnit, '(a)') particleGroups(n)%Name

                ! IDEA
                !read(inUnit, *) initialConditionFormat 1: concentration, 2: particles (classic)

                ! select case ( initialConditionFormat )
                !   case (1) ! Read concentrations in porosity format + mass parameters and transform to particles 
                
                        !! READ AS DENSITY/CONCENTRATION
                        !if((grid%GridType .eq. 1) .or. (grid%GridType .eq. 3)) then
                        !    call u3ddblmp(inUnit, outUnit, grid%LayerCount, grid%RowCount,      &
                        !      grid%ColumnCount, grid%CellCount, this%ICBound, ANAME(3))
                        !else if((grid%GridType .eq. 2) .or. (grid%GridType .eq. 4)) then
                        !    call u3ddblmpusg(inUnit, outUnit, grid%CellCount, grid%LayerCount,  &
                        !      this%ICBound, aname(3), cellsPerLayer)
                        !else
                        !      write(outUnit,*) 'Invalid grid type specified when reading ICBOUND array data.'
                        !      write(outUnit,*) 'Stopping.'
                        !      call ustop(' ')          
                        !end if

                !   case (2) ! Read particles, the classical way + mass parameters
                     
                        ! releaseOption
                        read(inUnit, *) releaseOption
                        select case (releaseOption)
                            case (1)
                                read(inUnit, *) initialReleaseTime
                                call particleGroups(n)%SetReleaseOption1(initialReleaseTime)
                            case (2)
                                read(inUnit, *) releaseTimeCount, initialReleaseTime, releaseInterval
                                call particleGroups(n)%SetReleaseOption2(initialReleaseTime, &
                                  releaseTimeCount, releaseInterval)
                            case (3)
                                read(inUnit, *) releaseTimeCount
                                if(allocated(releaseTimes)) deallocate(releaseTimes)
                                allocate(releaseTimes(releaseTimeCount))
                                read(inUnit, *) (releaseTimes(nn), nn = 1, releaseTimeCount)
                                call particleGroups(n)%SetReleaseOption3(releaseTimeCount,   &
                                  releaseTimes)
                            case default
                            ! write error message and stop
                        end select
                        
                        read(inUnit, '(a)') line
                        icol = 1
                        call urword(line,icol,istart,istop,1,n,r,0,0)
                        if(line(istart:istop) .eq. 'EXTERNAL') then
                            call urword(line,icol,istart,istop,0,n,r,0,0)
                            particleGroups(n)%LocationFile = line(istart:istop)
                            slocUnit = 0
                        else if(line(istart:istop) .eq. 'INTERNAL') then
                            particleGroups(n)%LocationFile = ''
                            slocUnit = inUnit
                        else
                            call ustop('Invalid starting locations file name. stop.')
                        end if

                ! end select

                ! THE LOCATION FUNCTION/DISTRIBUTION
                call ReadAndPrepareLocations(slocUnit, outUnit, particleGroups(n),   &
                  ibound, grid%CellCount, grid, seqNumber)

                ! THE HOLDER OF PARTICLES
                write(outUnit, '(a,i4,a,i10,a)') 'Initial condition ', n, ' contains ',   &
                  particleGroups(n)%TotalParticleCount, ' particles.'
                particleCount = particleCount + particleGroups(n)%TotalParticleCount
            end do

            write(outUnit, '(a,i10)') 'Total number of particles on initial conditions = ', particleCount
            write(outUnit, *)

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
        if(allocated(this%ICBoundTS)) deallocate( this%ICBoundTS )
        if(allocated(this%ICBound)) deallocate( this%ICBound )

    end subroutine pr_Reset



    subroutine pr_LoadTimeStep(this, stressPeriod, timeStep)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(TransportModelDataType) :: this
    integer,intent(in) :: stressPeriod, timeStep
    integer :: firstRecord, lastRecord, n, m, firstNonBlank, lastNonBlank, &
      trimmedLength
    integer :: spaceAssigned, status,cellCount, iface, index,              &
      boundaryFlowsOffset, listItemBufferSize, cellNumber, layer
    !type(BudgetRecordHeaderType) :: header
    character(len=16) :: textLabel
    doubleprecision :: top 
    real :: HDryTol, HDryDiff
    !---------------------------------------------------------------------------------------------------------------
     

        !! Readers of time variable data 
        !call this%ClearTimeStepBudgetData()
        !call this%BudgetReader%GetRecordHeaderRange(stressPeriod, timeStep, firstRecord, lastRecord)
        !if(firstRecord .eq. 0) return
    
        cellCount = this%Grid%CellCount
        !listItemBufferSize = size(this%ListItemBuffer)
        
        ! Set steady state = true, then change it if the budget file contains storage
        !this%SteadyState = .true.
        
        ! MASS TRANSPORT; RWPT: could be recycled for changes in boundary conditions
        ! Load heads for this time step
        !call this%HeadReader%FillTimeStepHeadBuffer(stressPeriod, timeStep, &
        !  this%Heads, cellCount, spaceAssigned)
        
        ! Fill ICBoundTS array and set the SaturatedTop array for the Grid.
        ! The saturated top is set equal to the top for confined cells and water table cells 
        ! where the head is above the top or below the bottom.
        !HDryTol = abs(epsilon(HDryTol)*sngl(this%HDry))
        if(this%Grid%GridType .gt. 2) then
            do n = 1, cellCount
                this%ICBoundTS(n) = this%ICBound(n)
                ! DO SOMETHING WITH IBOUNDTS 
                ! DEPENDNING ON HEADS OR OTHERS
                !layer = this%Grid%GetLayer(n)
                !if(this%Grid%CellType(n) .eq. 1) then
                !    HDryDiff = sngl(this%Heads(n)) - sngl(this%HDry)
                !    if(abs(HDryDiff) .lt. HDryTol) then
                !        this%IBoundTS(n) = 0
                !        if(this%Heads(n) .lt. this%Grid%Bottom(n)) then
                !            this%IBoundTS(n) = 0
                !            this%Grid%SaturatedTop(n) = this%Grid%Bottom(n)
                !        end if
                !    end if
                !    if(this%IBoundTS(n) .ne. 0) then
                !        if((this%Heads(n) .le. this%Grid%Top(n)) .and. &
                !          (this%Heads(n) .ge. this%Grid%Bottom(n))) then
                !            this%Grid%SaturatedTop(n) = this%Heads(n)
                !        end if
                !    end if
                !end if
            end do
            
        else
            do n = 1, cellCount
                this%ICBoundTS(n) = this%ICBound(n)
                !this%Grid%SaturatedTop(n) = this%Grid%Top(n)
                !this%StorageFlows(n) = 0.0
                !this%IBoundTS(n) = this%IBound(n)
                !layer = this%Grid%GetLayer(n)
                !if(this%Grid%CellType(n) .eq. 1) then
                !    HDryDiff = sngl(this%Heads(n)) - sngl(this%HDry)
                !    if((abs(HDryDiff) .lt. HDryTol) .or. (this%Heads(n) .gt. 1.0d+6)) then
                !        this%IBoundTS(n) = 0
                !    end if
                !    if(this%IBoundTS(n) .ne. 0) then
                !        if((this%Heads(n) .le. this%Grid%Top(n)) .and. &
                !          (this%Heads(n) .ge. this%Grid%Bottom(n))) then
                !            this%Grid%SaturatedTop(n) = this%Heads(n)
                !        end if
                !    end if
                !end if
            end do
        end if
       

        !! MASS TRANSPORT, RWPT: COULD BE USEFUL FOR READING SOMETHING FROM EXTERNAL FILE 
        !! Loop through record headers
        !do n = firstRecord, lastRecord
        !     header = this%BudgetReader%GetRecordHeader(n)
        !     textLabel = header%TextLabel
        !     call TrimAll(textLabel, firstNonBlank, lastNonBlank, trimmedLength)
        !     
        !     select case(textLabel(firstNonBlank:lastNonBlank))
        !     case('CONSTANT HEAD', 'CHD')
        !          ! Read constant head flows into the sinkFlows and sourceFlows arrays.
        !          ! For a standard budget file, Method = 0. For a compact budget file,
        !          ! Method = 2.
        !          if(header%Method .eq. 0) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ArrayBufferDbl, cellCount, spaceAssigned, status)
        !              if(cellCount .eq. spaceAssigned) then
        !                  do m = 1, spaceAssigned
        !                      if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
        !                          this%SourceFlows(m) = this%SourceFlows(m) +         &
        !                            this%ArrayBufferDbl(m)
        !                      else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
        !                          this%SinkFlows(m) = this%SinkFlows(m) +             &
        !                            this%ArrayBufferDbl(m)
        !                      end if
        !                  end do
        !              end if
        !          else if(header%Method .eq. 2) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ListItemBuffer, listItemBufferSize, spaceAssigned, status)
        !              if(spaceAssigned .gt. 0) then
        !                  do m = 1, spaceAssigned
        !                      cellNumber = this%ListItemBuffer(m)%CellNumber
        !                      if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
        !                          this%SourceFlows(cellNumber) =                      &
        !                            this%SourceFlows(cellNumber) + this%ListItemBuffer(m)%BudgetValue
        !                      else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
        !                          this%SinkFlows(cellNumber) =                        &
        !                            this%SinkFlows(cellNumber) + this%ListItemBuffer(m)%BudgetValue
        !                      end if
        !                  end do
        !              end if
        !          else if((header%Method .eq. 5) .or. (header%Method .eq. 6)) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
        !                status)
        !              if(spaceAssigned .gt. 0) then
        !                  do m = 1, spaceAssigned
        !                      call this%CheckForDefaultIface(header%TextLabel, iface)
        !                      index = header%FindAuxiliaryNameIndex('IFACE')
        !                      if(index .gt. 0) then
        !                          iface = int(this%ListItemBuffer(m)%AuxiliaryValues(index))
        !                      end if
        !                      
        !                      cellNumber = this%ListItemBuffer(m)%CellNumber
        !                      if(iface .gt. 0) then
        !                          boundaryFlowsOffset = 6 * (cellNumber - 1)
        !                          this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
        !                            this%BoundaryFlows(boundaryFlowsOffset + iface) + &
        !                            this%ListItemBuffer(m)%BudgetValue
        !                      else
        !                          if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
        !                              this%SourceFlows(cellNumber) =                  &
        !                                this%SourceFlows(cellNumber) +                &
        !                                this%ListItemBuffer(m)%BudgetValue
        !                          else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
        !                              this%SinkFlows(cellNumber) =                    &
        !                                this%SinkFlows(cellNumber) +                  &
        !                                this%ListItemBuffer(m)%BudgetValue
        !                          end if
        !                      end if
        !                  end do
        !              end if
        !          end if
        !         
        !     case('STORAGE', 'STO-SS', 'STO-SY')
        !          ! Read storage for all cells into the StorageFlows array.
        !          ! Method should always be 0 or 1, but check anyway to be sure.
        !          if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
        !              if(header%ArrayItemCount .eq. cellCount) then
        !                  call this%BudgetReader%FillRecordDataBuffer(header,         &
        !                    this%ArrayBufferDbl, cellCount, spaceAssigned, status)
        !                  if(cellCount .eq. spaceAssigned) then
        !                      do m = 1, spaceAssigned
        !                          this%StorageFlows(m) = this%StorageFlows(m) + this%ArrayBufferDbl(m)
        !                          if(this%StorageFlows(m) .ne. 0.0) this%SteadyState = .false.
        !                      end do
        !                  end if
        !              end if
        !          end if
        !         
        !     case('FLOW JA FACE', 'FLOW-JA-FACE')
        !          ! Read connected face flows into the FlowsJA array for unstructured grids.
        !          if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
        !              ! Method should always be 0 or 1 for flow between grid cells. 
        !              if(header%ArrayItemCount .eq. this%BudgetReader%GetFlowArraySize()) then
        !                  call this%BudgetReader%FillRecordDataBuffer(header,         &
        !                    this%FlowsJA, header%ArrayItemCount, spaceAssigned,       &
        !                    status)
        !              end if
        !          else if(header%Method .eq. 6) then
        !              ! Method code 6 indicates flow to or from cells in the current model grid
        !              ! and another connected model grid in a multi-model MODFLOW-6 simulation. 
        !              ! Treat flows to or from connected model grids as distributed source/sink flows 
        !              ! for the current grid.
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
        !                status)
        !              if(spaceAssigned .gt. 0) then
        !                  do m = 1, spaceAssigned
        !                      cellNumber = this%ListItemBuffer(m)%CellNumber
        !                      if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
        !                          this%SourceFlows(cellNumber) =                  &
        !                              this%SourceFlows(cellNumber) +                &
        !                              this%ListItemBuffer(m)%BudgetValue
        !                      else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
        !                          this%SinkFlows(cellNumber) =                    &
        !                              this%SinkFlows(cellNumber) +                  &
        !                              this%ListItemBuffer(m)%BudgetValue
        !                      end if
        !                  end do
        !              end if
        !          end if
        !         
        !     case('FLOW RIGHT FACE')
        !          ! Read flows across the right face for structured grids.
        !          ! Method should always be 0 or 1, but check anyway to be sure.
        !          if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
        !              if(header%ArrayItemCount .eq. this%BudgetReader%GetFlowArraySize()) then
        !                  call this%BudgetReader%FillRecordDataBuffer(header,         &
        !                    this%FlowsRightFace, header%ArrayItemCount, spaceAssigned,&
        !                    status)
        !              end if
        !          end if
        !         
        !     case('FLOW FRONT FACE')
        !          ! Read flows across the front face for structured grids.
        !          ! Method should always be 0 or 1, but check anyway to be sure.
        !          if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
        !              if(header%ArrayItemCount .eq. this%BudgetReader%GetFlowArraySize()) then
        !                  call this%BudgetReader%FillRecordDataBuffer(header,         &
        !                    this%FlowsFrontFace, header%ArrayItemCount, spaceAssigned,&
        !                    status)
        !              end if
        !          end if
        !         
        !     case('FLOW LOWER FACE')
        !          ! Read flows across the lower face for structured grids.
        !          ! Method should always be 0 or 1, but check anyway to be sure.
        !          if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
        !              if(header%ArrayItemCount .eq. this%BudgetReader%GetFlowArraySize()) then
        !                  call this%BudgetReader%FillRecordDataBuffer(header,         &
        !                    this%FlowsLowerFace, header%ArrayItemCount, spaceAssigned,&
        !                    status)
        !              end if
        !          end if
        !     
        !      case default
        !          ! Now handle any other records in the budget file.
        !           if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
        !              if(header%ArrayItemCount .eq. cellCount) then
        !                  call this%BudgetReader%FillRecordDataBuffer(header,         &
        !                    this%ArrayBufferDbl, cellCount, spaceAssigned, status)
        !                  if(cellCount .eq. spaceAssigned) then
        !                      call this%CheckForDefaultIface(header%TextLabel, iface)
        !                      if(iface .gt. 0) then
        !                          do m = 1, spaceAssigned
        !                              boundaryFlowsOffset = 6 * (m - 1)
        !                              this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
        !                                this%BoundaryFlows(boundaryFlowsOffset + iface) + &
        !                                this%ArrayBufferDbl(m)
        !                          end do
        !                      else
        !                          do m = 1, spaceAssigned
        !                              if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
        !                                  this%SourceFlows(m) = this%SourceFlows(m) +     &
        !                                    this%ArrayBufferDbl(m)
        !                              else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
        !                                  this%SinkFlows(m) = this%SinkFlows(m) +         &
        !                                    this%ArrayBufferDbl(m)
        !                              end if
        !                          end do
        !                      end if
        !                  end if
        !              end if
        !           else if(header%Method .eq. 3) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ArrayBufferDbl, this%ArrayBufferInt,                     &
        !                header%ArrayItemCount, spaceAssigned, status)
        !              if(header%ArrayItemCount .eq. spaceAssigned) then
        !                  call this%CheckForDefaultIface(header%TextLabel, iface)
        !                  if(iface .gt. 0) then
        !                      do m = 1, spaceAssigned
        !                          cellNumber = this%ArrayBufferInt(m)
        !                          boundaryFlowsOffset = 6 * (cellNumber - 1)
        !                          this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
        !                            this%BoundaryFlows(boundaryFlowsOffset + iface) + &
        !                            this%ArrayBufferDbl(m)
        !                      end do
        !                  else            
        !                      do m = 1, spaceAssigned
        !                          cellNumber = this%ArrayBufferInt(m)
        !                          if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
        !                              this%SourceFlows(cellNumber) =                  &
        !                                this%SourceFlows(cellNumber) +                &
        !                                this%ArrayBufferDbl(m)
        !                          else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
        !                              this%SinkFlows(cellNumber) =                    &
        !                                this%SinkFlows(cellNumber) +                  &
        !                                this%ArrayBufferDbl(m)
        !                          end if
        !                      end do
        !                  end if
        !              end if
        !           else if(header%Method .eq. 4) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ArrayBufferDbl, header%ArrayItemCount, spaceAssigned,    &
        !                status)
        !              if(header%ArrayItemCount .eq. spaceAssigned) then
        !                  call this%CheckForDefaultIface(header%TextLabel, iface)
        !                  if(iface .gt. 0) then
        !                      do m = 1, spaceAssigned
        !                          boundaryFlowsOffset = 6 * (m - 1)
        !                          this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
        !                            this%BoundaryFlows(boundaryFlowsOffset + iface) + &
        !                            this%ArrayBufferDbl(m)
        !                      end do
        !                  else            
        !                      do m = 1, spaceAssigned
        !                          if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
        !                              this%SourceFlows(m) = this%SourceFlows(m) +     &
        !                                this%ArrayBufferDbl(m)
        !                          else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
        !                              this%SinkFlows(m) = this%SinkFlows(m) +         &
        !                                this%ArrayBufferDbl(m)
        !                          end if
        !                      end do
        !                  end if
        !              end if
        !          else if(header%Method .eq. 2) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
        !                status)
        !              if(spaceAssigned .gt. 0) then
        !                  call this%CheckForDefaultIface(header%TextLabel, iface)
        !                  if(iface .gt. 0) then
        !                      do m = 1, spaceAssigned
        !                          cellNumber = this%ListItemBuffer(m)%CellNumber
        !                          boundaryFlowsOffset = 6 * (cellNumber - 1)
        !                          this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
        !                            this%BoundaryFlows(boundaryFlowsOffset + iface) + &
        !                            this%ListItemBuffer(m)%BudgetValue
        !                      end do
        !                  else            
        !                      do m = 1, spaceAssigned
        !                          cellNumber = this%ListItemBuffer(m)%CellNumber
        !                          if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
        !                              this%SourceFlows(cellNumber) =                  &
        !                                this%SourceFlows(cellNumber) +                &
        !                                this%ListItemBuffer(m)%BudgetValue
        !                          else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
        !                              this%SinkFlows(cellNumber) =                    &
        !                                this%SinkFlows(cellNumber) +                  &
        !                                this%ListItemBuffer(m)%BudgetValue
        !                          end if
        !                      end do
        !                  end if
        !              end if
        !          else if((header%Method .eq. 5) .or. (header%Method .eq. 6)) then
        !              call this%BudgetReader%FillRecordDataBuffer(header,             &
        !                this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
        !                status)
        !              if(spaceAssigned .gt. 0) then
        !                  do m = 1, spaceAssigned
        !                      call this%CheckForDefaultIface(header%TextLabel, iface)
        !                      index = header%FindAuxiliaryNameIndex('IFACE')
        !                      if(index .gt. 0) then
        !                          iface = int(this%ListItemBuffer(m)%AuxiliaryValues(index))
        !                      end if
        !                      
        !                      cellNumber = this%ListItemBuffer(m)%CellNumber
        !                      if(iface .gt. 0) then
        !                          boundaryFlowsOffset = 6 * (cellNumber - 1)
        !                          this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
        !                            this%BoundaryFlows(boundaryFlowsOffset + iface) + &
        !                            this%ListItemBuffer(m)%BudgetValue
        !                      else
        !                          if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
        !                              this%SourceFlows(cellNumber) =                  &
        !                                this%SourceFlows(cellNumber) +                &
        !                                this%ListItemBuffer(m)%BudgetValue
        !                          else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
        !                              this%SinkFlows(cellNumber) =                    &
        !                                this%SinkFlows(cellNumber) +                  &
        !                                this%ListItemBuffer(m)%BudgetValue
        !                          end if
        !                      end if
        !                  end do
        !              end if
        !          end if
        !     
        !     end select
        !     
        !end do
    
        this%CurrentStressPeriod = stressPeriod
        this%CurrentTimeStep = timeStep
    
    end subroutine pr_LoadTimeStep


end module TransportModelDataModule
