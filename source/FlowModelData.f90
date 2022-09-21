module FlowModelDataModule
  use BudgetReaderModule,only : BudgetReaderType
  use HeadReaderModule,only : HeadReaderType
  use BudgetListItemModule,only : BudgetListItemType
  use ModflowRectangularGridModule,only : ModflowRectangularGridType
  use BudgetRecordHeaderModule,only : BudgetRecordHeaderType
  use UtilMiscModule,only : TrimAll
  use UTL8MODULE,only : ustop

  implicit none
  !---------------------------------------------------------------------------------------------------------------

  ! Set default access status to private
  private


    type,public :: FlowModelDataType
      !doubleprecision :: ReferenceTime = 0d0
      !doubleprecision :: StoppingTime = 0d0
      doubleprecision :: HDry = 0d0
      doubleprecision :: HNoFlow = 0d0
      logical :: Initialized = .false.
      logical :: SteadyState = .true.
      integer :: DefaultIfaceCount
      character(len=16),dimension(20) :: DefaultIfaceLabels
      integer,dimension(20) :: DefaultIfaceValues
      integer,allocatable,dimension(:) :: IBoundTS
      integer,allocatable,dimension(:) :: ArrayBufferInt
      doubleprecision,allocatable,dimension(:) :: Heads
      doubleprecision,allocatable,dimension(:) :: FlowsJA
      doubleprecision,allocatable,dimension(:) :: FlowsRightFace
      doubleprecision,allocatable,dimension(:) :: FlowsFrontFace
      doubleprecision,allocatable,dimension(:) :: FlowsLowerFace
      doubleprecision,allocatable,dimension(:) :: SourceFlows
      doubleprecision,allocatable,dimension(:) :: SinkFlows
      doubleprecision,allocatable,dimension(:) :: StorageFlows
      doubleprecision,allocatable,dimension(:) :: BoundaryFlows
      doubleprecision,allocatable,dimension(:) :: SubFaceFlows
      doubleprecision,allocatable,dimension(:) :: ArrayBufferDbl
      ! Externally assigned arrays
      !integer,dimension(:),pointer :: LayerTypes
      integer,dimension(:),pointer :: IBound
      integer,dimension(:),pointer :: Zones
      doubleprecision,dimension(:),pointer :: Porosity
      doubleprecision,dimension(:),pointer :: Retardation
     
      type(HeadReadertype)  , pointer :: HeadReader   => null()
      type(BudgetReaderType), pointer :: BudgetReader => null()
      class(ModflowRectangularGridType),pointer :: Grid => null()
      type(BudgetListItemType),allocatable,dimension(:) :: ListItemBuffer
      logical,allocatable,dimension(:), private :: SubFaceFlowsComputed
      
      ! Private variables
      integer,private :: CurrentStressPeriod = 0
      integer,private :: CurrentTimeStep = 0
    

    contains

      procedure :: Initialize=>pr_Initialize
      procedure :: Reset=>pr_Reset
      procedure :: LoadTimeStep=>pr_LoadTimeStep
      procedure :: ClearTimeStepBudgetData=>pr_ClearTimeStepBudgetData
      procedure :: GetCurrentStressPeriod=>pr_GetCurrentStressPeriod
      procedure :: GetCurrentTimeStep=>pr_GetCurrentTimeStep
      procedure :: SetIBound=>pr_SetIBound
      procedure :: SetZones=>pr_SetZones
      procedure :: SetPorosity=>pr_SetPorosity
      procedure :: SetRetardation=>pr_SetRetardation
      procedure :: SetDefaultIface=>pr_SetDefaultIface
      procedure :: CheckForDefaultIface=>pr_CheckForDefaultIface
      !!procedure :: SetLayerTypes=>pr_SetLayerTypes

    end type


contains


    subroutine pr_Initialize(this,headReader, budgetReader, grid, hNoFlow, hDry)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(FlowModelDataType) :: this
    type(BudgetReaderType),intent(inout),target :: budgetReader
    type(HeadReaderType),intent(inout),target :: headReader
    class(ModflowRectangularGridType),intent(inout),pointer :: grid
    integer :: cellCount,gridType
    integer :: n, flowArraySize
    doubleprecision :: hNoFlow, hDry
    !---------------------------------------------------------------------------------------------------------------


        this%Initialized = .false.
        
        ! Call Reset to make sure that all arrays are initially unallocated
        call this%Reset()
        
        ! Return if the grid cell count equals 0
        cellCount = grid%CellCount
        if(cellCount .le. 0) return
        
        ! Check budget reader and grid data for compatibility and allocate appropriate cell-by-cell flow arrays
        gridType = grid%GridType
        select case (gridType)
            case (1)
                if((budgetReader%GetBudgetType() .ne. 1)) return
                if((headReader%GridStyle .ne. 1) .or. (headReader%CellCount .ne. cellCount)) return
                flowArraySize = budgetReader%GetFlowArraySize()
                if(flowArraySize .ne. cellCount) return
                allocate(this%FlowsRightFace(flowArraySize))
                allocate(this%FlowsFrontFace(flowArraySize))
                allocate(this%FlowsLowerFace(flowArraySize))
                allocate(this%FlowsJA(0))
            case (2)
                if((budgetReader%GetBudgetType() .ne. 2)) return
                if((headReader%GridStyle .ne. 2) .or. (headReader%CellCount .ne. cellCount)) return
                flowArraySize = budgetReader%GetFlowArraySize()
                if(flowArraySize .ne. grid%JaCount) return
                allocate(this%FlowsJA(flowArraySize))
                allocate(this%FlowsRightFace(0))
                allocate(this%FlowsFrontFace(0))
                allocate(this%FlowsLowerFace(0))
            case (3, 4)
                if((budgetReader%GetBudgetType() .ne. 2)) return
                if((headReader%GridStyle .ne. 1) .or. (headReader%CellCount .ne. cellCount)) return
                flowArraySize = budgetReader%GetFlowArraySize()
                if(flowArraySize .ne. grid%JaCount) return
                allocate(this%FlowsJA(flowArraySize))
                allocate(this%FlowsRightFace(0))
                allocate(this%FlowsFrontFace(0))
                allocate(this%FlowsLowerFace(0))          
            !case (4)
            !    ! Not implemented
            !    return
            case default
                return
        end select
        
        ! Set pointers to budgetReader and grid. Assign tracking options.
        this%HeadReader => headReader
        this%BudgetReader => budgetReader
        this%Grid => grid
        this%HNoFlow = hNoFlow
        this%HDry = hDry
        
        ! Allocate the rest of the arrays
        allocate(this%IBoundTS(cellCount))
        allocate(this%Heads(cellCount))
        allocate(this%SourceFlows(cellCount))
        allocate(this%SinkFlows(cellCount))
        allocate(this%StorageFlows(cellCount))
        allocate(this%SubFaceFlowsComputed(cellCount))
        allocate(this%BoundaryFlows(cellCount * 6))
        allocate(this%SubFaceFlows(cellCount * 4))
        
        ! Allocate buffers for reading array and list data
        allocate(this%ListItemBuffer(this%BudgetReader%GetMaximumListItemCount()))
        allocate(this%ArrayBufferDbl(this%BudgetReader%GetMaximumArrayItemCount()))  
        allocate(this%ArrayBufferInt(this%BudgetReader%GetMaximumArrayItemCount()))
        
        this%Initialized = .true.
  

    end subroutine pr_Initialize


    subroutine pr_Reset(this)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(FlowModelDataType) :: this
    !---------------------------------------------------------------------------------------------------------------
       
        !this%ReferenceTime = 0.0d0
        !this%StoppingTime = 0.0d0
        this%CurrentStressPeriod = 0
        this%CurrentTimeStep = 0
        this%HeadReader => null()
        this%BudgetReader => null()
        this%Grid => null()
        
        if(allocated(this%IBoundTS)) deallocate(this%IBoundTS)
        if(allocated(this%ArrayBufferInt)) deallocate(this%ArrayBufferInt)
        if(allocated(this%Heads)) deallocate(this%Heads)
        if(allocated(this%FlowsJA)) deallocate(this%FlowsJA)
        if(allocated(this%FlowsRightFace)) deallocate(this%FlowsRightFace)
        if(allocated(this%FlowsFrontFace)) deallocate(this%FlowsFrontFace)
        if(allocated(this%FlowsLowerFace)) deallocate(this%FlowsLowerFace)
        if(allocated(this%SourceFlows)) deallocate(this%SourceFlows)
        if(allocated(this%SinkFlows)) deallocate(this%SinkFlows)
        if(allocated(this%StorageFlows)) deallocate(this%StorageFlows)
        if(allocated(this%BoundaryFlows)) deallocate(this%BoundaryFlows)
        if(allocated(this%SubFaceFlows)) deallocate(this%SubFaceFlows)
        if(allocated(this%ArrayBufferDbl)) deallocate(this%ArrayBufferDbl)
        if(allocated(this%ListItemBuffer)) deallocate(this%ListItemBuffer)
        if(allocated(this%SubFaceFlowsComputed)) deallocate(this%SubFaceFlowsComputed)
        this%IBound => null()
        this%Porosity => null()
        this%Retardation => null()
        this%Zones => null()


    end subroutine pr_Reset


    subroutine pr_LoadTimeStep(this, stressPeriod, timeStep)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(FlowModelDataType) :: this
    integer,intent(in) :: stressPeriod, timeStep
    integer :: firstRecord, lastRecord, n, m, firstNonBlank, lastNonBlank, &
      trimmedLength
    integer :: spaceAssigned, status,cellCount, iface, index,              &
      boundaryFlowsOffset, listItemBufferSize, cellNumber, layer
    type(BudgetRecordHeaderType) :: header
    character(len=16) :: textLabel
    character(len=132) message
    doubleprecision :: top 
    real :: HDryTol, HDryDiff
    !---------------------------------------------------------------------------------------------------------------
      
        call this%ClearTimeStepBudgetData()
        call this%BudgetReader%GetRecordHeaderRange(stressPeriod, timeStep, firstRecord, lastRecord)
        if(firstRecord .eq. 0) then
          write(message,'(A,I5,A,I5,A)') ' Error loading Time Step ', timeStep, ' Period ', stressPeriod, '.'
          message = trim(message)
          write(*,'(A)') message
          call ustop('Missing budget information. Budget file must have output for every time step. Stop.')
        end if

        cellCount = this%Grid%CellCount
        listItemBufferSize = size(this%ListItemBuffer)
        
        ! Set steady state = true, then change it if the budget file contains storage
        this%SteadyState = .true.
        
        ! Load heads for this time step
        call this%HeadReader%FillTimeStepHeadBuffer(stressPeriod, timeStep, &
          this%Heads, cellCount, spaceAssigned)
        
        ! Fill IBoundTS array and set the SaturatedTop array for the Grid.
        ! The saturated top is set equal to the top for confined cells and water table cells 
        ! where the head is above the top or below the bottom.
        HDryTol = abs(epsilon(HDryTol)*sngl(this%HDry))
        if(this%Grid%GridType .gt. 2) then
            do n = 1, cellCount
                this%Grid%SaturatedTop(n) = this%Grid%Top(n)
                this%StorageFlows(n) = 0.0
                this%IBoundTS(n) = this%IBound(n)
                layer = this%Grid%GetLayer(n)
                if(this%Grid%CellType(n) .eq. 1) then
                    HDryDiff = sngl(this%Heads(n)) - sngl(this%HDry)
                    if(abs(HDryDiff) .lt. HDryTol) then
                        this%IBoundTS(n) = 0
                        if(this%Heads(n) .lt. this%Grid%Bottom(n)) then
                            this%IBoundTS(n) = 0
                            this%Grid%SaturatedTop(n) = this%Grid%Bottom(n)
                        end if
                    end if
                    if(this%IBoundTS(n) .ne. 0) then
                        if((this%Heads(n) .le. this%Grid%Top(n)) .and. &
                          (this%Heads(n) .ge. this%Grid%Bottom(n))) then
                            this%Grid%SaturatedTop(n) = this%Heads(n)
                        end if
                    end if
                end if
            end do
            
        else
            do n = 1, cellCount
                this%Grid%SaturatedTop(n) = this%Grid%Top(n)
                this%StorageFlows(n) = 0.0
                this%IBoundTS(n) = this%IBound(n)
                layer = this%Grid%GetLayer(n)
                if(this%Grid%CellType(n) .eq. 1) then
                    HDryDiff = sngl(this%Heads(n)) - sngl(this%HDry)
                    if((abs(HDryDiff) .lt. HDryTol) .or. (this%Heads(n) .gt. 1.0d+6)) then
                        this%IBoundTS(n) = 0
                    end if
                    if(this%IBoundTS(n) .ne. 0) then
                        if((this%Heads(n) .le. this%Grid%Top(n)) .and. &
                          (this%Heads(n) .ge. this%Grid%Bottom(n))) then
                            this%Grid%SaturatedTop(n) = this%Heads(n)
                        end if
                    end if
                end if
            end do
        end if
        
        ! Loop through record headers
        do n = firstRecord, lastRecord
             header = this%BudgetReader%GetRecordHeader(n)
             textLabel = header%TextLabel
             call TrimAll(textLabel, firstNonBlank, lastNonBlank, trimmedLength)
             
             select case(textLabel(firstNonBlank:lastNonBlank))
             case('CONSTANT HEAD', 'CHD')
                  ! Read constant head flows into the sinkFlows and sourceFlows arrays.
                  ! For a standard budget file, Method = 0. For a compact budget file,
                  ! Method = 2.
                  if(header%Method .eq. 0) then
                      call this%BudgetReader%FillRecordDataBuffer(header,             &
                        this%ArrayBufferDbl, cellCount, spaceAssigned, status)
                      if(cellCount .eq. spaceAssigned) then
                          do m = 1, spaceAssigned
                              if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
                                  this%SourceFlows(m) = this%SourceFlows(m) +         &
                                    this%ArrayBufferDbl(m)
                              else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
                                  this%SinkFlows(m) = this%SinkFlows(m) +             &
                                    this%ArrayBufferDbl(m)
                              end if
                          end do
                      end if
                  else if(header%Method .eq. 2) then
                      call this%BudgetReader%FillRecordDataBuffer(header,             &
                        this%ListItemBuffer, listItemBufferSize, spaceAssigned, status)
                      if(spaceAssigned .gt. 0) then
                          do m = 1, spaceAssigned
                              cellNumber = this%ListItemBuffer(m)%CellNumber
                              if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
                                  this%SourceFlows(cellNumber) =                      &
                                    this%SourceFlows(cellNumber) + this%ListItemBuffer(m)%BudgetValue
                              else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
                                  this%SinkFlows(cellNumber) =                        &
                                    this%SinkFlows(cellNumber) + this%ListItemBuffer(m)%BudgetValue
                              end if
                          end do
                      end if
                  else if((header%Method .eq. 5) .or. (header%Method .eq. 6)) then
                      call this%BudgetReader%FillRecordDataBuffer(header,             &
                        this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
                        status)
                      if(spaceAssigned .gt. 0) then
                          do m = 1, spaceAssigned
                              call this%CheckForDefaultIface(header%TextLabel, iface)
                              index = header%FindAuxiliaryNameIndex('IFACE')
                              if(index .gt. 0) then
                                  iface = int(this%ListItemBuffer(m)%AuxiliaryValues(index))
                              end if
                              
                              cellNumber = this%ListItemBuffer(m)%CellNumber
                              if(iface .gt. 0) then
                                  boundaryFlowsOffset = 6 * (cellNumber - 1)
                                  this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
                                    this%BoundaryFlows(boundaryFlowsOffset + iface) + &
                                    this%ListItemBuffer(m)%BudgetValue
                              else
                                  if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
                                      this%SourceFlows(cellNumber) =                  &
                                        this%SourceFlows(cellNumber) +                &
                                        this%ListItemBuffer(m)%BudgetValue
                                  else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
                                      this%SinkFlows(cellNumber) =                    &
                                        this%SinkFlows(cellNumber) +                  &
                                        this%ListItemBuffer(m)%BudgetValue
                                  end if
                              end if
                          end do
                      end if
                  end if
                 
             case('STORAGE', 'STO-SS', 'STO-SY')
                  ! Read storage for all cells into the StorageFlows array.
                  ! Method should always be 0 or 1, but check anyway to be sure.
                  if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
                      if(header%ArrayItemCount .eq. cellCount) then
                          call this%BudgetReader%FillRecordDataBuffer(header,         &
                            this%ArrayBufferDbl, cellCount, spaceAssigned, status)
                          if(cellCount .eq. spaceAssigned) then
                              do m = 1, spaceAssigned
                                  this%StorageFlows(m) = this%StorageFlows(m) + this%ArrayBufferDbl(m)
                                  if(this%StorageFlows(m) .ne. 0.0) this%SteadyState = .false.
                              end do
                          end if
                      end if
                  end if
                 
             case('FLOW JA FACE', 'FLOW-JA-FACE')
                  ! Read connected face flows into the FlowsJA array for unstructured grids.
                  if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
                      ! Method should always be 0 or 1 for flow between grid cells. 
                      if(header%ArrayItemCount .eq. this%BudgetReader%GetFlowArraySize()) then
                          call this%BudgetReader%FillRecordDataBuffer(header,         &
                            this%FlowsJA, header%ArrayItemCount, spaceAssigned,       &
                            status)
                      end if
                  else if(header%Method .eq. 6) then
                      ! Method code 6 indicates flow to or from cells in the current model grid
                      ! and another connected model grid in a multi-model MODFLOW-6 simulation. 
                      ! Treat flows to or from connected model grids as distributed source/sink flows 
                      ! for the current grid.
                      call this%BudgetReader%FillRecordDataBuffer(header,             &
                        this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
                        status)
                      if(spaceAssigned .gt. 0) then
                          do m = 1, spaceAssigned
                              cellNumber = this%ListItemBuffer(m)%CellNumber
                              if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
                                  this%SourceFlows(cellNumber) =                  &
                                      this%SourceFlows(cellNumber) +                &
                                      this%ListItemBuffer(m)%BudgetValue
                              else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
                                  this%SinkFlows(cellNumber) =                    &
                                      this%SinkFlows(cellNumber) +                  &
                                      this%ListItemBuffer(m)%BudgetValue
                              end if
                          end do
                      end if
                  end if
                 
             case('FLOW RIGHT FACE')
                  ! Read flows across the right face for structured grids.
                  ! Method should always be 0 or 1, but check anyway to be sure.
                  if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
                      if(header%ArrayItemCount .eq. this%BudgetReader%GetFlowArraySize()) then
                          call this%BudgetReader%FillRecordDataBuffer(header,         &
                            this%FlowsRightFace, header%ArrayItemCount, spaceAssigned,&
                            status)
                      end if
                  end if
                 
             case('FLOW FRONT FACE')
                  ! Read flows across the front face for structured grids.
                  ! Method should always be 0 or 1, but check anyway to be sure.
                  if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
                      if(header%ArrayItemCount .eq. this%BudgetReader%GetFlowArraySize()) then
                          call this%BudgetReader%FillRecordDataBuffer(header,         &
                            this%FlowsFrontFace, header%ArrayItemCount, spaceAssigned,&
                            status)
                      end if
                  end if
                 
             case('FLOW LOWER FACE')
                  ! Read flows across the lower face for structured grids.
                  ! Method should always be 0 or 1, but check anyway to be sure.
                  if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
                      if(header%ArrayItemCount .eq. this%BudgetReader%GetFlowArraySize()) then
                          call this%BudgetReader%FillRecordDataBuffer(header,         &
                            this%FlowsLowerFace, header%ArrayItemCount, spaceAssigned,&
                            status)
                      end if
                  end if
             
              case default
                  ! Now handle any other records in the budget file.
                   if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
                      if(header%ArrayItemCount .eq. cellCount) then
                          call this%BudgetReader%FillRecordDataBuffer(header,         &
                            this%ArrayBufferDbl, cellCount, spaceAssigned, status)
                          if(cellCount .eq. spaceAssigned) then
                              call this%CheckForDefaultIface(header%TextLabel, iface)
                              if(iface .gt. 0) then
                                  do m = 1, spaceAssigned
                                      boundaryFlowsOffset = 6 * (m - 1)
                                      this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
                                        this%BoundaryFlows(boundaryFlowsOffset + iface) + &
                                        this%ArrayBufferDbl(m)
                                  end do
                              else
                                  do m = 1, spaceAssigned
                                      if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
                                          this%SourceFlows(m) = this%SourceFlows(m) +     &
                                            this%ArrayBufferDbl(m)
                                      else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
                                          this%SinkFlows(m) = this%SinkFlows(m) +         &
                                            this%ArrayBufferDbl(m)
                                      end if
                                  end do
                              end if
                          end if
                      end if
                   else if(header%Method .eq. 3) then
                      call this%BudgetReader%FillRecordDataBuffer(header,             &
                        this%ArrayBufferDbl, this%ArrayBufferInt,                     &
                        header%ArrayItemCount, spaceAssigned, status)
                      if(header%ArrayItemCount .eq. spaceAssigned) then
                          call this%CheckForDefaultIface(header%TextLabel, iface)
                          if(iface .gt. 0) then
                              do m = 1, spaceAssigned
                                  if(this%Grid%GridType .ne. 2) then
                                    ! structured
                                    layer = this%ArrayBufferInt(m)
                                    cellNumber = (layer - 1) * spaceAssigned + m
                                  else
                                    ! mfusg unstructured
                                    cellNumber = this%ArrayBufferInt(m)
                                  end if
                                  boundaryFlowsOffset = 6 * (cellNumber - 1)
                                  this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
                                    this%BoundaryFlows(boundaryFlowsOffset + iface) + &
                                    this%ArrayBufferDbl(m)
                              end do
                          else            
                              do m = 1, spaceAssigned
                                  cellNumber = this%ArrayBufferInt(m)
                                  if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
                                      this%SourceFlows(cellNumber) =                  &
                                        this%SourceFlows(cellNumber) +                &
                                        this%ArrayBufferDbl(m)
                                  else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
                                      this%SinkFlows(cellNumber) =                    &
                                        this%SinkFlows(cellNumber) +                  &
                                        this%ArrayBufferDbl(m)
                                  end if
                              end do
                          end if
                      end if
                   else if(header%Method .eq. 4) then
                      call this%BudgetReader%FillRecordDataBuffer(header,             &
                        this%ArrayBufferDbl, header%ArrayItemCount, spaceAssigned,    &
                        status)
                      if(header%ArrayItemCount .eq. spaceAssigned) then
                          call this%CheckForDefaultIface(header%TextLabel, iface)
                          if(iface .gt. 0) then
                              do m = 1, spaceAssigned
                                  boundaryFlowsOffset = 6 * (m - 1)
                                  this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
                                    this%BoundaryFlows(boundaryFlowsOffset + iface) + &
                                    this%ArrayBufferDbl(m)
                              end do
                          else            
                              do m = 1, spaceAssigned
                                  if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
                                      this%SourceFlows(m) = this%SourceFlows(m) +     &
                                        this%ArrayBufferDbl(m)
                                  else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
                                      this%SinkFlows(m) = this%SinkFlows(m) +         &
                                        this%ArrayBufferDbl(m)
                                  end if
                              end do
                          end if
                      end if
                  else if(header%Method .eq. 2) then
                      call this%BudgetReader%FillRecordDataBuffer(header,             &
                        this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
                        status)
                      if(spaceAssigned .gt. 0) then
                          call this%CheckForDefaultIface(header%TextLabel, iface)
                          if(iface .gt. 0) then
                              do m = 1, spaceAssigned
                                  cellNumber = this%ListItemBuffer(m)%CellNumber
                                  boundaryFlowsOffset = 6 * (cellNumber - 1)
                                  this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
                                    this%BoundaryFlows(boundaryFlowsOffset + iface) + &
                                    this%ListItemBuffer(m)%BudgetValue
                              end do
                          else            
                              do m = 1, spaceAssigned
                                  cellNumber = this%ListItemBuffer(m)%CellNumber
                                  if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
                                      this%SourceFlows(cellNumber) =                  &
                                        this%SourceFlows(cellNumber) +                &
                                        this%ListItemBuffer(m)%BudgetValue
                                  else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
                                      this%SinkFlows(cellNumber) =                    &
                                        this%SinkFlows(cellNumber) +                  &
                                        this%ListItemBuffer(m)%BudgetValue
                                  end if
                              end do
                          end if
                      end if
                  else if((header%Method .eq. 5) .or. (header%Method .eq. 6)) then
                      call this%BudgetReader%FillRecordDataBuffer(header,             &
                        this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
                        status)
                      if(spaceAssigned .gt. 0) then
                          do m = 1, spaceAssigned
                              call this%CheckForDefaultIface(header%TextLabel, iface)
                              index = header%FindAuxiliaryNameIndex('IFACE')
                              if(index .gt. 0) then
                                  iface = int(this%ListItemBuffer(m)%AuxiliaryValues(index))
                              end if
                              
                              cellNumber = this%ListItemBuffer(m)%CellNumber
                              if(iface .gt. 0) then
                                  boundaryFlowsOffset = 6 * (cellNumber - 1)
                                  this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
                                    this%BoundaryFlows(boundaryFlowsOffset + iface) + &
                                    this%ListItemBuffer(m)%BudgetValue
                              else
                                  if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
                                      this%SourceFlows(cellNumber) =                  &
                                        this%SourceFlows(cellNumber) +                &
                                        this%ListItemBuffer(m)%BudgetValue
                                  else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
                                      this%SinkFlows(cellNumber) =                    &
                                        this%SinkFlows(cellNumber) +                  &
                                        this%ListItemBuffer(m)%BudgetValue
                                  end if
                              end if
                          end do
                      end if
                  end if
             
             end select
             
        end do
    
        this%CurrentStressPeriod = stressPeriod
        this%CurrentTimeStep = timeStep
    
    end subroutine pr_LoadTimeStep


    subroutine pr_ClearTimeStepBudgetData(this)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    !
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(FlowModelDataType) :: this
    integer :: cellCount, n, arraySize
    !---------------------------------------------------------------------------------------------------------------
      
        this%CurrentStressPeriod = 0
        this%CurrentTimeStep = 0
        
        if(allocated(this%SinkFlows)) then
            cellCount = this%Grid%CellCount
            do n = 1, cellCount
                this%IBoundTS(n) = this%IBound(n)
                this%Heads(n) = 0.0d0
                this%SourceFlows(n) = 0.0d0
                this%SinkFlows(n) = 0.0d0
                this%StorageFlows(n) = 0.0d0
                this%SubFaceFlowsComputed(n) = .false.
            end do
            
            arraySize = cellCount * 6
            do n = 1, arraySize
                this%BoundaryFlows(n) = 0.0d0
            end do
            
            arraySize = cellCount * 4
            do n = 1, arraySize
                this%SubFaceFlows(n) = 0.0d0
            end do
        
            arraySize = this%BudgetReader%GetFlowArraySize()
            if(this%Grid%GridType .eq. 1) then
                do n = 1, arraySize
                this%FlowsRightFace(n) = 0.0d0
                this%FlowsFrontFace(n) = 0.0d0
                this%FlowsLowerFace(n) = 0.0d0
                end do
            else if(this%Grid%GridType .eq. 2) then
                do n = 1, arraySize
                    this%FlowsJA(n) = 0.0d0
                end do
            end if
           
        end if
 

    end subroutine pr_ClearTimeStepBudgetData


    function pr_GetCurrentStressPeriod(this) result(stressPeriod)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(FlowModelDataType) :: this
    integer :: stressPeriod
    !---------------------------------------------------------------------------------------------------------------
     

        stressPeriod = this%CurrentStressPeriod
     

    end function pr_GetCurrentStressPeriod


    function pr_GetCurrentTimeStep(this) result(timeStep)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(FlowModelDataType) :: this
    integer :: timeStep
    !---------------------------------------------------------------------------------------------------------------
    

        timeStep = this%CurrentTimeStep
      

    end function pr_GetCurrentTimeStep


    subroutine pr_SetIbound(this, ibound, arraySize)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(FlowModelDataType) :: this
    integer,intent(in) :: arraySize
    integer :: n
    integer,dimension(arraySize),intent(in),target :: ibound
    !---------------------------------------------------------------------------------------------------------------
      

        if(arraySize .ne. this%Grid%CellCount) then
            write(*,*) "FlowModelDataType: The IBound array size does not match the cell count for the grid. stop"
            stop
        end if
        
        this%IBound => ibound
        ! Initialize the IBoundTS array to the same values as IBound whenever the IBound array is set.
        ! The values of IBoundTS will be updated for dry cells every time that data for a time step is loaded.
        do n = 1, arraySize
            this%IBoundTS(n) = this%IBound(n)
        end do
  

    end subroutine pr_SetIbound


    subroutine pr_SetZones(this, zones, arraySize)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(FlowModelDataType) :: this
    integer,intent(in) :: arraySize
    integer,dimension(arraySize),intent(in),target :: zones
    !---------------------------------------------------------------------------------------------------------------
      

        if(arraySize .ne. this%Grid%CellCount) then
            write(*,*) "FlowModelDataType: The Zones array size does not match the cell count for the grid. stop"
            stop
        end if
        
        this%Zones => zones
    

    end subroutine pr_SetZones


    subroutine pr_SetPorosity(this, porosity, arraySize)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(FlowModelDataType) :: this
    integer,intent(in) :: arraySize
    doubleprecision,dimension(arraySize),intent(in),target :: porosity
    !---------------------------------------------------------------------------------------------------------------


        if(arraySize .ne. this%Grid%CellCount) then
            write(*,*) "FlowModelDataType: The Porosity array size does not match the cell count for the grid. stop"
            stop
        end if
        
        this%Porosity => porosity


    end subroutine pr_SetPorosity


    subroutine pr_SetRetardation(this, retardation, arraySize)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(FlowModelDataType) :: this
    integer,intent(in) :: arraySize
    doubleprecision,dimension(arraySize),intent(in),target :: retardation
    !---------------------------------------------------------------------------------------------------------------
  

        if(arraySize .ne. this%Grid%CellCount) then
            write(*,*) "FlowModelDataType: The Retardation array size does not match the cell count for the grid. stop"
            stop
        end if
        
        this%Retardation => retardation
 

    end subroutine pr_SetRetardation

  
    subroutine pr_SetDefaultIface(this, defaultIfaceLabels, defaultIfaceValues, arraySize)
    !***************************************************************************************************************
    ! Description goes here
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    use UtilMiscModule,only : utrimall
    implicit none
    class(FlowModelDataType) :: this
    integer,intent(in) :: arraySize
    integer,dimension(arraySize),intent(in) :: defaultIfaceValues
    character(len=16),dimension(arraySize),intent(in) :: defaultIfaceLabels
    integer :: n, firstNonBlank, lastNonBlank, trimmedLength
    character(len=16) :: label
    !---------------------------------------------------------------------------------------------------------------
      

        this%DefaultIfaceCount = 0
        do n = 1, 20
            this%DefaultIfaceValues(n) = 0
            this%DefaultIfaceLabels(n) = '                '
        end do
        
        do n = 1, arraySize
            this%DefaultIfaceValues(n) = defaultIfaceValues(n)
            label = defaultIfaceLabels(n)
            call utrimall(label)
            this%DefaultIfaceLabels(n) = label
        end do
        this%DefaultIfaceCount = arraySize
   

    end subroutine pr_SetDefaultIface


    subroutine pr_CheckForDefaultIface(this, textLabel, iface)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    !
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    use UtilMiscModule,only : utrimall
    implicit none
    class(FlowModelDataType) :: this
    character*(*), intent(in) :: textLabel
    integer,intent(inout) :: iface
    integer :: n
    character(len=16) :: label
    !---------------------------------------------------------------------------------------------------------------
      
        iface = 0
        label = textLabel
        call utrimall(label)
        do n = 1, this%DefaultIfaceCount
            if(label .eq. this%DefaultIfaceLabels(n)) then
                iface = this%DefaultIfaceValues(n)
                return
            end if
        end do
      
    end subroutine pr_CheckForDefaultIface


    !subroutine pr_SetLayerTypes(this, layerTypes, arraySize)
    !!***************************************************************************************************************
    !!
    !!***************************************************************************************************************
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !  implicit none
    !  class(FlowModelDataType) :: this
    !  integer,intent(in) :: arraySize
    !  integer,dimension(arraySize),intent(in),target :: layerTypes
    !!---------------------------------------------------------------------------------------------------------------
    !
    !  
    !  if(arraySize .ne. this%Grid%LayerCount) then
    !      write(*,*) "FlowModelDataType: The LayerTypes array size does not match the layer count for the grid. stop"
    !      stop
    !  end if
    !  
    !  this%LayerTypes => layerTypes
    !
    !
    !end subroutine pr_SetLayerTypes


end module FlowModelDataModule
