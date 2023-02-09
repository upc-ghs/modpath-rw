module FlowModelDataModule
  use BudgetReaderModule,only : BudgetReaderType
  use HeadReaderModule,only : HeadReaderType
  use BudgetListItemModule,only : BudgetListItemType
  use ModflowRectangularGridModule,only : ModflowRectangularGridType
  use BudgetRecordHeaderModule,only : BudgetRecordHeaderType
  use UtilMiscModule,only : TrimAll
  use UTL8MODULE,only : ustop
  implicit none
  !--------------------------------------------------------------------------------------

  ! Set default access status to private
  private

    type,public :: FlowModelDataType
      logical :: Initialized = .false.
      logical :: SteadyState = .true.
      integer :: DefaultIfaceCount
      doubleprecision :: HDry = 0d0
      doubleprecision :: HNoFlow = 0d0
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

      ! DEV
      procedure :: LoadFlowTimeSeries => pr_LoadFlowTimeSeries
      procedure :: ValidateAuxVarNames => pr_ValidateAuxVarNames
      procedure :: LoadFlowAndAuxTimeseries => pr_LoadFlowAndAuxTimeseries

    end type


contains


    subroutine pr_Initialize(this,headReader, budgetReader, grid, hNoFlow, hDry)
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    ! Specifications
    !---------------------------------------------------------------------------
    implicit none
    class(FlowModelDataType) :: this
    type(BudgetReaderType),intent(inout),target :: budgetReader
    type(HeadReaderType),intent(inout),target :: headReader
    class(ModflowRectangularGridType),intent(inout),pointer :: grid
    integer :: cellCount,gridType
    integer :: n, flowArraySize
    doubleprecision :: hNoFlow, hDry
    !---------------------------------------------------------------------------

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
          ! RWPT: In case budget and/or head were not saved 
          ! for all times, there is an error triggered 
          ! downstream because of the early exits at this point.
          ! Specifically, the headReader%CellCount is not properly defined.
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
    !---------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    ! Specifications
    !---------------------------------------------------------------------------
    implicit none
    class(FlowModelDataType) :: this
    !---------------------------------------------------------------------------
       
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
    !------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------
    ! Specifications
    !------------------------------------------------------------------------
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
    !------------------------------------------------------------------------
      
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
      if(this%Grid%GridType .gt. 2) then ! MODFLOW-6 DIS, DISV
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
            call this%BudgetReader%FillRecordDataBuffer(header,       &
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
            call this%BudgetReader%FillRecordDataBuffer(header,       &
              this%ListItemBuffer, listItemBufferSize, spaceAssigned, &
              status)
            if(spaceAssigned .gt. 0) then
              do m = 1, spaceAssigned
                cellNumber = this%ListItemBuffer(m)%CellNumber
                if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
                  this%SourceFlows(cellNumber) =                  &
                      this%SourceFlows(cellNumber) +              &
                      this%ListItemBuffer(m)%BudgetValue
                else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
                  this%SinkFlows(cellNumber) =                    &
                      this%SinkFlows(cellNumber) +                &
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

        case('DATA-SPDIS', 'DATA-SAT')
          ! Skip
          ! The “DATA” prefix on the text identifier can
          ! be used by post-processors to recognize that the record
          ! does not contain a cell flow budget term.
          cycle
          
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


    subroutine pr_SetIBound(this, ibound, arraySize)
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
  

    end subroutine pr_SetIBound


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


    subroutine pr_LoadFlowAndAuxTimeseries(this, sourcePkgName, auxVarNames,& 
                                    isMF6, initialTime, finalTime, tdisData,&
                        flowTimeseries, auxTimeseries, timeIntervals, times,&
                              cellNumbers, outUnit, iFaceOption, iFaceCells )
    use TimeDiscretizationDataModule,only : TimeDiscretizationDataType
    !------------------------------------------------------------------------
    ! Given range of times, extract timeseries for flow and aux vars related 
    ! to package name/budget header 
    !------------------------------------------------------------------------
    ! Specifications
    !------------------------------------------------------------------------
    implicit none
    ! input
    class(FlowModelDataType) :: this
    character(len=16), intent(in) :: sourcePkgName
    character(len=16), dimension(:), intent(in) :: auxVarNames
    logical, intent(in) :: isMF6
    doubleprecision, intent(in) :: initialTime, finalTime
    type( TimeDiscretizationDataType ), intent(in) :: tdisData
    integer, optional, intent(in) :: outUnit
    logical, optional, intent(in) :: iFaceOption
    ! out
    doubleprecision, allocatable, dimension(:,:)  , intent(inout) :: flowTimeseries ! nt x ncells
    doubleprecision, allocatable, dimension(:,:,:), intent(inout) :: auxTimeseries  ! nt x ncells x nauxvars
    doubleprecision, allocatable, dimension(:)    , intent(inout) :: timeIntervals  ! nt
    doubleprecision, allocatable, dimension(:)    , intent(inout) :: times          ! nt + 1
    integer        , allocatable, dimension(:)    , intent(inout) :: cellNumbers    ! nCells
    integer, allocatable, dimension(:), optional  , intent(inout) :: iFaceCells     ! nCells
    ! local
    logical :: lookForIFace = .false.
    integer :: ifaceindex
    integer :: n, m, naux
    integer :: stressPeriod, timeStep
    integer :: firstRecord,lastRecord
    integer :: firstNonBlank,lastNonBlank,trimmedLength
    integer :: firstNonBlankIn,lastNonBlankIn,trimmedLengthIn
    integer :: firstNonBlankLoc,lastNonBlankLoc,trimmedLengthLoc
    integer :: spaceAssigned, status, cellCount, iface, index, auxindex, cellindex
    integer :: listItemBufferSize, cellNumber
    type(BudgetRecordHeaderType) :: header
    character(len=16)  :: textLabel
    character(len=16)  :: textNameLabel
    character(len=132) :: message
    integer :: nCells, newcounter
    integer :: kinitial, kfinal, ktime, kcounter
    integer :: nTimes, nTimeIntervals, nAuxVars
    integer :: spInit, tsInit, spEnd, tsEnd, nStressPeriods, nsp 
    integer, allocatable, dimension(:) :: tempCellNumbers
    integer, allocatable, dimension(:) :: spCellNumbers
    integer, allocatable, dimension(:) :: tempSPIFaces
    integer, allocatable, dimension(:) :: spIFaces
    ! For identifying pkgs with aux vars from modflow != mf6
    logical :: foundTheSource = .false.
    integer :: nb, nbindex
    integer :: nbmax = 5
    character(len=16)  :: anamebud(5)
    DATA anamebud(1) /'           WELLS'/ ! WEL
    DATA anamebud(2) /'    DRAINS (DRT)'/ ! DRT
    DATA anamebud(3) /'          DRAINS'/ ! DRN
    DATA anamebud(4) /'   RIVER LEAKAGE'/ ! RIV
    DATA anamebud(5) /' HEAD DEP BOUNDS'/ ! GHB
    character(len=16)  :: anameid(5)
    DATA anameid(1)  /'             WEL'/ ! WEL
    DATA anameid(2)  /'             DRT'/ ! DRT
    DATA anameid(3)  /'             DRN'/ ! DRN
    DATA anameid(4)  /'             RIV'/ ! RIV
    DATA anameid(5)  /'             GHB'/ ! GHB
    !------------------------------------------------------------------------

      ! Check the iface flag
      if ( present( iFaceOption) ) then 
        lookForIface = .false.
        ifaceindex   = 0
        lookForIFace = iFaceOption
      end if 

      ! Trim input pkg name
      call TrimAll(sourcePkgName, firstNonBlankIn, lastNonBlankIn, trimmedLengthIn)
      listItemBufferSize = size(this%ListItemBuffer)

      ! Given initial and final times, 
      ! compute the initial and final time step indexes
      ! Note 1: TotalTimes starts from dt, not zero.
      ! Note 2: FindContainingTimeStep returns the index 
      ! correspoding to the upper limit of the time interval 
      ! in TotalTimes. Meaning, if on a TotalTimes vector 
      ! [dt,2dt,...] the time 1.5dt is requested, then 
      ! the function will return the index 2, corresponding to 
      ! TotalTimes=2dt. 
      kinitial = tdisData%FindContainingTimeStep(initialTime)
      kfinal   = tdisData%FindContainingTimeStep(finalTime)
    
      ! There are cases in which, depending on the stoptimeoption, 
      ! finalTime may have the value 1.0d+30. Something really 
      ! big to track particles until all of them get to a stop
      ! conditio for steady state models. For such value, 
      ! FindContainingTimeStep will return, assuming that this large
      ! number is higher than the length of the modflow simulation stoptime.
      ! In such case, it is enforced that kfinal adopt the value of 
      ! tdisData%CumulativeTimeStepCount, the highest possible value.
      !if ( (kfinal .eq. 0) .and. (this%SteadyState)  ) then
      if ( (kfinal .eq. 0) ) then
       if ( present(outUnit) ) then       
        write(outUnit,'(a)') 'FlowModelData: LoadFlowAndAuxTimeseries: kfinal is assumed to be CumulativeTimeStepCount'
        write(outUnit,'(a,e15.7)') 'FlowModelData: LoadFlowAndAuxTimeseries: final time is ', finalTime
       end if
       kfinal = tdisData%CumulativeTimeStepCount
      end if

      ! The number of intervals
      nTimeIntervals = kfinal - kinitial + 1
      nTimes = nTimeIntervals + 1 
      ! Something wrong with times 
      if ( nTimeIntervals .lt. 1 ) then 
         write(message,'(A)') 'Error: the number of times is .lt. 1. Check definition of reference and stoptimes. Stop.'
         message = trim(message)
         call ustop(message)
      end if  

      ! times: includes intial and final times ( reference, stoptime )
      if ( allocated( times ) ) deallocate( times ) 
      allocate( times(nTimes) )
      if ( allocated( timeIntervals ) ) deallocate( timeIntervals ) 
      allocate( timeIntervals(nTimeIntervals) )  

      ! Fill times and time intervals
      do ktime=1,nTimeIntervals
        if ( ktime .eq. 1 ) times(1) = initialTime
        if ( ktime .lt. nTimeIntervals ) times(ktime+1) = tdisData%TotalTimes(ktime+kinitial-1)
        if ( ktime .eq. nTimeIntervals ) times(nTimes)  = finalTime 
        timeIntervals(ktime) = times(ktime+1) - times(ktime)
      end do 

      nAuxVars = size(auxVarNames)
      if ( nAuxVars .lt. 1 ) then 
         write(message,'(A)') 'Error: number of aux variables for timeseries should be at least 1. Stop.'
         message = trim(message)
         call ustop(message)
      end if  

      ! In order to extract the pkg cells, it needs to 
      ! run over stress periods

      ! Get the initial and final stress
      call tdisData%GetPeriodAndStep(kinitial, spInit, tsInit)
      call tdisData%GetPeriodAndStep(kfinal  , spEnd , tsEnd )

      nStressPeriods = spEnd - spInit + 1
      timeStep = 1

      ! Loop over range of stress periods
      do nsp=1, nStressPeriods

        ! Determine record range for stressPeriod and timeStep
        call this%BudgetReader%GetRecordHeaderRange(nsp, timeStep, firstRecord, lastRecord)

        if(firstRecord .eq. 0) then
          write(message,'(A,I5,A,I5,A)') ' Error loading Time Step ', timeStep, ' Period ', nsp, '.'
          message = trim(message)
          write(*,'(A)') message
          call ustop('Missing budget information. Budget file must have output for every time step. Stop.')
        end if

        ! Loop through record headers
        do n = firstRecord, lastRecord
          header    = this%BudgetReader%GetRecordHeader(n)
          if ( ( header%Method .eq. 5 ) .or. ( header%Method .eq. 6 ) ) then

            ! Is the requested pkg ?
            foundTheSource = .false.

            ! MF6
            if ( isMF6 ) then 
              textLabel = header%TXT2ID2
              call TrimAll(textLabel, firstNonBlank, lastNonBlank, trimmedLength)
              if (&
                textLabel(firstNonBlank:lastNonBlank) .eq. & 
                sourcePkgName(firstNonBlankIn:lastNonBlankIn) ) then
                foundTheSource = .true.
              end if 
            ! OTHER MODFLOW
            else
              textLabel = header%TextLabel
              call TrimAll(textLabel, firstNonBlank, lastNonBlank, trimmedLength)

              ! Needs to verify relation
              ! Find equivalence
              nbindex = 0
              do nb=1,nbmax
                textNameLabel = anamebud(nb) 
                call TrimAll(textNameLabel, firstNonBlankLoc, lastNonBlankLoc, trimmedLengthLoc)
                if (&
                  textNameLabel(firstNonBlankLoc:lastNonBlankLoc) .eq. & 
                  textLabel(firstNonBlank:lastNonBlank) ) then
                  ! Found, continue
                  nbindex = nb
                  exit
                end if   
              end do

              ! Not found in the list of known budgets supporting
              ! aux variables, try next budget header. It might be useful 
              ! to report the header text label for validation.
              if ( nbindex .eq. 0 ) then
                exit
              end if

              ! Compare the id/ftype (e.g. WEL) against the given src name,
              ! and if not, give it another chance by comparing against the 
              ! budget label itself (e.g. WELLS)
              textNameLabel = anameid(nbindex) 
              call TrimAll(textNameLabel, firstNonBlankLoc, lastNonBlankLoc, trimmedLengthLoc)
              if (&
                textNameLabel(firstNonBlankLoc:lastNonBlankLoc) .eq. & 
                sourcePkgName(firstNonBlankIn:lastNonBlankIn) ) then
                foundTheSource = .true.
              else
                textNameLabel = anamebud(nbindex) 
                call TrimAll(textNameLabel, firstNonBlankLoc, lastNonBlankLoc, trimmedLengthLoc)
                if (&
                  textNameLabel(firstNonBlankLoc:lastNonBlankLoc) .eq. & 
                  sourcePkgName(firstNonBlankIn:lastNonBlankIn) ) then
                  foundTheSource = .true.
                end if
              end if
            end if ! isMF6

            if ( foundTheSource ) then 
              ! Check cells
              call this%BudgetReader%FillRecordDataBuffer(header,       &
                this%ListItemBuffer, listItemBufferSize, spaceAssigned, &
                status)
              if(spaceAssigned .gt. 0) then
                      
                ! If allocated with different size, reallocate 
                ! else restart indexes
                if ( allocated(spCellNumbers) ) then 
                  if ( size(spCellNumbers) .ne. spaceAssigned ) then 
                    deallocate( spCellNumbers )
                    allocate(spCellNumbers(spaceAssigned))
                    if( lookForIFace ) then 
                      deallocate( spIFaces )
                      allocate(spIFaces(spaceAssigned))
                      spIFaces(:) = 0
                    end if 
                  else
                    spCellNumbers(:) = 0
                    if( lookForIFace ) then 
                      spIFaces(:) = 0
                    end if 
                  end if
                else
                  allocate(spCellNumbers(spaceAssigned))
                  if( lookForIFace ) then 
                    allocate(spIFaces(spaceAssigned))
                    spIFaces(:) = 0
                  end if 
                end if

                ! Assign to stress period cell numbers
                do m = 1, spaceAssigned
                  cellNumber = this%ListItemBuffer(m)%CellNumber
                  spCellNumbers(m) = cellNumber
                end do
                
                ! Again, it is considered that the existence of IFACE 
                ! has been already established, as every other aux var, 
                ! regardless of which value it may have
                if ( lookForIFace ) then 
                  ifaceindex = header%FindAuxiliaryNameIndex('IFACE')
                  if ( ifaceindex .gt. 0 ) then 
                   do m = 1, spaceAssigned
                    spIFaces(m) = int(this%ListItemBuffer(m)%AuxiliaryValues(ifaceindex))
                   end do
                  end if
                end if 

                ! Assign cellNumbers
                if ( .not. allocated( cellNumbers ) ) then
                  ! First initialization
                  allocate( cellNumbers(spaceAssigned) )
                  cellNumbers(:) = spCellNumbers(:)
                  if ( lookForIFace ) then 
                    allocate( iFaceCells(spaceAssigned) )
                    iFaceCells(:) = spIFaces(:)
                  end if 
                  ! Break the records loop and continue to next stress period
                  exit
                else
                  ! If allocated, verify if any new cell
                  newcounter = 0
                  do m =1, spaceAssigned
                    cellindex = findloc( cellNumbers, spCellNumbers(m), 1 ) 
                    if ( cellindex .eq. 0 ) newcounter = newcounter + 1 ! is new cell
                  end do 
                  ! If any new, add it to cellNumbers
                  if ( newcounter .gt. 0 ) then 
                    if ( allocated( tempCellNumbers ) ) deallocate( tempCellNumbers ) 
                    allocate( tempCellNumbers(size(cellNumbers)+newcounter) )
                    tempCellNumbers(1:size(cellNumbers)) = cellNumbers(:) ! save the old
                    if ( lookForIFace ) then 
                      if ( allocated( tempSPIFaces ) ) deallocate( tempSPIFaces ) 
                      allocate( tempSPIFaces(size(iFaceCells)+newcounter) )
                      tempSPIFaces(1:size(iFaceCells)) = iFaceCells(:) ! save the old
                    end if 
                    newcounter = 0
                    do m =1, spaceAssigned
                      cellindex = findloc( cellNumbers, spCellNumbers(m), 1 )
                      if ( cellindex .eq. 0 ) then 
                        newcounter = newcounter + 1 ! is new cell
                        tempCellNumbers(size(cellNumbers)+newcounter) = spCellNumbers(m)
                        if( lookForIFace ) then 
                          tempSPIFaces(size(iFaceCells)+newcounter) = spIFaces(m)
                        end if 
                      end if
                    end do
                    call move_alloc( tempCellNumbers, cellNumbers )
                    if ( lookForIFace ) then 
                      call move_alloc( tempSPIFaces, iFaceCells )
                    end if 
                  end if
                  ! Break the records loop and continue to next stress period
                  exit
                end if !Assign cellNumbers

              end if !if(spaceAssigned .gt. 0)

            end if ! foundTheSource

          end if ! ( header%Method .eq. 5 ) .or. ( header%Method .eq. 6 )

        end do !n = firstRecord, lastRecord

      end do !nsp=1, nStressPeriods

      
      print *, 'VERIFY IFACEs'
      print *, cellNumbers
      print *, iFaceCells
      call exit(0)


     
      ! No cells found, something wrong 
      if ( .not. allocated( cellNumbers ) ) then 
         write(message,'(A,A,A)') 'Error: no cells were found for source package ', trim(adjustl(sourcePkgName)), '. Stop.'
         message = trim(message)
         call ustop(message)
      end if  
      nCells = size(cellNumbers)

      ! Allocate arrays for storing timeseries
      if( allocated( flowTimeseries ) ) deallocate( flowTimeseries ) 
      allocate( flowTimeseries( nTimeIntervals, nCells ) )
      flowTimeseries(:,:) = 0d0
      if ( allocated( auxTimeseries ) ) deallocate( auxTimeseries ) 
      allocate( auxTimeseries( nTimeIntervals, nCells, nAuxVars ) )
      auxTimeseries(:,:,:) = 0d0

      ! Supposedly, previous to run this function 
      ! aux variables were already validated with 
      ! call this%ValidateAuxVarNames

      ! Use the determined steps (kinitial,kfinal) to build the timeseries
      kcounter = 0
      do ktime=kinitial,kfinal

        ! Get the stress period and time step from the cummulative time steps
        call tdisData%GetPeriodAndStep(ktime, stressPeriod, timeStep)
        kcounter = kcounter + 1 

        ! Determine record range for stressPeriod and timeStep
        call this%BudgetReader%GetRecordHeaderRange(stressPeriod, timeStep, firstRecord, lastRecord)
        if(firstRecord .eq. 0) then
          write(message,'(A,I5,A,I5,A)') ' Error loading Time Step ', timeStep, ' Period ', stressPeriod, '.'
          message = trim(message)
          write(*,'(A)') message
          call ustop('Missing budget information. Budget file must have output for every time step. Stop.')
        end if


        ! Loop through record headers
        do n = firstRecord, lastRecord
          header    = this%BudgetReader%GetRecordHeader(n)
          ! Only methods 5,6 support aux variables
          if ( ( header%Method .eq. 5 ) .or. ( header%Method .eq. 6 ) ) then

            ! Is the requested pkg ?
            foundTheSource = .false.

            ! MF6
            if ( isMF6 ) then 
              textLabel = header%TXT2ID2
              call TrimAll(textLabel, firstNonBlank, lastNonBlank, trimmedLength)
              if (&
                textLabel(firstNonBlank:lastNonBlank) .eq. & 
                sourcePkgName(firstNonBlankIn:lastNonBlankIn) ) then
                foundTheSource = .true.
              end if 
            ! OTHER MODFLOW
            else
              textLabel = header%TextLabel
              call TrimAll(textLabel, firstNonBlank, lastNonBlank, trimmedLength)

              ! Needs to verify relation
              ! Find equivalence
              nbindex = 0
              do nb=1,nbmax
                textNameLabel = anamebud(nb) 
                call TrimAll(textNameLabel, firstNonBlankLoc, lastNonBlankLoc, trimmedLengthLoc)
                if (&
                  textNameLabel(firstNonBlankLoc:lastNonBlankLoc) .eq. & 
                  textLabel(firstNonBlank:lastNonBlank) ) then
                  ! Found, continue
                  nbindex = nb
                  exit
                end if   
              end do

              ! Not found in the list of known budgets supporting
              ! aux variables, try next budget header. It might be useful 
              ! to report the header text label for validation.
              if ( nbindex .eq. 0 ) then
                exit
              end if

              ! Compare the id/ftype (e.g. WEL) against the given src name,
              ! and if not, give it another chance by comparing against the 
              ! budget label itself (e.g. WELLS)
              textNameLabel = anameid(nbindex) 
              call TrimAll(textNameLabel, firstNonBlankLoc, lastNonBlankLoc, trimmedLengthLoc)
              if (&
                textNameLabel(firstNonBlankLoc:lastNonBlankLoc) .eq. & 
                sourcePkgName(firstNonBlankIn:lastNonBlankIn) ) then
                foundTheSource = .true.
              else
                textNameLabel = anamebud(nbindex) 
                call TrimAll(textNameLabel, firstNonBlankLoc, lastNonBlankLoc, trimmedLengthLoc)
                if (&
                  textNameLabel(firstNonBlankLoc:lastNonBlankLoc) .eq. & 
                  sourcePkgName(firstNonBlankIn:lastNonBlankIn) ) then
                  foundTheSource = .true.
                end if
              end if
            end if ! isMF6

            if ( foundTheSource ) then 
              ! Found the pkg
              call this%BudgetReader%FillRecordDataBuffer(header,             &
                this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
                status)
              if(spaceAssigned .gt. 0) then
                do m = 1, spaceAssigned
                  cellNumber = this%ListItemBuffer(m)%CellNumber

                  ! Determine the index of cellNumber in the list of cells 
                  ! requested for timeseseries
                  cellindex = findloc( cellNumbers, cellNumber, 1 ) 
                  if ( cellindex .eq. 0 ) cycle ! Not found, but it should not be the case

                  !call this%CheckForDefaultIface(header%TextLabel, iface)
                  !index = header%FindAuxiliaryNameIndex('IFACE')
                  !if(index .gt. 0) then
                  !  iface = int(this%ListItemBuffer(m)%AuxiliaryValues(index))
                  !end if
                  !if(iface .gt. 0) ifaces(cellindex) = iface

                  ! Load into flow rates timeseries only if positive, 
                  ! otherwise leave as zero. Notice that the same 
                  ! applies to concentration. However, because aux var
                  ! could be something else (?), it will save whatever 
                  ! it finds
                  if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
                    flowTimeseries( kcounter, cellindex )  = this%ListItemBuffer(m)%BudgetValue
                  end if
                  
                  ! Load aux vars
                  do naux=1, nAuxVars
                    auxindex = header%FindAuxiliaryNameIndex(auxVarNames(naux))
                    if(auxindex .gt. 0) then
                      auxTimeseries( kcounter, cellindex, naux ) = this%ListItemBuffer(m)%AuxiliaryValues(auxindex)
                    end if
                  end do

                end do ! spaceAssigned
              end if
              
              ! Break the records loop and continue to next ktime
              exit

            end if ! foundTheSource

          end if ! ( header%Method .eq. 5 ) .or. ( header%Method .eq. 6 ) 

        end do ! n = firstRecord, lastRecord

      end do ! ktime=kinitial,kfinal


      if( allocated(tempCellNumbers))deallocate(tempCellNumbers)
      if( allocated(spCellNumbers)  )deallocate(spCellNumbers)


    end subroutine pr_LoadFlowAndAuxTimeseries



    function pr_ValidateAuxVarNames( this, sourcePkgName, auxVarNames, isMF6,& 
                                              iFaceOption ) result ( isValid )
    !------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------
    ! Specifications
    !------------------------------------------------------------------------
    implicit none
    ! input
    class(FlowModelDataType) :: this
    character(len=16), intent(in) :: sourcePkgName
    character(len=16), dimension(:), intent(in) :: auxVarNames
    logical, intent(in) :: isMF6
    logical, optional, intent(in) :: iFaceOption
    ! output
    logical :: isValid
    ! local
    integer :: timeStep = 1 ! to look for aux vars, use as ref...
    integer :: stressPeriod = 1 ! the first stperiod, first tstep
    integer :: n, m, naux, nx, nval, auxindex
    integer :: cellNumber, status
    integer :: firstRecord,lastRecord
    integer :: listItemBufferSize, spaceAssigned
    integer :: firstNonBlank,lastNonBlank,trimmedLength
    integer :: firstNonBlankIn,lastNonBlankIn,trimmedLengthIn
    integer :: firstNonBlankLoc,lastNonBlankLoc,trimmedLengthLoc
    type(BudgetRecordHeaderType) :: header
    character(len=16)  :: textLabel
    character(len=16)  :: textNameLabel
    character(len=132) :: message
    logical :: lookForIFace = .false.
    integer :: ifaceindex
    integer :: nb, nbindex
    integer :: nbmax = 5
    character(len=16)  :: anamebud(5)
    DATA anamebud(1) /'           WELLS'/ ! WEL
    DATA anamebud(2) /'    DRAINS (DRT)'/ ! DRT
    DATA anamebud(3) /'          DRAINS'/ ! DRN
    DATA anamebud(4) /'   RIVER LEAKAGE'/ ! RIV
    DATA anamebud(5) /' HEAD DEP BOUNDS'/ ! GHB
    character(len=16)  :: anameid(5)
    DATA anameid(1)  /'             WEL'/ ! WEL
    DATA anameid(2)  /'             DRT'/ ! DRT
    DATA anameid(3)  /'             DRN'/ ! DRN
    DATA anameid(4)  /'             RIV'/ ! RIV
    DATA anameid(5)  /'             GHB'/ ! GHB
    !------------------------------------------------------------------------

      ! Check the iface flag
      if ( present( iFaceOption) ) then 
        lookForIface = .false.
        ifaceindex   = 0
        lookForIFace = iFaceOption
      end if 

      ! Initialize output
      isValid = .false.

      ! Trim input pkg type
      call TrimAll(sourcePkgName, firstNonBlankIn, lastNonBlankIn, trimmedLengthIn)

      ! Determine record range for stressPeriod and timeStep
      call this%BudgetReader%GetRecordHeaderRange(stressPeriod, timeStep, firstRecord, lastRecord)
      if(firstRecord .eq. 0) then
        write(message,'(A,I5,A,I5,A)') ' Error loading Time Step ', timeStep, ' Period ', stressPeriod, '.'
        message = trim(message)
        write(*,'(A)') message
        call ustop('Missing budget information. Budget file must have output for every time step. Stop.')
      end if

      naux = size( auxVarNames )
      listItemBufferSize = size(this%ListItemBuffer) 

      ! For MF6, the easiest validation is to 
      ! verify auxiliary var names by identifying the 
      ! package against TXT2ID2
      if ( isMF6 ) then
        ! Loop through record headers
        do n = firstRecord, lastRecord
          header    = this%BudgetReader%GetRecordHeader(n)
          ! Only methods 5,6 support aux variables
          if ( ( header%Method .eq. 5 ) .or. ( header%Method .eq. 6 ) ) then
            textLabel = header%TXT2ID2
            call TrimAll(textLabel, firstNonBlank, lastNonBlank, trimmedLength)
            if (&
              textLabel(firstNonBlank:lastNonBlank) .eq. & 
              sourcePkgName(firstNonBlankIn:lastNonBlankIn) ) then
              nval = 0
              do nx =1, naux 
                auxindex = header%FindAuxiliaryNameIndex(auxVarNames(nx)) ! it does a trim
                if(auxindex .gt. 0) then
                  nval = nval + 1
                end if
              end do
              if ( lookForIFace ) then
                ifaceindex = header%FindAuxiliaryNameIndex('IFACE')
                ! Is valid
                if ( (nval .eq. naux) .and. (ifaceindex.gt.0) ) then 
                  isValid = .true.
                  ! Leave 
                  return
                end if
              else
                ! Is valid
                if ( nval .eq. naux ) then 
                  isValid = .true.
                  ! Leave 
                  return
                end if
              end if 
            end if
          end if
        end do
      ! Other MODFLOW flavors
      ! Compare against the list of known headers that accept
      ! aux variables
      else
        ! Loop through record headers
        do n = firstRecord, lastRecord
          header    = this%BudgetReader%GetRecordHeader(n)
          ! Only methods 5,6 support aux variables
          if ( ( header%Method .eq. 5 ) .or. ( header%Method .eq. 6 ) ) then
            textLabel = header%TextLabel
            call TrimAll(textLabel, firstNonBlank, lastNonBlank, trimmedLength)

            ! Needs to verify relation
            ! Find equivalence
            nbindex = 0
            do nb=1,nbmax
              textNameLabel = anamebud(nb) 
              call TrimAll(textNameLabel, firstNonBlankLoc, lastNonBlankLoc, trimmedLengthLoc)
              if (&
                textNameLabel(firstNonBlankLoc:lastNonBlankLoc) .eq. & 
                textLabel(firstNonBlank:lastNonBlank) ) then
                ! Found, continue
                nbindex = nb
                exit
              end if   
            end do

            ! Not found in the list of known budgets supporting
            ! aux variables, try next budget header. It might be useful 
            ! to report the header text label for validation.
            if ( nbindex .eq. 0 ) then
              exit
            end if

            ! Compare the id/ftype (e.g. WEL) against the given src name,
            ! and if not, give it another chance by comparing against the 
            ! budget label itself (e.g. WELLS)
            textNameLabel = anameid(nbindex) 
            call TrimAll(textNameLabel, firstNonBlankLoc, lastNonBlankLoc, trimmedLengthLoc)
            if (&
              textNameLabel(firstNonBlankLoc:lastNonBlankLoc) .eq. & 
              sourcePkgName(firstNonBlankIn:lastNonBlankIn) ) then
              nval = 0
              do nx =1, naux 
                auxindex = header%FindAuxiliaryNameIndex(auxVarNames(nx)) ! it does a trim
                if(auxindex .gt. 0) then
                  nval = nval + 1
                end if
              end do
              if ( lookForIFace ) then
                ifaceindex = header%FindAuxiliaryNameIndex('IFACE')
                print *, 'case1', ifaceindex
                ! Is valid
                if ( (nval .eq. naux) .and. (ifaceindex.gt.0) ) then 
                  isValid = .true.
                  ! Leave 
                  return
                end if
              else
                ! Is valid
                if ( nval .eq. naux ) then 
                  isValid = .true.
                  ! Leave 
                  return
                end if
              end if
            else
              textNameLabel = anamebud(nbindex) 
              call TrimAll(textNameLabel, firstNonBlankLoc, lastNonBlankLoc, trimmedLengthLoc)
              if (&
                textNameLabel(firstNonBlankLoc:lastNonBlankLoc) .eq. & 
                sourcePkgName(firstNonBlankIn:lastNonBlankIn) ) then
                nval = 0
                do nx =1, naux 
                  auxindex = header%FindAuxiliaryNameIndex(auxVarNames(nx)) ! it does a trim
                  if(auxindex .gt. 0) then
                    nval = nval + 1
                  end if
                end do
                if ( lookForIFace ) then
                  ifaceindex = header%FindAuxiliaryNameIndex('IFACE')
                  print *, 'case2', ifaceindex
                  ! Is valid
                  if ( (nval .eq. naux) .and. (ifaceindex.gt.0) ) then 
                    isValid = .true.
                    ! Leave 
                    return
                  end if
                else
                  ! Is valid
                  if ( nval .eq. naux ) then 
                    isValid = .true.
                    ! Leave 
                    return
                  end if
                end if
              end if
            end if
          end if ! ( header%Method .eq. 5 ) .or. ( header%Method .eq. 6 )

        end do ! n = firstRecord, lastRecord

      end if 


      ! Done
      return


    end function pr_ValidateAuxVarNames






    !! DEPRECATION WARNING !!


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




    ! TO BE DEPRECATED!
    subroutine pr_LoadFlowTimeSeries(this, initialTime, finalTime, cellNumbers, tdisData )
    use TimeDiscretizationDataModule,only : TimeDiscretizationDataType
    !------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------
    ! Specifications
    !------------------------------------------------------------------------
    implicit none
    ! input
    class(FlowModelDataType) :: this
    doubleprecision, intent(in) :: initialTime, finalTime
    integer, dimension(:), intent(in) :: cellNumbers
    type( TimeDiscretizationDataType ), intent(in) :: tdisData
    ! should be obtained from boundary
    character(len=16), allocatable, dimension(:) :: auxnames  
    ! local
    integer :: stressPeriod, timeStep
    integer :: firstRecord,lastRecord,n,m,firstNonBlank,lastNonBlank,trimmedLength
    integer :: spaceAssigned, status, cellCount, iface, index, auxindex, cellindex
    integer :: boundaryFlowsOffset, listItemBufferSize, cellNumber, layer
    type(BudgetRecordHeaderType) :: header
    character(len=16) :: textLabel
    character(len=132) message
    integer :: nCells
    integer :: kinitial, kfinal, ktime, kcounter
    integer :: nTimes 
    doubleprecision, allocatable, dimension(:,:) :: timeseries
    doubleprecision, allocatable, dimension(:,:) :: ctimeseries
    integer, allocatable, dimension(:) :: ifaces
    !------------------------------------------------------------------------

        print *, '-------------------------------------------------------'
        print *, 'FLOWMODELDATA:: LOADFLOWTIMESERIES'
        print *, 'FLOWMODELDATA:: LOADFLOWTIMESERIES ( OR LOADBUDGET...) '
        print *, '-------------------------------------------------------'

      ! Seems needed
      listItemBufferSize = size(this%ListItemBuffer)

      ! While defining cellNumbers, shall we determine 
      ! in advance on which BudgetRecords are we intererested ?

      nCells = size(cellNumbers)

      if ( allocated( ifaces ) ) deallocate( ifaces )
      allocate( ifaces( nCells ) )
      if ( allocated( auxnames ) ) deallocate( auxnames )
      allocate( auxnames( nCells ) )
      auxnames(1) = 'CONCENTRATION'
      ifaces(:) = 0

      ! Given initial and final times, 
      ! compute an array of times from which 
      ! determine and extract stress periods/time steps
      kinitial = tdisData%FindContainingTimeStep(initialTime)
      kfinal   = tdisData%FindContainingTimeStep(finalTime)

      do ktime=kinitial,kfinal
        if ( ktime .eq. kinitial ) print *, 'TINITIAL::',initialTime
        if ( ktime .lt. kfinal   ) print *, 'KTIME   ::',tdisData%TotalTimes(ktime)
        if ( ktime .eq. kfinal   ) print *, 'TEND    ::',finalTime 
      end do 

      ! nTimes: This number represent the number of intervals
      nTimes = kfinal - kinitial + 1 

      print *, '-------------------------------------------------'
      if ( allocated( timeseries ) ) deallocate( timeseries ) 
      allocate( timeseries( nTimes, nCells ) )
      timeseries(:,:) = 0d0
      if ( allocated( ctimeseries ) ) deallocate( ctimeseries ) 
      allocate( ctimeseries( nTimes, nCells ) )
      ctimeseries(:,:) = 0d0

      print *, 'STRESS PERIOD AND TIME STEP: '
      kcounter = 0
      do ktime=kinitial,kfinal
        ! Get the stress period and time step from the cummulative time step
        call tdisData%GetPeriodAndStep(ktime, stressPeriod, timeStep)
        kcounter = kcounter + 1 
        print *, 'KINFO:: ', kcounter, ktime, stressPeriod, timeStep
        print *, '-------------------------------------------------'

        ! Determine record range for stressPeriod and timeStep
        call this%BudgetReader%GetRecordHeaderRange(stressPeriod, timeStep, firstRecord, lastRecord)
        if(firstRecord .eq. 0) then
          write(message,'(A,I5,A,I5,A)') ' Error loading Time Step ', timeStep, ' Period ', stressPeriod, '.'
          message = trim(message)
          write(*,'(A)') message
          call ustop('Missing budget information. Budget file must have output for every time step. Stop.')
        end if

        print *, '----------------------------------------'
        print *, ' LOOP OVER BUDGET HEADERS'
        print *, '----------------------------------------'

        !! Loop through record headers
        do n = firstRecord, lastRecord
          header    = this%BudgetReader%GetRecordHeader(n)
          textLabel = header%TextLabel
          call TrimAll(textLabel, firstNonBlank, lastNonBlank, trimmedLength)
          select case(textLabel(firstNonBlank:lastNonBlank))
          case('CONSTANT HEAD', 'CHD')
            ! Do not skip ( some understanding of IFACE? )
            print *, '----------TEXTLABEL: ', textLabel(firstNonBlank:lastNonBlank), ' - HEADERMETHOD', header%Method

      !    ! Read constant head flows into the sinkFlows and sourceFlows arrays.
      !    ! For a standard budget file, Method = 0. For a compact budget file,
      !    ! Method = 2.
      !    if(header%Method .eq. 0) then
      !      call this%BudgetReader%FillRecordDataBuffer(header,       &
      !        this%ArrayBufferDbl, cellCount, spaceAssigned, status)
      !      if(cellCount .eq. spaceAssigned) then
      !        do m = 1, spaceAssigned
      !          if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
      !            this%SourceFlows(m) = this%SourceFlows(m) +         &
      !              this%ArrayBufferDbl(m)
      !          else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
      !            this%SinkFlows(m) = this%SinkFlows(m) +             &
      !              this%ArrayBufferDbl(m)
      !          end if
      !        end do
      !      end if
      !    else if(header%Method .eq. 2) then
      !      call this%BudgetReader%FillRecordDataBuffer(header,             &
      !        this%ListItemBuffer, listItemBufferSize, spaceAssigned, status)
      !      if(spaceAssigned .gt. 0) then
      !        do m = 1, spaceAssigned
      !          cellNumber = this%ListItemBuffer(m)%CellNumber
      !          if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
      !            this%SourceFlows(cellNumber) =                      &
      !              this%SourceFlows(cellNumber) + this%ListItemBuffer(m)%BudgetValue
      !          else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
      !            this%SinkFlows(cellNumber) =                        &
      !              this%SinkFlows(cellNumber) + this%ListItemBuffer(m)%BudgetValue
      !          end if
      !        end do
      !      end if
      !    else if((header%Method .eq. 5) .or. (header%Method .eq. 6)) then
      !      call this%BudgetReader%FillRecordDataBuffer(header,             &
      !        this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
      !        status)
      !      if(spaceAssigned .gt. 0) then
      !        do m = 1, spaceAssigned
      !          call this%CheckForDefaultIface(header%TextLabel, iface)
      !          index = header%FindAuxiliaryNameIndex('IFACE')
      !          if(index .gt. 0) then
      !            iface = int(this%ListItemBuffer(m)%AuxiliaryValues(index))
      !          end if
      !          
      !          cellNumber = this%ListItemBuffer(m)%CellNumber
      !          if(iface .gt. 0) then
      !            boundaryFlowsOffset = 6 * (cellNumber - 1)
      !            this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
      !              this%BoundaryFlows(boundaryFlowsOffset + iface) + &
      !              this%ListItemBuffer(m)%BudgetValue
      !          else
      !            if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
      !              this%SourceFlows(cellNumber) =                  &
      !                this%SourceFlows(cellNumber) +                &
      !                this%ListItemBuffer(m)%BudgetValue
      !            else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
      !              this%SinkFlows(cellNumber) =                    &
      !                this%SinkFlows(cellNumber) +                  &
      !                this%ListItemBuffer(m)%BudgetValue
      !            end if
      !          end if
      !        end do
      !      end if
      !    end if

            continue

          case('STORAGE', 'STO-SS', 'STO-SY')
            ! Skip
            cycle
          case('FLOW JA FACE', 'FLOW-JA-FACE')
            ! Skip
            cycle
          case('DATA-SPDIS', 'DATA-SAT')
            ! Skip
            cycle
          case('FLOW RIGHT FACE', 'FLOW FRONT FACE', 'FLOW LOWER FACE')
            ! Skip
            cycle
          case default
            ! Do not skip ( some understanding of IFACE? )
            print *, '----------TEXTLABEL: ', textLabel(firstNonBlank:lastNonBlank), ' - HEADERMETHOD', header%Method

       !    ! Now handle any other records in the budget file.
       !    if((header%Method .eq. 0) .or. (header%Method .eq. 1)) then
       !      if(header%ArrayItemCount .eq. cellCount) then
       !        call this%BudgetReader%FillRecordDataBuffer(header,         &
       !          this%ArrayBufferDbl, cellCount, spaceAssigned, status)
       !        if(cellCount .eq. spaceAssigned) then
       !          call this%CheckForDefaultIface(header%TextLabel, iface)
       !          if(iface .gt. 0) then
       !            do m = 1, spaceAssigned
       !              boundaryFlowsOffset = 6 * (m - 1)
       !              this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
       !                this%BoundaryFlows(boundaryFlowsOffset + iface) + &
       !                this%ArrayBufferDbl(m)
       !            end do
       !          else
       !            do m = 1, spaceAssigned
       !              if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
       !                this%SourceFlows(m) = this%SourceFlows(m) +     &
       !                  this%ArrayBufferDbl(m)
       !              else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
       !                this%SinkFlows(m) = this%SinkFlows(m) +         &
       !                  this%ArrayBufferDbl(m)
       !              end if
       !            end do
       !          end if
       !        end if
       !      end if
       !    else if(header%Method .eq. 3) then
       !      call this%BudgetReader%FillRecordDataBuffer(header,             &
       !        this%ArrayBufferDbl, this%ArrayBufferInt,                     &
       !        header%ArrayItemCount, spaceAssigned, status)
       !      if(header%ArrayItemCount .eq. spaceAssigned) then
       !        call this%CheckForDefaultIface(header%TextLabel, iface)
       !        if(iface .gt. 0) then
       !          do m = 1, spaceAssigned
       !            if(this%Grid%GridType .ne. 2) then
       !              ! structured
       !              layer = this%ArrayBufferInt(m)
       !              cellNumber = (layer - 1) * spaceAssigned + m
       !            else
       !              ! mfusg unstructured
       !              cellNumber = this%ArrayBufferInt(m)
       !            end if
       !            boundaryFlowsOffset = 6 * (cellNumber - 1)
       !            this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
       !              this%BoundaryFlows(boundaryFlowsOffset + iface) + &
       !              this%ArrayBufferDbl(m)
       !          end do
       !        else            
       !          do m = 1, spaceAssigned
       !            cellNumber = this%ArrayBufferInt(m)
       !            if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
       !              this%SourceFlows(cellNumber) =                  &
       !                this%SourceFlows(cellNumber) +                &
       !                this%ArrayBufferDbl(m)
       !            else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
       !              this%SinkFlows(cellNumber) =                    &
       !                this%SinkFlows(cellNumber) +                  &
       !                this%ArrayBufferDbl(m)
       !            end if
       !          end do
       !        end if
       !      end if
       !    else if(header%Method .eq. 4) then
       !      call this%BudgetReader%FillRecordDataBuffer(header,             &
       !        this%ArrayBufferDbl, header%ArrayItemCount, spaceAssigned,    &
       !        status)
       !      if(header%ArrayItemCount .eq. spaceAssigned) then
       !        call this%CheckForDefaultIface(header%TextLabel, iface)
       !        if(iface .gt. 0) then
       !          do m = 1, spaceAssigned
       !            boundaryFlowsOffset = 6 * (m - 1)
       !            this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
       !              this%BoundaryFlows(boundaryFlowsOffset + iface) + &
       !              this%ArrayBufferDbl(m)
       !          end do
       !        else            
       !          do m = 1, spaceAssigned
       !            if(this%ArrayBufferDbl(m) .gt. 0.0d0) then
       !              this%SourceFlows(m) = this%SourceFlows(m) +     &
       !                this%ArrayBufferDbl(m)
       !            else if(this%ArrayBufferDbl(m) .lt. 0.0d0) then
       !              this%SinkFlows(m) = this%SinkFlows(m) +         &
       !                this%ArrayBufferDbl(m)
       !            end if
       !          end do
       !        end if
       !      end if
       !    else if(header%Method .eq. 2) then
       !      call this%BudgetReader%FillRecordDataBuffer(header,             &
       !        this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
       !        status)
       !      if(spaceAssigned .gt. 0) then
       !        call this%CheckForDefaultIface(header%TextLabel, iface)
       !        if(iface .gt. 0) then
       !          do m = 1, spaceAssigned
       !              cellNumber = this%ListItemBuffer(m)%CellNumber
       !              boundaryFlowsOffset = 6 * (cellNumber - 1)
       !              this%BoundaryFlows(boundaryFlowsOffset + iface) =   &
       !                this%BoundaryFlows(boundaryFlowsOffset + iface) + &
       !                this%ListItemBuffer(m)%BudgetValue
       !          end do
       !        else            
       !          do m = 1, spaceAssigned
       !            cellNumber = this%ListItemBuffer(m)%CellNumber
       !            if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
       !              this%SourceFlows(cellNumber) =                  &
       !                this%SourceFlows(cellNumber) +                &
       !                this%ListItemBuffer(m)%BudgetValue
       !            else if(this%ListItemBuffer(m)%BudgetValue .lt. 0.0d0) then
       !              this%SinkFlows(cellNumber) =                    &
       !                this%SinkFlows(cellNumber) +                  &
       !                this%ListItemBuffer(m)%BudgetValue
       !            end if
       !          end do
       !        end if
       !      end if
       !    else if((header%Method .eq. 5) .or. (header%Method .eq. 6)) then
            if((header%Method .eq. 5) .or. (header%Method .eq. 6)) then
              call this%BudgetReader%FillRecordDataBuffer(header,             &
                this%ListItemBuffer, listItemBufferSize, spaceAssigned,       &
                status)
              if(spaceAssigned .gt. 0) then
                do m = 1, spaceAssigned



                  cellNumber = this%ListItemBuffer(m)%CellNumber
                  ! Determine the index of cellNumber in the list of cells 
                  ! requested for timeseseries
                  cellindex = findloc( cellNumbers, cellNumber, 1 ) 
                  if ( cellindex .eq. 0 ) cycle ! Not requested
                  print *, cellindex, cellNumber
                  print *, header%AuxiliaryNames 
                  print *, this%ListItemBuffer(m)%AuxiliaryValues


                  call this%CheckForDefaultIface(header%TextLabel, iface)
                  index = header%FindAuxiliaryNameIndex('IFACE')
                  if(index .gt. 0) then
                    iface = int(this%ListItemBuffer(m)%AuxiliaryValues(index))
                  end if
                  if(iface .gt. 0) ifaces(cellindex) = iface

                  ! Load into flow rates timeseries only if positive, 
                  ! otherwise leave as zero
                  if(this%ListItemBuffer(m)%BudgetValue .gt. 0.0d0) then
                    timeseries( kcounter, cellindex )  = this%ListItemBuffer(m)%BudgetValue
                  end if

                  ! Load into concentration timeseries only if positive, 
                  ! otherwise leave as zero

                  ! NOTE: this should load all the aux variables of this cell
                  auxindex = header%FindAuxiliaryNameIndex(auxnames(cellindex))
                  if(auxindex .gt. 0) then
                    ctimeseries( kcounter, cellindex ) = this%ListItemBuffer(m)%AuxiliaryValues(cellindex)
                  end if

                end do
              end if
            end if

            continue

          end select
        end do

      end do


      print *, '-----------------------------------'
      print *, ' THE TIMESERIES !!                 '
      print *, '-----------------------------------'
      kcounter = 0
      do ktime = kinitial, kfinal
        kcounter = kcounter + 1 
        print *, kcounter, ktime, timeseries(kcounter, 1), ctimeseries(kcounter, 1) 
      end do


      print *, 'LEAVING!!'
      call exit(0)


    end subroutine pr_LoadFlowTimeSeries


end module FlowModelDataModule
