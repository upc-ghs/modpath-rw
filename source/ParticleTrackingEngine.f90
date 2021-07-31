module ParticleTrackingEngineModule
  use TrackPathResultModule,only : TrackPathResultType
  use ParticleLocationModule,only : ParticleLocationType
  use ParticleLocationListModule,only : ParticleLocationListType
  use ParticleCoordinateModule,only : ParticleCoordinateType
  use TrackCellModule,only : TrackCellType
  use TrackCellResultModule,only : TrackCellResultType
  use BudgetReaderModule,only : BudgetReaderType
  use HeadReaderModule,only : HeadReaderType
  use BudgetListItemModule,only : BudgetListItemType
  use ModflowRectangularGridModule,only : ModflowRectangularGridType
  use ModpathCellDataModule,only : ModpathCellDataType
  use ModpathSubCellDataModule,only : ModpathSubCellDataType
  use ParticleTrackingOptionsModule,only : ParticleTrackingOptionsType
  use BudgetRecordHeaderModule,only : BudgetRecordHeaderType
  use UtilMiscModule,only : TrimAll
!  use ModpathUnstructuredBasicDataModule,only : ModpathUnstructuredBasicDataType
  implicit none
  
! Set default access status to private
  !private

  type,public :: ParticleTrackingEngineType
!    doubleprecision :: ReferenceTime = 0d0
!    doubleprecision :: StoppingTime = 0d0
    doubleprecision :: HDry = 0d0
    doubleprecision :: HNoFlow = 0d0
    type(ParticleTrackingOptionsType) :: TrackingOptions
    !type(ModpathCellDataType) :: CellDataBuffer
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
    
    
    ! Private variables
    type(HeadReadertype),pointer :: HeadReader => null()
    type(BudgetReaderType),pointer :: BudgetReader => null()
    class(ModflowRectangularGridType),pointer :: Grid => null()
    type(TrackCellType),private :: TrackCell
    type(TrackCellResultType),private :: TrackCellResult
    integer,private :: CurrentStressPeriod = 0
    integer,private :: CurrentTimeStep = 0
    integer,private :: MaxReducedCellConnectionsCount = 17
    doubleprecision,dimension(17),private :: CellFlowsBuffer
    integer,dimension(17),private :: ReducedCellConnectionsBuffer
    type(BudgetListItemType),allocatable,dimension(:) :: ListItemBuffer
    logical,allocatable,dimension(:),private :: SubFaceFlowsComputed
    type(ParticleLocationListType) :: LocBuffP
    type(ParticleLocationListType) :: LocBuffTS

    ! RWPT
    type(ModpathSubCellDataType), dimension(18) :: NeighborSubCellData
    type(ModpathCellDataType), dimension(18) :: NeighborCellData

  contains
    procedure :: Initialize=>pr_Initialize
    procedure :: Reset=>pr_Reset
    procedure :: ClearTimeStepBudgetData=>pr_ClearTimeStepBudgetData
    procedure :: TrackPath=>pr_TrackPath
    procedure :: FillCellBuffer=>pr_FillCellBuffer
    procedure :: LoadTimeStep=>pr_LoadTimeStep
    procedure :: FillFaceFlowsBuffer=>pr_FillFaceFlowsBuffer
    procedure :: GetCurrentStressPeriod=>pr_GetCurrentStressPeriod
    procedure :: GetCurrentTimeStep=>pr_GetCurrentTimeStep
    procedure :: FillCellFlowsBuffer=>pr_FillCellFlowsBuffer
    procedure :: SetIBound=>pr_SetIBound
    procedure :: SetZones=>pr_SetZones
    procedure :: SetPorosity=>pr_SetPorosity
    procedure :: SetRetardation=>pr_SetRetardation
    !procedure :: SetLayerTypes=>pr_SetLayerTypes
    procedure :: SetDefaultIface=>pr_SetDefaultIface
    procedure :: CheckForDefaultIface=>pr_CheckForDefaultIface
    procedure :: GetVolumetricBalanceSummary=>pr_GetVolumetricBalanceSummary
    procedure :: WriteCellBuffer=>pr_WriteCellBuffer
    procedure :: GetTopMostActiveCell=>pr_GetTopMostActiveCell

    ! RWPT
    procedure :: FillNeighborSubCellData=>pr_FillNeighborSubCellData
    procedure :: FillNeighborCellData=>pr_FillNeighborCellData

  end type
  
contains

!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
!  implicit none
!---------------------------------------------------------------------------------------------------------------
function pr_GetTopMostActiveCell(this, cellNumber) result(n)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
  integer,intent(in) :: cellNumber
  integer :: n
  
  n = cellNumber
  do while(.true.)
      if(n .eq. 0) return
      if(this%IboundTS(n) .ne. 0) return
      n = this%Grid%GetFaceConnection(n, 5, 1)
  end do
      
end function pr_GetTopMostActiveCell

subroutine pr_WriteCellBuffer(this, unit, cellData, backwardTracking)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
  integer,intent(in) :: unit
  logical,intent(in) :: backwardTracking
  type(ModpathCellDataType),intent(in) :: cellData
!---------------------------------------------------------------------------------------------------------------
  
  write(unit, '(1X,A)')     '-------------------------------------'
  write(unit, '(1X,A,I10)') '      Cell', cellData%CellNumber
  write(unit, '(1X,A)')     '-------------------------------------'
  call WriteCellData(unit, cellData, backwardTracking)
  
end subroutine pr_WriteCellBuffer

subroutine WriteCellData(unit, cellData, backwardTracking)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  integer,intent(in) :: unit
  logical,intent(in) :: backwardTracking
  type(ModpathCellDataType),intent(in) :: cellData
  type(ModpathSubCellDataType) :: subCellData
  integer :: m, n, count, row, column, subRowCount, subColumnCount
  doubleprecision :: balance, totalFaceInflow, totalFaceOutflow, sourceInflow,  &
    sinkOutflow, storageInflow, storageOutflow, netFaceInflow
  doubleprecision :: totalInflow, totalOutflow, netInflow
!---------------------------------------------------------------------------------------------------------------
  write(unit, '(1X,A,5I10)')   'Layer, Ibound, IboundTS, Zone, LayerType: ',    &
    cellData%Layer, cellData%Ibound, cellData%IboundTS, cellData%Zone,          &
    cellData%LayerType
  write(unit, '(1X,A,4E15.7)') 'DX, DY, MinX, MinY:    ', cellData%DX,          &
    cellData%DY, cellData%MinX, cellData%MinY
  write(unit, '(1X,A,3E15.7)') 'Bottom, Top, Head:     ', cellData%Bottom,      &
    cellData%Top, cellData%Head
  write(unit, '(1x,A,2E15.7)') 'Porosity, Retardation: ', cellData%Porosity,    &
    cellData%Retardation
  
  write(unit, *)
  write(unit, '(1X,A)') 'Volumetric Face Flows (L**3/T):'
  write(unit, '(17X,4(15X,A))') 'sub-face 1', 'sub-face 2', 'sub-face 3', 'sub-face 4'
  write(unit, '(1X,A,4E25.15)') 'Left   (face 1):',(cellData%GetFaceFlow(1,n), n = 1, cellData%GetSubFaceCount(1))
  write(unit, '(1X,A,4E25.15)') 'Right  (face 2):',(cellData%GetFaceFlow(2,n), n = 1, cellData%GetSubFaceCount(2))
  write(unit, '(1X,A,4E25.15)') 'Front  (face 3):',(cellData%GetFaceFlow(3,n), n = 1, cellData%GetSubFaceCount(3))
  write(unit, '(1X,A,4E25.15)') 'Back   (face 4):',(cellData%GetFaceFlow(4,n), n = 1, cellData%GetSubFaceCount(4))
  write(unit, '(1X,A,4E25.15)') 'Bottom (face 5):',(cellData%GetFaceFlow(5,n), n = 1, cellData%GetSubFaceCount(5))
  write(unit, '(1X,A,4E25.15)') 'Top    (face 6):',(cellData%GetFaceFlow(6,n), n = 1, cellData%GetSubFaceCount(6))
  
  call cellData%GetVolumetricBalanceComponents(totalFaceInflow,                 &
    totalFaceOutflow, sourceInflow, sinkOutflow, storageInflow, storageOutflow, balance)
  totalInflow = totalFaceInflow + sourceInflow + storageInflow
  totalOutflow = totalFaceOutflow + sinkOutflow + storageOutflow
  netInflow = totalInflow - totalOutflow
  write(unit, *)
  write(unit, '(1X,A)') 'Water balance components:'
  write(unit, '(1X,A)') 'Inflow (L**3/T)'
  write(unit, '(1X,A,E25.15)')    '     Total face inflow =', totalFaceInflow
  write(unit, '(1X,A,E25.15)')    '         Source inflow =', sourceInflow
  write(unit, '(1X,A,E25.15)')    '        Storage inflow =', storageInflow
  write(unit, '(27X,A)')          '-----------------------'
  write(unit, '(1X,A,E25.15)')    '               Inflow  =', totalInflow
  write(unit, *)
  write(unit, '(1X,A)') 'Outflow (L**3/T)'
  write(unit, '(1X,A,E25.15)')    '    Total face outflow =', totalFaceOutflow
  write(unit, '(1X,A,E25.15)')    '          Sink outflow =', sinkOutflow
  write(unit, '(1X,A,E25.15)')    '       Storage outflow =', storageOutflow
  write(unit, '(27X,A)')          '-----------------------'
  write(unit, '(1X,A,E25.15)')    '                Ouflow =', totalOutflow
  write(unit, *)
  write(unit, '(1X,A,E25.15)')    '      Inflow - Outflow =', netInflow
  write(unit, *)
  write(unit, '(1X,A,E25.15)')    '   Percent discrepancy =', balance
  
  write(unit, *)
  if(backwardTracking) then
      write(unit, '(1X,A)')                                                     &
        'Face velocities (Backward tracking. Velocity components has been reversed to represent backward tracking.)'
  else
      write(unit, '(1X,A)') 'Face velocities (Forward tracking)'
  end if
  subRowCount = cellData%GetSubCellRowCount()
  subColumnCount = cellData%GetSubCellColumnCount()
  do row = 1, subRowCount
      do column = 1, subColumnCount
          call cellData%FillSubCellDataBuffer(subCellData, row, column, backwardTracking)
          write(unit, '(1X,A,I2,A,I2,A)') 'SubCell (', row, ', ', column, ')'
          write(unit, '(1X,A,3E15.7)')    'DX, DY, DZ: ', subCellData%DX, subCellData%DY, subCellData%DZ
          write(unit, '(23X,A,5X,A)') 'Face Velocity (L/T)', 'Connection'  
          write(unit, '(1X,A,E25.15,I15)') 'Left   (face 1):', subCellData%VX1, subCellData%Connection(1)
          write(unit, '(1X,A,E25.15,I15)') 'Right  (face 2):', subCellData%VX2, subCellData%Connection(2)
          write(unit, '(1X,A,E25.15,I15)') 'Front  (face 3):', subCellData%VY1, subCellData%Connection(3)
          write(unit, '(1X,A,E25.15,I15)') 'Back   (face 4):', subCellData%VY2, subCellData%Connection(4)
          write(unit, '(1X,A,E25.15,I15)') 'Bottom (face 5):', subCellData%VZ1, subCellData%Connection(5)
          write(unit, '(1X,A,E25.15,I15)') 'Top    (face 6):', subCellData%VZ2, subCellData%Connection(6)
          write(unit, *)
      end do
  end do

end subroutine WriteCellData

subroutine WriteTraceData(unit, trackCell, tcResult, stressPeriod, timeStep)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  type(TrackCellType),intent(in),target :: trackCell
  type(TrackCellResultType),intent(in),target :: tcResult
  type(ModpathCellDataType),pointer :: cellData
  type(ModpathSubCellDataType) :: subCellData
  integer,intent(in) :: unit, stressPeriod, timeStep
  integer :: m, n, count, row, column, subRowCount, subColumnCount
  character(len=28) :: statusText
!---------------------------------------------------------------------------------------------------------------
  cellData => trackCell%CellData
  count = tcResult%TrackingPoints%GetItemCount()          
  write(unit, *)
  write(unit, '(1X,A,I10,A,I6,A,I6)') '----- Call TrackCell: Cell',             &
    cellData%CellNumber, ',  Stress period',stressPeriod,                       &
    ',  Time step', timeStep

  select case (tcResult%Status)
      case (1)
          statusText = '  (Reached stopping time)'
      case (2)
          statusText = '  (Exit at cell face)'
      case (3)
          statusText = '  (Stop at weak sink)'
      case (4)
          statusText = '  (Stop at weak source)'
      case (5)
          statusText = '  (No exit possible)'
      case (6)
          statusText = '  (Stop zone cell)'
      case (7)
          statusText = '  (Inactive cell)'
      case (8)
          statusText = '  (Inactive cell)'
      case default
          statusText = '  (Undefined)'
  end select
      
  write(unit, '(1X,A,I3,A)') 'Exit status =', tcResult%Status, statusText
  write(unit, '(1X,A)')         'Particle locations: Local X, Local Y, Local Z, Tracking time'
  do n = 1, count
      write(unit, '(1X,I10,4E25.15,I10)')                                       &
        tcResult%TrackingPoints%Items(n)%CellNumber,                            &
        tcResult%TrackingPoints%Items(n)%LocalX,                                &
        tcResult%TrackingPoints%Items(n)%LocalY,                                &
        tcResult%TrackingPoints%Items(n)%LocalZ,                                &
        tcResult%TrackingPoints%Items(n)%TrackingTime,                          &
        tcResult%TrackingPoints%Items(n)%Layer              
  end do
  write(unit, *)
  
  call WriteCellData(unit, cellData, trackCell%TrackingOptions%BackwardTracking)
  
end subroutine WriteTraceData

subroutine pr_GetVolumetricBalanceSummary(this, intervalCount, intervalBreaks,  &
  balanceCounts, maxError, maxErrorCell)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
  integer, intent(in) :: intervalCount
  doubleprecision, dimension(intervalCount), intent(in) :: intervalBreaks
  integer, dimension(intervalCount + 1), intent(inout) :: balanceCounts
  integer, intent(inout) :: maxErrorCell
  doubleprecision, intent(inout) :: maxError
  type(ModpathCellDataType) :: cellBuffer
  integer :: cellCount, n, m
  doubleprecision :: balance, absBalance
  
  
  cellCount = this%Grid%CellCount
  
  maxErrorCell = 0
  maxError = 0.0d0
  do m = 1, intervalCount + 1
      balanceCounts(m) = 0
  end do
  
  do n = 1, cellCount
      call this%FillCellBuffer(n, cellBuffer)
      if(cellBuffer%IboundTS .gt. 0) then
          balance = cellBuffer%GetVolumetricBalance()
          absBalance = dabs(balance)
          if((maxErrorCell .eq. 0) .or. (absBalance .gt. maxError) ) then
              maxError = absBalance
              maxErrorCell = n
          end if
          do m = 1, intervalCount
              if(absBalance .le. intervalBreaks(m)) then
                  balanceCounts(m) = balanceCounts(m) + 1
                  exit
              end if
              if(m .eq. intervalCount) balanceCounts(intervalCount + 1) =       &
                balanceCounts(intervalCount + 1) + 1
          end do
      end if
  end do

end subroutine pr_GetVolumetricBalanceSummary

subroutine pr_SetDefaultIface(this, defaultIfaceLabels, defaultIfaceValues,     &
  arraySize)
!***************************************************************************************************************
! Description goes here
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  use UtilMiscModule,only : utrimall
  implicit none
  class(ParticleTrackingEngineType) :: this
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

!subroutine pr_SetLayerTypes(this, layerTypes, arraySize)
!!***************************************************************************************************************
!!
!!***************************************************************************************************************
!!
!! Specifications
!!---------------------------------------------------------------------------------------------------------------
!  implicit none
!  class(ParticleTrackingEngineType) :: this
!  integer,intent(in) :: arraySize
!  integer,dimension(arraySize),intent(in),target :: layerTypes
!  
!  if(arraySize .ne. this%Grid%LayerCount) then
!      write(*,*) "ParticleTrackingEngine: The LayerTypes array size does not match the layer count for the grid. stop"
!      stop
!  end if
!  
!  this%LayerTypes => layerTypes
!
!end subroutine pr_SetLayerTypes

subroutine pr_SetIbound(this, ibound, arraySize)
!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
  integer,intent(in) :: arraySize
  integer :: n
  integer,dimension(arraySize),intent(in),target :: ibound
  
  if(arraySize .ne. this%Grid%CellCount) then
      write(*,*) "ParticleTrackingEngine: The IBound array size does not match the cell count for the grid. stop"
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
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
  integer,intent(in) :: arraySize
  integer,dimension(arraySize),intent(in),target :: zones
  
  if(arraySize .ne. this%Grid%CellCount) then
      write(*,*) "ParticleTrackingEngine: The Zones array size does not match the cell count for the grid. stop"
      stop
  end if
  
  this%Zones => zones

end subroutine pr_SetZones

subroutine pr_SetPorosity(this, porosity, arraySize)
!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
  integer,intent(in) :: arraySize
  doubleprecision,dimension(arraySize),intent(in),target :: porosity
  
  if(arraySize .ne. this%Grid%CellCount) then
      write(*,*) "ParticleTrackingEngine: The Porosity array size does not match the cell count for the grid. stop"
      stop
  end if
  
  this%Porosity => porosity

end subroutine pr_SetPorosity

subroutine pr_SetRetardation(this, retardation, arraySize)
!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
  integer,intent(in) :: arraySize
  doubleprecision,dimension(arraySize),intent(in),target :: retardation
  
  if(arraySize .ne. this%Grid%CellCount) then
      write(*,*) "ParticleTrackingEngine: The Retardation array size does not match the cell count for the grid. stop"
      stop
  end if
  
  this%Retardation => retardation

end subroutine pr_SetRetardation

function pr_FindTimeIndex(timeSeriesPoints, currentTime, maximumTime, timePointsCount) result(index)
!***************************************************************************************************************
! Find the index in the timeSeriesPoints array of the next stopping time after the currentTime value. 
! Return -1 if none is found or if maximumTime is less than the stopping time found in the timeSeriesPoints 
! array.
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  integer,intent(in) :: timePointsCount
  doubleprecision,dimension(timePointsCount),intent(in) :: timeSeriesPoints
  doubleprecision,intent(in) :: currentTime, maximumTime
  integer :: index, n
  doubleprecision :: t
!---------------------------------------------------------------------------------------------------------------
  index = -1
  
  if(timePointsCount .lt. 1) return
  do n = 1, timePointsCount
      if((timeSeriesPoints(n) .gt. currentTime) .and. (timeSeriesPoints(n) .le. maximumTime)) then
          index = n
          return
      end if
  end do

end function pr_FindTimeIndex

subroutine pr_Initialize(this,headReader, budgetReader, grid, hNoFlow, hDry, trackingOptions)
!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
  type(BudgetReaderType),intent(inout),target :: budgetReader
  type(HeadReaderType),intent(inout),target :: headReader
  class(ModflowRectangularGridType),intent(inout),pointer :: grid
  type(ParticleTrackingOptionsType),intent(in) :: trackingOptions
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
  this%TrackingOptions = trackingOptions
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

subroutine  pr_FillFaceFlowsBuffer(this,buffer,bufferSize,count)
!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
  integer,intent(in) :: bufferSize
  doubleprecision,intent(inout),dimension(bufferSize) :: buffer
  integer,intent(inout) :: count
  integer :: n,offset
!---------------------------------------------------------------------------------------------------------------
  
  do n = 1, bufferSize
      buffer(n) = 0.0d0
  end do
  
  count = size(this%FlowsJA)
  do n = 1, count
      buffer(n) = this%FlowsJA(n)
  end do
  
end subroutine pr_FillFaceFlowsBuffer
  
subroutine pr_FillCellFlowsBuffer(this,cellNumber,buffer,bufferSize,count)
!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
  integer,intent(in) :: cellNumber,bufferSize
  doubleprecision,intent(inout),dimension(bufferSize) :: buffer
  integer,intent(inout) :: count
  integer :: n,offset
!---------------------------------------------------------------------------------------------------------------
  
  do n = 1, bufferSize
      buffer(n) = 0.0d0
  end do
  
  !offset = this%Grid%GetOffsetJa(cellNumber)
  !count = this%Grid%GetOffsetJa(cellNumber + 1) - offset
  offset = this%Grid%JaOffsets(cellNumber)
  count = this%Grid%JaOffsets(cellNumber + 1) - offset
  do n = 1, count
      buffer(n) = this%FlowsJA(offset + n)
  end do
  
end subroutine pr_FillCellFlowsBuffer
  
function pr_GetCurrentStressPeriod(this) result(stressPeriod)
!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
  integer :: stressPeriod
!---------------------------------------------------------------------------------------------------------------
 
  stressPeriod = this%CurrentStressPeriod
  
end function pr_GetCurrentStressPeriod

function pr_GetCurrentTimeStep(this) result(timeStep)
 !***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
 implicit none
  class(ParticleTrackingEngineType) :: this
  integer :: timeStep
!---------------------------------------------------------------------------------------------------------------
  
  timeStep = this%CurrentTimeStep
  
end function pr_GetCurrentTimeStep

subroutine pr_LoadTimeStep(this, stressPeriod, timeStep)
!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
  integer,intent(in) :: stressPeriod, timeStep
  integer :: firstRecord, lastRecord, n, m, firstNonBlank, lastNonBlank,        &
    trimmedLength
  integer :: spaceAssigned, status,cellCount, iface, index,                     &
    boundaryFlowsOffset, listItemBufferSize, cellNumber, layer
  type(BudgetRecordHeaderType) :: header
  character(len=16) :: textLabel
  doubleprecision :: top 
  real :: HDryTol, HDryDiff
!---------------------------------------------------------------------------------------------------------------
  
  call this%ClearTimeStepBudgetData()
  call this%BudgetReader%GetRecordHeaderRange(stressPeriod, timeStep, firstRecord, lastRecord)
  if(firstRecord .eq. 0) return

  cellCount = this%Grid%CellCount
  listItemBufferSize = size(this%ListItemBuffer)
  
  ! Set steady state = true, then change it if the budget file contains storage
  this%SteadyState = .true.
  
  ! Load heads for this time step
  call this%HeadReader%FillTimeStepHeadBuffer(stressPeriod, timeStep,           &
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
                  if((this%Heads(n) .le. this%Grid%Top(n)) .and.                                 &
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
                  if((this%Heads(n) .le. this%Grid%Top(n)) .and.                                 &
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
                            cellNumber = this%ArrayBufferInt(m)
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

subroutine pr_CheckForDefaultIface(this, textLabel, iface)
!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  use UtilMiscModule,only : utrimall
  implicit none
  class(ParticleTrackingEngineType) :: this
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

subroutine pr_ClearTimeStepBudgetData(this)
!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
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
  
subroutine pr_FillCellBuffer(this, cellNumber, cellBuffer)
!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType) :: this
  integer,intent(in) :: cellNumber
  type(ModpathCellDataType),intent(inout) :: cellBuffer
  doubleprecision,dimension(6) :: boundaryFlows
  integer :: n, layer, boundaryFlowsOffset, gridType, cellType
!---------------------------------------------------------------------------------------------------------------
  
  boundaryFlowsOffset = 6 * (cellNumber - 1)
  do n = 1, 6
      boundaryFlows(n) = this%BoundaryFlows(boundaryFlowsOffset + n)
  end do
  
  layer = this%Grid%GetLayer(cellNumber)
  
  gridType = this%Grid%GridType
  cellType = this%Grid%CellType(cellNumber)
  select case (gridType)
      case (1)
          ! Set cell buffer data for a structured grid
          call cellBuffer%SetDataStructured(cellNumber,this%Grid%CellCount,     &
            this%Grid,this%IBound,this%IBoundTS,                                &
            this%Porosity(cellNumber),this%Retardation(cellNumber),             & 
            this%StorageFlows(cellNumber),this%SourceFlows(cellNumber),         &
            this%SinkFlows(cellNumber), this%FlowsRightFace,                    &
            this%FlowsFrontFace, this%FlowsLowerFace, boundaryFlows,            & 
            this%Heads(cellNumber), cellType,                                   &
            this%Zones(cellNumber))
      case (2)
          ! Set cell buffer data for a MODFLOW-USG unstructured grid
          call cellBuffer%SetDataUnstructured(cellNumber,this%Grid%CellCount,   &
            this%Grid%JaCount,this%Grid,                                        &
            this%IBound,this%IBoundTS,                                          &
            this%Porosity(cellNumber),this%Retardation(cellNumber),             & 
            this%StorageFlows(cellNumber),this%SourceFlows(cellNumber),         &
            this%SinkFlows(cellNumber), this%FlowsJA, boundaryFlows,            & 
            this%Heads(cellNumber), cellType,                                   &
            this%Zones(cellNumber))
          ! Compute internal sub-cell face flows for cells with multiple sub-cells
          if(cellBuffer%GetSubCellCount() .gt. 1) then
              call cellBuffer%ComputeSubCellFlows()
          end if
      case (3)
          ! Set cell buffer data for a MODFLOW-6 structured grid (DIS)
          call cellBuffer%SetDataUnstructured(cellNumber,this%Grid%CellCount,   &
            this%Grid%JaCount,this%Grid,                                        &
            this%IBound,this%IBoundTS,                                          &
            this%Porosity(cellNumber),this%Retardation(cellNumber),             & 
            this%StorageFlows(cellNumber),this%SourceFlows(cellNumber),         &
            this%SinkFlows(cellNumber), this%FlowsJA, boundaryFlows,            & 
            this%Heads(cellNumber), cellType,                                   &
            this%Zones(cellNumber))
      case (4)
          ! Set cell buffer data for a MODFLOW-6 unstructured grid (DISV)
          call cellBuffer%SetDataUnstructured(cellNumber,this%Grid%CellCount,   &
            this%Grid%JaCount,this%Grid,                                        &
            this%IBound,this%IBoundTS,                                          &
            this%Porosity(cellNumber),this%Retardation(cellNumber),             & 
            this%StorageFlows(cellNumber),this%SourceFlows(cellNumber),         &
            this%SinkFlows(cellNumber), this%FlowsJA, boundaryFlows,            & 
            this%Heads(cellNumber), cellType,                                   &
            this%Zones(cellNumber))
           ! Compute internal sub-cell face flows for cells with multiple sub-cells
          if(cellBuffer%GetSubCellCount() .gt. 1) then
              call cellBuffer%ComputeSubCellFlows()
          end if
         
      case default
      ! Write error message and stop
  end select
  
end subroutine pr_FillCellBuffer

subroutine pr_Reset(this)
!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
!---------------------------------------------------------------------------------------------------------------
  class(ParticleTrackingEngineType) :: this
   
!   this%ReferenceTime = 0.0d0
!   this%StoppingTime = 0.0d0
   this%CurrentStressPeriod = 0
   this%CurrentTimeStep = 0
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

subroutine pr_TrackPath(this, trackPathResult, traceModeOn, traceModeUnit,      &
  group, particleID, seqNumber, location, maximumTrackingTime, timeseriesPoints,&
  timeseriesPointCount)
!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
  implicit none
  class(ParticleTrackingEngineType),target :: this
  type(TrackPathResultType),target,intent(out) :: trackPathResult
  type(ParticleLocationType),intent(in) :: location
  integer,intent(in) :: group, particleID, seqNumber, timeseriesPointCount,     &
    traceModeUnit
  logical,intent(in) :: traceModeOn
  type(ParticleLocationType) :: loc
  type(ParticleCoordinateType) :: pCoord
  type(ModpathCellDataType),pointer :: cellData
  type(TrackCellResultType),pointer :: tcResult
  doubleprecision,intent(in) :: maximumTrackingTime
  doubleprecision,dimension(timeseriesPointCount),intent(in) :: timeseriesPoints
  doubleprecision :: stopTime, fromLocalX, fromLocalY, fromLocalZ, globalX,     &
    globalY, globalZ
  integer :: timeIndex, n, count, nextCell
  logical :: continueLoop, isTimeSeriesPoint, isMaximumTime


  type(ModpathCellDataType), dimension( 2, 18 ) :: neighborCellData
  !type(ModpathCellDataType), dimension(18) :: neighborCellData

!---------------------------------------------------------------------------------------------------------------
  
  ! Reset trackPathResult and initialize particleID
  call trackPathResult%Reset()
  trackPathResult%ParticleID = particleID
  trackPathResult%Group = group
  trackPathResult%SequenceNumber = seqNumber
  
  ! Reset LocBuffP and LocBuffTS and initialize location data
  call this%LocBuffP%Clear()
  call this%LocBuffTS%Clear()
  call loc%SetData(location)
  call this%LocBuffP%AddItem(loc)
  
  ! Initialize loc
  call loc%SetData(location)
  
  ! Initialize TrackCell
  this%TrackCell%SteadyState = this%SteadyState
  this%TrackCell%TrackingOptions = this%TrackingOptions
  call this%FillCellBuffer(loc%CellNumber,  this%TrackCell%CellData)


  print *, '** TrackPath: Will fill neighbor cell data'

  ! RWPT
  if (this%TrackingOptions%RandomWalkParticleTracking) then
      !call this%FillNeighborSubCellData( this%NeighborSubCellData ) 
      ! DEV
      !call this%FillNeighborCellData( this%NeighborCellData )
      call this%FillNeighborCellData( neighborCellData )

  end if


  continueLoop = .true.
  isTimeSeriesPoint = .false.
  isMaximumTime = .false.
  
  do while(continueLoop)
      ! Check to see if the particle has moved to another cell. If so, load the new cell data
      if(loc%CellNumber .ne. this%TrackCell%CellData%CellNumber) then
          call this%FillCellBuffer(loc%CellNumber, this%TrackCell%CellData)
          ! RWPT
          if (this%TrackingOptions%RandomWalkParticleTracking) then
              !call this%FillNeighborSubCellData( this%NeighborSubCellData )
              ! DEV
              !call this%FillNeighborCellData( this%NeighborCellData )
              call this%FillNeighborCellData( neighborCellData )
          end if
      end if
      
      ! Find the next stopping time value (tmax), then track the particle through the cell starting at location loc.
      timeIndex = -1
      if(timeseriesPointCount .gt. 0) then
          timeIndex = pr_FindTimeIndex(timeseriesPoints, loc%TrackingTime,      &
            maximumTrackingTime, timeseriesPointCount)
      end if
      stopTime = maximumTrackingTime
      isTimeSeriesPoint = .false.
      if(timeIndex .ne. -1) then
          stopTime = timeseriesPoints(timeIndex)
          isTimeSeriesPoint = .true.
      end if
      isMaximumTime = .false.
      if(stopTime .eq. maximumTrackingTime) isMaximumTime = .true.


      ! RWPT
      ! Start with the particle loacion loc and track it through the cell until it reaches
      ! an exit face or the tracking time reaches the value specified by stopTime
      if ( .not. this%TrackingOptions%RandomWalkParticleTracking ) then 
          call this%TrackCell%ExecuteTracking(loc, stopTime, this%TrackCellResult)
      else
          !call this%TrackCell%ExecuteRandomWalkParticleTracking(loc, stopTime, this%TrackCellResult, this%NeighborSubCellData)
          ! DEV
          call this%TrackCell%ExecuteRandomWalkParticleTracking(loc, stopTime, this%TrackCellResult,&
              this%NeighborSubCellData, neighborCellData)

          call exit(0)

          !call this%TrackCell%ExecuteRandomWalkParticleTracking(loc, stopTime, this%TrackCellResult, &
          !                                             this%NeighborSubCellData, this%NeighborCellData)
          !call this%TrackCell%ExecuteRandomWalkParticleTracking(loc, stopTime, this%TrackCellResult, &
          !                                             this%NeighborSubCellData, neighborCellData)
      end if

      ! Check the status flag of the result to find out what to do next
      if(this%TrackCellResult%Status .eq. this%TrackCellResult%Status_Undefined()) then
          continueLoop = .false.
          trackPathResult%Status = this%TrackCellResult%Status
      else if(this%TrackCellResult%Status .eq. this%TrackCellResult%Status_ExitAtCellFace()) then
          count = this%TrackCellResult%TrackingPoints%GetItemCount()
          if(count .gt. 1) then
              do n = 2, count
                  call this%LocBuffP%AddItem(this%TrackCellResult%TrackingPoints%Items(n))
              end do 
          end if
          
          ! If NextCellNumber is > 0, it means the particle has moved to another cell. 
          ! If so, convert loc from the current cell coordinates to the equivalent location in the new cell.
          nextCell = this%TrackCellResult%NextCellNumber
          if(nextCell .gt. 0) then
              if(this%IBoundTS(nextCell) .ne. 0) then
                  ! The next cell is active
                  fromLocalX = this%TrackCellResult%TrackingPoints%Items(count)%LocalX
                  fromLocalY = this%TrackCellResult%TrackingPoints%Items(count)%LocalY
                  fromLocalZ = this%TrackCellResult%TrackingPoints%Items(count)%LocalZ         
                  call this%Grid%ConvertFromNeighbor(                           &
                    this%TrackCellResult%NextCellNumber,                        &
                    this%TrackCellResult%CellNumber, fromLocalX, fromLocalY,    &
                    fromLocalZ, loc)
                  loc%TrackingTime = this%TrackCellResult%TrackingPoints%Items(count)%TrackingTime
              else
                  ! If next cell is inactive, it implies that a boundary face has been reached. 
                  ! Set status and return.
                  continueLoop = .false.
                  trackPathResult%Status = trackPathResult%Status_ReachedBoundaryFace()        
              end if
          else
              ! If next cell number = 0, the boundary of the grid has been reached. 
              ! Set status and return.
              continueLoop = .false.
              trackPathResult%Status = trackPathResult%Status_ReachedBoundaryFace()
          end if
          
      else if(this%TrackCellResult%Status .eq. this%TrackCellResult%Status_ReachedStoppingTime()) then
          count = this%TrackCellResult%TrackingPoints%GetItemCount()
          if(count .gt. 1) then
              do n = 2, count
                  call this%LocBuffP%AddItem(this%TrackCellResult%TrackingPoints%Items(n))
              end do 
              if(isTimeSeriesPoint) then
                  call this%LocBuffTS%AddItem(this%TrackCellResult%TrackingPoints%Items(count))
              end if
          end if
          
          call loc%SetData(this%TrackCellResult%TrackingPoints%Items(count))
          if(isMaximumTime) then
              continueLoop = .false.
              trackPathResult%Status = trackPathResult%Status_ReachedStoppingTime()
          end if
          
      else
          continueLoop = .false.
          if(this%TrackCellResult%Status .eq. this%TrackCellResult%Status_NoExitPossible()) then
              trackPathResult%Status = this%TrackCellResult%Status_NoExitPossible()         
          else if(this%TrackCellResult%Status .eq. this%TrackCellResult%Status_StopZoneCell()) then
              trackPathResult%Status = this%TrackCellResult%Status_StopZoneCell()                  
          else if(this%TrackCellResult%Status .eq. this%TrackCellResult%Status_StopAtWeakSink()) then
              trackPathResult%Status = this%TrackCellResult%Status_StopAtWeakSink()                   
          else if(this%TrackCellResult%Status .eq. this%TrackCellResult%Status_StopAtWeakSource()) then
              trackPathResult%Status = this%TrackCellResult%Status_StopAtWeakSource()
          else if(this%TrackCellResult%Status .eq. this%TrackCellResult%Status_InactiveCell()) then
              trackPathResult%Status = this%TrackCellResult%Status_InactiveCell()      
          else
              trackPathResult%Status = this%TrackCellResult%Status_Undefined()
          end if
          
          ! If the trackPathResult status is anything except Undefined, add the last tracking point to
          ! the trackPathResult tracking points
          if(trackPathResult%Status .ne. this%TrackCellResult%Status_Undefined()) then
              call this%LocBuffP%AddItem(this%TrackCellResult%TrackingPoints%Items(1))              
          end if
          
      end if
  
      ! Write trace mode data if the trace mode is on for this particle
      if(traceModeOn) then
         !$omp critical (tracedata)
         call WriteTraceData(traceModeUnit, this%TrackCell,                     &
           this%TrackCellResult, this%GetCurrentStressPeriod(),                 &
           this%GetCurrentTimeStep())
         !$omp end critical (tracedata)
      end if
      
      ! If continueLoop is still set to true, go through the loop again. If set to false, exit the loop now.
  end do
  
  ! Generate global coordinates and finish initializing the result data
  count = this%LocBuffP%GetItemCount()
  do n = 1, count
      pCoord%CellNumber = this%LocBuffP%Items(n)%CellNumber
      pCoord%Layer = this%LocBuffP%Items(n)%Layer
      pCoord%LocalX = this%LocBuffP%Items(n)%LocalX
      pCoord%LocalY = this%LocBuffP%Items(n)%LocalY
      pCoord%LocalZ = this%LocBuffP%Items(n)%LocalZ
      pCoord%TrackingTime = this%LocBuffP%Items(n)%TrackingTime
      call this%Grid%ConvertToModelXYZ(pCoord%CellNumber, pCoord%LocalX,       &
        pCoord%LocalY, pCoord%LocalZ, pCoord%GlobalX, pCoord%GlobalY,           &
        pCoord%GlobalZ)
      call trackPathResult%ParticlePath%Pathline%AddItem(pCoord)
  end do
  
  do n = 1, this%LocBuffTS%GetItemCount()
      pCoord%CellNumber = this%LocBuffTS%Items(n)%CellNumber
      pCoord%Layer = this%LocBuffTS%Items(n)%Layer
      pCoord%LocalX = this%LocBuffTS%Items(n)%LocalX
      pCoord%LocalY = this%LocBuffTS%Items(n)%LocalY
      pCoord%LocalZ = this%LocBuffTS%Items(n)%LocalZ
      pCoord%TrackingTime = this%LocBuffTS%Items(n)%TrackingTime
      call this%Grid%ConvertToModelXYZ(pCoord%CellNumber, pCoord%LocalX,       &
        pCoord%LocalY, pCoord%LocalZ, pCoord%GlobalX, pCoord%GlobalY,           &
        pCoord%GlobalZ)
      call trackPathResult%ParticlePath%Timeseries%AddItem(pCoord)
  end do
  
end subroutine pr_TrackPath


! RWPT
! DEPRECATION WARNING
subroutine pr_FillNeighborSubCellData( this, neighborSubCellData )
    !------------------------------------------------------------------
    ! 
    ! Populates array of neighborSubCellData, with
    ! initialized ModpathSubCellData objects required for computation
    ! of corner velocities for the given TrackCell
    !
    ! Format of ordering into the neighborSubCellData array is (scd: subCellData):
    !
    !   ( scd1, scd13, scd14, sc2, scd23, scd24, 
    !     scd3, scd35, scd36, ... scd6, scd61, scd62)
    !
    ! Note:
    !   This code should be executed after each update of the current TrackCell.
    !   It would be desirable something to organize indexes involved into computation 
    !   of corner velocities, to avoid hardcoding or cumbersome access
    !------------------------------------------------------------------
    !
    ! Specifications
    !------------------------------------------------------------------
    implicit none
    class(ParticleTrackingEngineType),target :: this
    type(ModpathSubCellDataType), dimension(18) :: neighborSubCellData
    type(ModpathCellDataType) :: neighborCellData, neighborConnectionCellData
    integer :: n, m, subCellCounter
    !------------------------------------------------------------------


    ! RWPT
    ! Once cellBuffer has the data, 
    ! populate the array with neighbor cells
    ! Assuming only one connection per cell, that is
    ! fully structured grid.
    ! Only done for cells with a connection,
    ! stills remains to be verified how to handle the case
    ! without connection
    subCellCounter = 0
    do n = 1, 6
        subCellCounter = subCellCounter + 1
  
        ! If connection exists
        if ( this%TrackCell%CellData%GetFaceConnection(n,1) .gt. 0) then
  
            ! Fill the data buffer
            call this%FillCellBuffer( this%TrackCell%CellData%GetFaceConnection(n,1) , neighborCellData )
  
            ! Fill the subCell data buffer
            call neighborCellData%FillSubCellDataBuffer( neighborSubCellData( subCellCounter ), &
               1, 1, this%TrackingOptions%BackwardTracking ) 
  
            ! Extract surrounding cells from connection
            ! to complete information required for 
            ! computation of cell corner velocities
            select case (n)
                case (1)
                    do m = 3,4
                        subCellCounter = subCellCounter + 1
                        if ( neighborCellData%GetFaceConnection(m,1) .gt. 0 ) then
                            call this%FillCellBuffer( neighborCellData%GetFaceConnection(m,1), neighborConnectionCellData )
                            call neighborConnectionCellData%FillSubCellDataBuffer( neighborSubCellData( subCellCounter ), &
                                1, 1, this%TrackingOptions%BackwardTracking )
                        else
                           call neighborSubCellData( subCellCounter )%Reset() 
                        end if
                    end do
                case (2)
                    do m = 3,4
                        subCellCounter = subCellCounter + 1
                        if ( neighborCellData%GetFaceConnection(m,1) .gt. 0 ) then
                            call this%FillCellBuffer( neighborCellData%GetFaceConnection(m,1), neighborConnectionCellData )
                            call neighborConnectionCellData%FillSubCellDataBuffer( neighborSubCellData( subCellCounter ), &
                                1, 1, this%TrackingOptions%BackwardTracking )
                        else
                           call neighborSubCellData( subCellCounter )%Reset() 
                        end if
                    end do
                case (3)
                    do m = 5,6
                        subCellCounter = subCellCounter + 1
                        if ( neighborCellData%GetFaceConnection(m,1) .gt. 0 ) then
                            call this%FillCellBuffer( neighborCellData%GetFaceConnection(m,1), neighborConnectionCellData )
                            call neighborConnectionCellData%FillSubCellDataBuffer( neighborSubCellData( subCellCounter ), &
                                1, 1, this%TrackingOptions%BackwardTracking )
                        else
                           call neighborSubCellData( subCellCounter )%Reset() 
                        end if
                    end do
                case (4)
                    do m = 5,6
                        subCellCounter = subCellCounter + 1
                        if ( neighborCellData%GetFaceConnection(m,1) .gt. 0 ) then
                            call this%FillCellBuffer( neighborCellData%GetFaceConnection(m,1), neighborConnectionCellData )
                            call neighborConnectionCellData%FillSubCellDataBuffer( neighborSubCellData( subCellCounter ), &
                                1, 1, this%TrackingOptions%BackwardTracking )
                        else
                           call neighborSubCellData( subCellCounter )%Reset() 
                        end if
                    end do
                case (5)
                    do m = 1,2
                        subCellCounter = subCellCounter + 1
                        if ( neighborCellData%GetFaceConnection(m,1) .gt. 0 ) then
                            call this%FillCellBuffer( neighborCellData%GetFaceConnection(m,1), neighborConnectionCellData )
                            call neighborConnectionCellData%FillSubCellDataBuffer( neighborSubCellData( subCellCounter ), &
                                1, 1, this%TrackingOptions%BackwardTracking )
                        else
                           call neighborSubCellData( subCellCounter )%Reset() 
                        end if
                    end do
                case (6)
                    do m = 1,2
                        subCellCounter = subCellCounter + 1
                        if ( neighborCellData%GetFaceConnection(m,1) .gt. 0 ) then
                            call this%FillCellBuffer( neighborCellData%GetFaceConnection(m,1), neighborConnectionCellData )
                            call neighborConnectionCellData%FillSubCellDataBuffer( neighborSubCellData( subCellCounter ), &
                                1, 1, this%TrackingOptions%BackwardTracking )
                        else
                           call neighborSubCellData( subCellCounter )%Reset()
                        end if
                    end do
            end select
        else
            ! If no connection a cell face n, 
            ! suppose no connections initialize velocities 
            ! to be exactly zero and move counter
            call neighborSubCellData( subCellCounter )%Reset()
            call neighborSubCellData( subCellCounter + 1 )%Reset()
            call neighborSubCellData( subCellCounter + 2 )%Reset()
            subCellCounter = subCellCounter + 2
        endif
    end do


end subroutine pr_FillNeighborSubCellData


! RWPT
subroutine pr_FillNeighborCellData( this, neighborCellData )
    !------------------------------------------------------------------
    !------------------------------------------------------------------
    !
    ! Specifications
    !------------------------------------------------------------------
    implicit none
    class(ParticleTrackingEngineType),target :: this
    type(ModpathCellDataType), dimension(2, 18) :: neighborCellData
    integer :: n, m, id, cellCounter, firstNeighborFaceNumber
    integer, dimension(2)   :: orthogonalFaceNumbers, connectedSubCellIndexes
    integer, dimension(2,2) :: subCellIndexes
    logical :: forceCellRefinement = .false.
    !------------------------------------------------------------------

    ! Reset all buffers
    do m = 1, 18 
        call neighborCellData( 1, m )%Reset() 
        call neighborCellData( 2, m )%Reset() 
    end do 



    ! Verify what is going on with the current TrackCell
    ! Is it refined ?
    print *, '*******************************************************************'
    if ( this%TrackCell%CellData%GetSubCellCount() .gt. 1 ) then 
        print *, '** ParticleTrackingEngine.FillNeighborCellData: TrackCell is refined ', &
        this%TrackCell%CellData%CellNumber, this%TrackCell%CellData%GetSubCellCount(), ' subcells.'
        forceCellRefinement = .true.
    else
        print *, '** ParticleTrackingEngine.FillNeighborCellData: TrackCell NOT REFINED ', &
        this%TrackCell%CellData%CellNumber, this%TrackCell%CellData%GetSubCellCount(), ' subcells.'
    end if 


    ! Loop through cell faces and fill the neighborhood
    cellCounter = 0
    do n = 1, 6

        ! Each cell face gives three neighbor cells. 
        ! Itself, and two other orthogonal connections
        cellCounter = ( n - 1 )*3 + 1

        ! From the cell connected at face n,
        ! create TrackCells from two connections 
        ! starting from firstNeighborFaceNumber  
        if ( n .le. 2 ) then
            firstNeighborFaceNumber = 3
        else if ( n .le. 4 ) then
            firstNeighborFaceNumber = 5
        else
            firstNeighborFaceNumber = 1
        end if

        ! Thus, broadly speaking, what this loop 
        ! should do is fill buffers corresponding 
        ! to: 
        !   - cellCounter    : direct connection 
        !   - cellCounter + 1: connected neighbor 1, (neighborFaceNumber)
        !   - cellCounter + 2: connected neighbor 2, (neighborFaceNumber+1)
        call pr_FillNeighborCellsSubBuffer( this, this%TrackCell%CellData, n, & 
            firstNeighborFaceNumber, neighborCellData, cellCounter, forceCellRefinement ) 

    end do ! End loop through cell faces

    !! DEBUG/DEV
    ! Print a report
    print *, '********************************************************************************************'
    print *, '** FillNeighborCellData: neighbors information for TrackCell', this%TrackCell%CellData%CellNumber

    do n = 1, 18
        print *, '** FillNeighborCellData: neighbor cell buffer ', n
        do m = 1, 2
            if ( neighborCellData(m, n)%CellNumber .gt. 0 ) then 
                print *, '     -- Face Connection  ', m, ' with cell ', neighborCellData(m, n)%CellNumber
            else if ( neighborCellData(m, n)%fromSubSubCell ) then 
                print *, '     -- Face Connection  ', m, ' with cell filled from subcell ', &
                    neighborCellData(m, n)%parentSubRow, neighborCellData(m, n)%parentSubColumn, & 
                    'from grandParentCellNumber ', neighborCellData(m, n)%parentCellNumber
            else if ( neighborCellData(m, n)%fromSubCell ) then 
                print *, '     -- Face Connection  ', m, ' with cell filled from subcell ', &
                    neighborCellData(m, n)%parentSubRow, neighborCellData(m, n)%parentSubColumn, & 
                    'from parentCellNumber ', neighborCellData(m, n)%parentCellNumber
            end if 
        end do 
    end do

    !call exit(0)

end subroutine pr_FillNeighborCellData


!RWPT
subroutine pr_FillNeighborCellsSubBuffer( this, &
        centerCellDataBuffer, faceNumber, firstNeighborFaceNumber, & 
        neighborCellsDataBuffer, directConnectionBufferIndex, forceCellRefinement )
    !----------------------------------------------------------------------------------
    !
    ! NEEDS REVIEW
    !
    ! Objective of this function is to fill portion of neighborCellsDataBuffer,
    ! related to cells connected to centerCellDataBuffer through faceNumber.
    ! In order to build a complete neighborhood for interpolation purposes, 
    ! it is necessary to request information about the cell connected through faceNumber
    ! and then from this cell build two additional neighbors, for each faceNumber. 
    ! The latter yields a method to get information from indirect connections 
    ! like those cells diagonally connected to the centerCellDataBuffer.
    !
    ! Method assumes that centerCellDataBuffer will always be a real grid cell, 
    ! such that the function FillCellBuffer can be employed without issues and from 
    ! there it is possible to extract reliable connection information.
    !
    ! Function will request information about connection state through faceNumber.
    ! This leads to the following possible scenarios:
    !
    !   - If faceNumber is horizontal (x-y plane)
    !
    !       - centerCellDataBuffer has multiple connections through faceNumber 
    !           
    !       - centerCellDataBuffer has a single connection through faceNumber
    !
    !           - Is this connection shared with another cell ? (connected to a bigger cell)
    !           - Is this connection unique (same size cells)
    !
    !   - If faceNumber is vertical (z axis)
    !       
    !       - Fill direct connection buffer with connection through faceNumber
    !       - Fill further connections (horizontal) 
    !
    ! Assumptions:
    !
    !   - No vertical refinement 
    !   - Smoothed unstructured grid
    !
    !----------------------------------------------------------------------------------
    class(ParticleTrackingEngineType) :: this
    ! input
    type(ModpathCellDataType), intent(in) :: centerCellDataBuffer
    integer, intent(in) :: faceNumber, firstNeighborFaceNumber, directConnectionBufferIndex
    logical, intent(in) :: forceCellRefinement
    ! output
    type(ModpathCellDataType), dimension(2, 18), intent(inout)  :: neighborCellsDataBuffer
    ! local
    type(ModpathCellDataType) :: auxCellDataBuffer
    integer, dimension(2) :: faceConnectionIndexes
    integer :: sameParentSubRow
    integer :: sameParentSubColumn
    integer :: sameParentCellBufferIndex
    integer :: notSameParentCellBufferIndex
    integer :: notSameParentCellFaceConnection
    integer :: n, m 
    integer :: directionId
    !---------------------------------------------------------------------


    ! Handle the case of a z-axis faceNumber
    if ( faceNumber .ge. 5 ) then
        ! IT SHOULD NOT DO THIS IF THE SIMULATION IS 2D
        ! DEBUG/DEV 
        print *, '*** FillNeighborCellsSubBuffer: will initialize from a vertical faceNumber', faceNumber
        
        ! By design, indirect connections requested 
        ! from a direct connection initialized from 
        ! a vertical faceNumber are done through 
        ! horizontal faces

        ! Fill direct connection buffer
        call pr_FillCellFromRealData( this, centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ), &
                      neighborCellsDataBuffer( 1, directConnectionBufferIndex ), forceCellRefinement )
        neighborCellsDataBuffer( 1, directConnectionBufferIndex )%requestedFromDirection = 3

        ! Fill indirect connection buffers
        call pr_FillNeighborCellsConnectionFromHorizontalFace( this,                          &
                                   neighborCellsDataBuffer( 1, directConnectionBufferIndex ), &
                                                                     firstNeighborFaceNumber, &
                               neighborCellsDataBuffer( :, directConnectionBufferIndex + 1 ), &
                                                                          forceCellRefinement )

        call pr_FillNeighborCellsConnectionFromHorizontalFace( this,                              &
                                       neighborCellsDataBuffer( 1, directConnectionBufferIndex ), &
                                                                     firstNeighborFaceNumber + 1, &
                                   neighborCellsDataBuffer( :, directConnectionBufferIndex + 2 ), &
                                                                              forceCellRefinement )
        ! Done
        return

    end if 


    ! DEBUG/DEV
    print *, '*** FillNeighborCellsSubBuffer: will initialize from a horizontal faceNumber', faceNumber


    ! If horizontal faceNumber
    
    ! Connected to one or more cells ?
    if ( centerCellDataBuffer%SubFaceCounts( faceNumber ) .gt. 1 ) then
        ! DEBUG/DEV
        print *, '*** FillNeighborCellsSubBuffer: found more than one connection'

        ! Fill direct connection buffer 
        call pr_FillNeighborCellsConnectionFromHorizontalFace( this, centerCellDataBuffer, &
                    faceNumber, neighborCellsDataBuffer( :, directConnectionBufferIndex ), &
                                                                       forceCellRefinement )

        ! Fill indirect connection buffers
        if ( firstNeighborFaceNumber .ge. 5 ) then 
            ! IT SHOULD NOT DO THIS FOR 2D (SINGLE LAYER)
            ! DEBUG/DEV 
            print *, '*** FillNeighborCellsSubBuffer: will fill remaining from vertical faces '

            ! Under the assumption of same horizontal 
            ! refinement distribution, vertical connections
            ! are filled directly

            ! First indirect connection buffer
            do m =1,2
                call pr_FillCellFromRealData( this,    &
                    neighborCellsDataBuffer(           &
                        m, directConnectionBufferIndex &
                    )%GetFaceConnection( firstNeighborFaceNumber, 1 ), &
                    neighborCellsDataBuffer( m, directConnectionBufferIndex + 1 ), &
                                                               forceCellRefinement )
                neighborCellsDataBuffer( m, directConnectionBufferIndex + 1 )%requestedFromDirection = 3
            end do 

            ! Second indirect connection buffer
            do m =1,2
                call pr_FillCellFromRealData( this,    &
                    neighborCellsDataBuffer(           &
                        m, directConnectionBufferIndex &
                    )%GetFaceConnection( firstNeighborFaceNumber + 1, 1 ), &
                    neighborCellsDataBuffer( m, directConnectionBufferIndex + 2 ), &
                                                               forceCellRefinement )
                neighborCellsDataBuffer( m, directConnectionBufferIndex + 2 )%requestedFromDirection = 3
            end do 

            ! Done
            return

        else
            ! DEBUG/DEV 
            print *, '*** FillNeighborCellsSubBuffer: will fill remaining from horizontal faces '

            ! If the others requested neighbor cells are through 
            ! horizontal faces, change centerCellDataBuffer to the 
            ! recently filled buffer and request info 

            ! Verify direction of requested faces, 
            ! to identify in which order access the 
            ! previously filled buffers
            if ( firstNeighborFaceNumber .le. 2 ) then
                ! In the x-axis, order in which neighbor faces are requested
                ! coincides with their indexation in buffer (faceConnections) 
                faceConnectionIndexes = (/1,2/)
            else
                faceConnectionIndexes = (/2,1/)
            end if
           

            call pr_FillNeighborCellsConnectionFromHorizontalFace( this,                          &
                neighborCellsDataBuffer( faceConnectionIndexes(1), directConnectionBufferIndex ), &
                                                                         firstNeighborFaceNumber, &
                                   neighborCellsDataBuffer( :, directConnectionBufferIndex + 1 ), &
                                                                              forceCellRefinement )

            call pr_FillNeighborCellsConnectionFromHorizontalFace( this,                              &
                    neighborCellsDataBuffer( faceConnectionIndexes(2), directConnectionBufferIndex ), &
                                                                         firstNeighborFaceNumber + 1, &
                                       neighborCellsDataBuffer( :, directConnectionBufferIndex + 2 ), &
                                                                                  forceCellRefinement )

            ! Done 
            return


        end if


    else if ( centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ) .gt. 0 ) then
        print *, '*** FillNeighborCellsSubBuffer: there is one connection', &
            centerCellDataBuffer%GetFaceConnection( faceNumber, 1 )


        ! Fill direct connection buffer 
        call pr_FillNeighborCellsConnectionFromHorizontalFace( this, centerCellDataBuffer, &
                    faceNumber, neighborCellsDataBuffer( :, directConnectionBufferIndex ), &
                                                                       forceCellRefinement )

        ! If the buffer was filled with the one of the subcells
        ! from a bigger cell
        if ( neighborCellsDataBuffer( 1, directConnectionBufferIndex )%fromSubCell ) then 
            print *, '*** FillNeighborCellsSubBuffer: buffer was filled from a BIG cell ', &
            neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentCellNumber, &
            neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentSubRow, &
            neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentSubColumn


            ! Fill indirect connection buffers
            if ( firstNeighborFaceNumber .ge. 5 ) then 
                print *, '*** FillNeighborCellsSubBuffer: will fill remaining from vertical faces '

                ! Fill auxCellDataBuffer with vertical connections
                ! and from there extract subcell
                call pr_FillCellFromRealData( this, &
                    this%Grid%GetFaceConnection( &
                        neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentCellNumber, &
                                            firstNeighborFaceNumber, 1 ), auxCellDataBuffer, .true. )

                if ( auxCellDataBuffer%CellNumber .gt. 0 ) then
                    ! Fill connection from subcell
                    call pr_FillCellBufferFromSubCell( this, auxCellDataBuffer, &
                        neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentSubRow,    &
                        neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentSubColumn, & 
                                     neighborCellsDataBuffer( 1, directConnectionBufferIndex + 1), &
                                                                               forceCellRefinement )
                    neighborCellsDataBuffer( 1, directConnectionBufferIndex + 1)%requestedFromDirection = 3
                end if

                call pr_FillCellFromRealData( this, &
                    this%Grid%GetFaceConnection( &
                        neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentCellNumber, &
                                        firstNeighborFaceNumber + 1, 1 ), auxCellDataBuffer, .true. )

                if ( auxCellDataBuffer%CellNumber .gt. 0 ) then
                    ! Fill connection from subcell
                    call pr_FillCellBufferFromSubCell( this, auxCellDataBuffer, &
                        neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentSubRow,    &
                        neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentSubColumn, & 
                                     neighborCellsDataBuffer( 1, directConnectionBufferIndex + 2), &
                                                                               forceCellRefinement )
                    neighborCellsDataBuffer( 1, directConnectionBufferIndex + 2)%requestedFromDirection = 3
                end if

            else 
                print *, '*** FillNeighborCellsSubBuffer: will fill remaining from horizontal faces '

                ! Requires  a parent cellDataBuffer
                call pr_FillCellFromRealData( this, &
                    neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentCellNumber, &
                                                                      auxCellDataBuffer, .true. )
                auxCellDataBuffer%isParentCell    = .true. 
                auxCellDataBuffer%parentSubRow    = neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentSubRow
                auxCellDataBuffer%parentSubColumn = neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentSubColumn

                ! Verify direction of requested faces
                if ( firstNeighborFaceNumber .ge. 3 ) then
                    print *, '*** FillNeighborCellsSubBuffer: Y AXIS FACES '
                    directionId = 2
                    if ( neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentSubRow .eq. 1 ) then
                        sameParentCellBufferIndex       = directConnectionBufferIndex + 1
                        notSameParentCellBufferIndex    = directConnectionBufferIndex + 2
                        notSameParentCellFaceConnection = firstNeighborFaceNumber + 1 
                        sameParentSubRow    = 2 
                        sameParentSubColumn = neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentSubColumn
                    else
                        sameParentCellBufferIndex       = directConnectionBufferIndex + 2
                        notSameParentCellBufferIndex    = directConnectionBufferIndex + 1
                        notSameParentCellFaceConnection = firstNeighborFaceNumber 
                        sameParentSubRow    = 1 
                        sameParentSubColumn = neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentSubColumn
                    end if 
                else
                    print *, '*** FillNeighborCellsSubBuffer: X AXIS FACES '
                    directionId = 1
                    ! This case does not happens under the current scheme. 
                    ! That is because x axis faces are requested after 
                    ! a vertical faceNumber. As the grid cell is assumed
                    ! to have same refinement features for all layers, 
                    ! the vertical connection is never filled from a 
                    ! a bigger cell, only from cells same size as the 
                    ! current/center cell.
                end if

                ! Fill connection, not same parent
                call pr_FillNeighborCellsConnectionFromHorizontalFace( this,    & 
                            auxCellDataBuffer, notSameParentCellFaceConnection, &
                    neighborCellsDataBuffer( :, notSameParentCellBufferIndex ), &
                                                            forceCellRefinement )


                ! Fill connection from subcell, same parent
                call pr_FillCellBufferFromSubCell( this, auxCellDataBuffer, &
                                                          sameParentSubRow, &
                                                       sameParentSubColumn, & 
                   neighborCellsDataBuffer( 1, sameParentCellBufferIndex ), &
                                                        forceCellRefinement )
                neighborCellsDataBuffer( 1, sameParentCellBufferIndex)%requestedFromDirection = directionId

            end if 


            ! Done
            return


        else


            print *, '*** FillNeighborCellsSubBuffer: buffer was filled from an equal size cell '

            ! Fill indirect connection buffers
            ! taking as center the already filled buffer 
            if ( firstNeighborFaceNumber .ge. 5 ) then 
                ! IT SHOULD NOT DO THIS FOR 2D (SINGLE LAYER)
                ! DEBUG/DEV 
                print *, '*** FillNeighborCellsSubBuffer: will fill remaining from vertical faces '
        
                ! This is solved by filling cell from real cell data
                ! requesting directly vertical faces as they have the same 
                ! size that center cell, due to assumption of same refinement distribution
                ! for all layers.

                call pr_FillCellFromRealData( this,    &
                    neighborCellsDataBuffer(           &
                        1, directConnectionBufferIndex &
                    )%GetFaceConnection( firstNeighborFaceNumber, 1 ), &
                    neighborCellsDataBuffer( 1, directConnectionBufferIndex + 1 ), &
                                                               forceCellRefinement )
                neighborCellsDataBuffer( 1, directConnectionBufferIndex + 1 )%requestedFromDirection = 3
                call pr_FillCellFromRealData( this, &
                    neighborCellsDataBuffer(           &
                        1, directConnectionBufferIndex &
                    )%GetFaceConnection( firstNeighborFaceNumber + 1, 1 ), &
                    neighborCellsDataBuffer( 1, directConnectionBufferIndex + 2 ), &
                                                               forceCellRefinement )
                neighborCellsDataBuffer( 1, directConnectionBufferIndex + 2 )%requestedFromDirection = 3
            else

                ! DEBUG/DEV
                print *, '*** FillNeighborCellsSubBuffer: will fill remaining from horizontal faces '

                call pr_FillNeighborCellsConnectionFromHorizontalFace( this,                          &
                                           neighborCellsDataBuffer( 1, directConnectionBufferIndex ), &
                                                                             firstNeighborFaceNumber, &
                                       neighborCellsDataBuffer( :, directConnectionBufferIndex + 1 ), &
                                                                                  forceCellRefinement )

                call pr_FillNeighborCellsConnectionFromHorizontalFace( this,                              &
                                               neighborCellsDataBuffer( 1, directConnectionBufferIndex ), &
                                                                             firstNeighborFaceNumber + 1, &
                                           neighborCellsDataBuffer( :, directConnectionBufferIndex + 2 ), &
                                                                                      forceCellRefinement )
            end if

            
            ! Done
            return


        end if 
                

    else
        ! No connection 

        ! Done
        return

    end if



end subroutine pr_FillNeighborCellsSubBuffer


! RWPT
subroutine pr_FillNeighborCellsConnectionFromHorizontalFace(this, centerCellDataBuffer, &
                                   faceNumber, neighborCellsBuffer, forceCellRefinement )
    !--------------------------------------------------------------------------------------
    ! 
    !--------------------------------------------------------------------------------------
    class(ParticleTrackingEngineType) :: this
    ! input
    type(ModpathCellDataType), intent(in)    :: centerCellDataBuffer
    integer, intent(in)                      :: faceNumber
    logical, intent(in)                      :: forceCellRefinement
    ! output
    type(ModpathCellDataType), dimension(2), intent(inout) :: neighborCellsBuffer
    ! local 
    type(ModpathCellDataType)                :: parentCellDataBuffer
    type(ModpathCellDataType)                :: auxCellDataBuffer
    integer, dimension(2)                    :: orthogonalFaceNumbers
    integer, dimension(2,2)                  :: subCellIndexes
    integer                                  :: directConnectionSubRow
    integer                                  :: subCellIndexesRow
    integer                                  :: directionId
    integer                                  :: subRow, subColumn
    integer                                  :: m
    !--------------------------------------------------------------------------------------


    ! If invalid center cell leave 
    if ( centerCellDataBuffer%CellNumber .le. 0 ) then
        ! Done
        return
    end if 

    
    ! Determine direction of faceNumber 
    if ( faceNumber .le. 2 ) then 
        directionId = 1 
    else if ( faceNumber .le. 4 ) then 
        directionId = 2 
    else 
        directionId = 3
    end if


    ! Verify horizontal connection state
    if ( centerCellDataBuffer%SubFaceCounts( faceNumber ) .gt. 1 ) then

        if ( .not. centerCellDataBuffer%isParentCell ) then 

            ! Loop over face connections and extract 
            ! their connected cell through same faceNumber
            do m = 1, centerCellDataBuffer%SubFaceCounts( faceNumber )
                call this%FillCellBuffer(&
                    centerCellDataBuffer%GetFaceConnection( faceNumber, m ), neighborCellsBuffer( m ) )
                neighborCellsBuffer(m)%requestedFromDirection = directionId
            end do

        else 

            ! Which face connection ?
            if ( directionId .eq. 2 ) then 
                m = centerCellDataBuffer%parentSubColumn  
            else if ( directionId .eq. 1 ) then
                ! It should not happen under the current neighborhood generation scheme
                print *, 'IT HAPPENED'
                m = centerCellDataBuffer%parentSubRow 
            end if

            call this%FillCellBuffer(&
                centerCellDataBuffer%GetFaceConnection( faceNumber, m ), neighborCellsBuffer( 1 ) )
            neighborCellsBuffer( 1 )%requestedFromDirection = directionId
            
            if ( forceCellRefinement ) then 
                ! Forces refinement and sub cell flows computation
                neighborCellsBuffer( 1 )%SubCellRowCount    = 2
                neighborCellsBuffer( 1 )%SubCellColumnCount = 2
                call neighborCellsBuffer( 1 )%ComputeSubCellFlows()
            end if 

        end if 


        ! Done
        return

    else if ( centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ) .gt. 0 ) then

        ! Only one connection, verify if 
        ! same size or bigger size cell

        ! Info for verification
        select case ( faceNumber )
            case (1)
                orthogonalFaceNumbers(1)   = 3 
                orthogonalFaceNumbers(2)   = 4
                subCellIndexes(1,:)        = (/1,2/) 
                subCellIndexes(2,:)        = (/2,2/) 
            case (2)
                orthogonalFaceNumbers(1)   = 3 
                orthogonalFaceNumbers(2)   = 4
                subCellIndexes(1,:)        = (/1,1/) 
                subCellIndexes(2,:)        = (/2,1/) 
            case (3)
                orthogonalFaceNumbers(1)   = 1 
                orthogonalFaceNumbers(2)   = 2
                subCellIndexes(1,:)        = (/1,2/) 
                subCellIndexes(2,:)        = (/1,1/) 
            case (4)
                orthogonalFaceNumbers(1)   = 1 
                orthogonalFaceNumbers(2)   = 2
                subCellIndexes(1,:)        = (/2,2/) 
                subCellIndexes(2,:)        = (/2,1/) 
        end select


        ! Detect connection state
        ! Verify response when no connection 
        if ( &
            ( centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ) .eq.                                 &
              this%Grid%GetFaceConnection(                                                                 & 
                  centerCellDataBuffer%GetFaceConnection( orthogonalFaceNumbers(1), 1 ), faceNumber, 1 ) ) &
            .or.                                                                                           &
            ( centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ) .eq.                                 &
              this%Grid%GetFaceConnection(                                                                 & 
                  centerCellDataBuffer%GetFaceConnection( orthogonalFaceNumbers(2), 1 ), faceNumber, 1 ) ) &
        ) then 

            ! If bigger cell
            ! Fill the parentCellDataBuffer
            ! and force refinement
            call pr_FillCellFromRealData( this, centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ), &
                                                                            parentCellDataBuffer, .true. )

            ! Location of current TrackCell relative to bigger cell
            ! defines indexes employed for filling buffers. 
            if ( &
                ( centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ) .eq.                               &
                  this%Grid%GetFaceConnection(                                                               & 
                    centerCellDataBuffer%GetFaceConnection( orthogonalFaceNumbers(1), 1 ), faceNumber, 1 ) ) &
            ) then 
                subCellIndexesRow = 1
            else
                subCellIndexesRow = 2
            end if 
           

            if ( .not. centerCellDataBuffer%isParentCell ) then
                ! Fill Connection
                call pr_FillCellBufferFromSubCell( this, parentCellDataBuffer, &
                                       subCellIndexes( subCellIndexesRow, 1 ), & 
                                       subCellIndexes( subCellIndexesRow, 2 ), & 
                                   neighborCellsBuffer(1), forceCellRefinement )
                neighborCellsBuffer(1)%requestedFromDirection = directionId
            else
                ! If center cell is also parent, then 
                ! it detected a grand parent, so it should
                ! fill the final buffer with a double refinement

                ! Fill aux data buffer, with the same sub cell indexes
                call pr_FillCellBufferFromSubCell( this, parentCellDataBuffer, &
                                       subCellIndexes( subCellIndexesRow, 1 ), & 
                                       subCellIndexes( subCellIndexesRow, 2 ), & 
                                                     auxCellDataBuffer, .true. )
                auxCellDataBuffer%isParentCell = .true.
                auxCellDataBuffer%CellNumber   = parentCellDataBuffer%CellNumber

                ! Fill definitive buffer
                call pr_FillCellBufferFromSubCell( this, auxCellDataBuffer, &
                                       subCellIndexes( subCellIndexesRow, 1 ), & 
                                       subCellIndexes( subCellIndexesRow, 2 ), & 
                                   neighborCellsBuffer(1), forceCellRefinement )
                neighborCellsBuffer(1)%requestedFromDirection = directionId


            end if

        else

            if ( .not. centerCellDataBuffer%isParentCell ) then

                ! If equal size cell, fill buffer and leave
                call pr_FillCellFromRealData( this,                          &
                    centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ), &
                                 neighborCellsBuffer(1), forceCellRefinement )
                neighborCellsBuffer(1)%requestedFromDirection = directionId

            else

                ! Which face connection ?
                if ( directionId .eq. 2 ) then 
                    subColumn = centerCellDataBuffer%parentSubColumn
                    select case ( faceNumber )
                        case(3)
                            subRow = 1
                        case(4)
                            subRow = 2
                    end select
                else if ( directionId .eq. 1 ) then
                    subRow = centerCellDataBuffer%parentSubRow
                    select case ( faceNumber )
                        case(1)
                            subColumn = 2
                        case(2)
                            subColumn = 1
                    end select
                end if

                call pr_FillCellFromRealData( this,                          &
                    centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ), &
                                                parentCellDataBuffer, .true. )

                call pr_FillCellBufferFromSubCell( this, parentCellDataBuffer, &
                                                           subRow,  subColumn, & 
                                   neighborCellsBuffer(1), forceCellRefinement )
                neighborCellsBuffer(1)%requestedFromDirection = directionId

            end if
        end if 


        ! Done
        return


    else
        ! No connection

        ! Done
        return

    end if 


end subroutine pr_FillNeighborCellsConnectionFromHorizontalFace


! RWPT
subroutine pr_FillCellBufferFromSubCell( this, parentCellBuffer, subRow, subColumn, &
                                                    cellBuffer, forceCellRefinement )
    !----------------------------------------------------------------------------- 
    !----------------------------------------------------------------------------- 
    class(ParticleTrackingEngineType), target :: this
    ! input
    type(ModpathCellDataType), intent(in)  :: parentCellBuffer
    integer, intent(in) :: subRow, subColumn 
    logical, intent(in) :: forceCellRefinement
    ! output
    type(ModpathCellDataType), intent(out)  :: cellBuffer
    ! local
    doubleprecision, dimension(6) :: faceFlows
    !-----------------------------------------------------------------------------

    ! Reset the output buffer
    call cellBuffer%Reset()

    ! Fill faceFlows
    call parentCellBuffer%FillSubCellFaceFlowsBuffer( subRow, subColumn, faceFlows )

    ! Given faceFlows, fill cellBuffer faceFlows
    cellBuffer%Q1(1) = faceFlows(1)
    cellBuffer%Q2(1) = faceFlows(2)
    cellBuffer%Q3(1) = faceFlows(3)
    cellBuffer%Q4(1) = faceFlows(4)
    cellBuffer%Q5(1) = faceFlows(5)
    cellBuffer%Q6(1) = faceFlows(6)

    if ( forceCellRefinement ) then 
        ! Fill sources and sinks with the same method 
        ! employed when computing sub cell flows, using 
        ! parent source/sink/storage flows
        cellBuffer%SourceFlow  = parentCellBuffer%SourceFlow  / 4d0
        cellBuffer%SinkFlow    = parentCellBuffer%SinkFlow    / 4d0
        cellBuffer%StorageFlow = parentCellBuffer%StorageFlow / 4d0

        cellBuffer%SubCellRowCount    = 2
        cellBuffer%SubCellColumnCount = 2
        call cellBuffer%ComputeSubCellFlows()
    end if 

    ! "subCell" properties
    cellBuffer%fromSubCell = .true.
    if ( parentCellBuffer%fromSubCell .and. parentCellBuffer%isParentCell ) then 
        cellBuffer%fromSubSubCell = .true.
    end if 
    cellBuffer%parentCellNumber = parentCellBuffer%CellNumber
    cellBuffer%parentSubRow     = subRow
    cellBuffer%parentSubColumn  = subColumn 

    ! Fill dimensions:
    ! Half horizontal, same vertical
    cellBuffer%DX     = 0.5*parentCellBuffer%DX 
    cellBuffer%DY     = 0.5*parentCellBuffer%DY 
    cellBuffer%Bottom = parentCellBuffer%Bottom
    cellBuffer%Top    = parentCellBuffer%Top

    ! Transport properties
    cellBuffer%Porosity    = parentCellBuffer%Porosity
    cellBuffer%Retardation = parentCellBuffer%Retardation


    ! Done
    return 


end subroutine pr_FillCellBufferFromSubCell


! RWPT
subroutine pr_FillCellFromRealData( this, cellNumber, cellDataBuffer, forceCellRefinement )
    !----------------------------------------------------------
    !----------------------------------------------------------
    class(ParticleTrackingEngineType) :: this
    ! input
    integer, intent(in) :: cellNumber
    logical, intent(in) :: forceCellRefinement
    ! output
    type(ModpathCellDataType), intent(inout)  :: cellDataBuffer
    !----------------------------------------------------------

    ! Reset the buffer
    call cellDataBuffer%Reset()

    ! If no valid cellNumber, leave
    if ( cellNumber .le. 0 ) then
        ! Done
        return
    end if  

    ! Fill buffer from real (grid) cellNumber
    call this%FillCellBuffer( cellNumber, cellDataBuffer )

    if ( forceCellRefinement ) then
        ! Forces refinement and sub cell flows computation
        cellDataBuffer%SubCellRowCount    = 2
        cellDataBuffer%SubCellColumnCount = 2
        call cellDataBuffer%ComputeSubCellFlows()
    end if


    return


end subroutine pr_FillCellFromRealData



end module ParticleTrackingEngineModule
