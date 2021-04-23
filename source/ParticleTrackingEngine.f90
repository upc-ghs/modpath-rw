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
    type(ModpathCellDataType) :: CellDataBuffer
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

  ! RWPT
  if (this%TrackingOptions%RandomWalkParticleTracking) then
      call this%FillNeighborSubCellData( this%NeighborSubCellData ) 
      ! DEV
      call this%FillNeighborCellData( this%NeighborCellData )
  end if


      ! DEV
      !call this%FillNeighborCellData( this%NeighborCellData )

  continueLoop = .true.
  isTimeSeriesPoint = .false.
  isMaximumTime = .false.
  
  do while(continueLoop)
      ! Check to see if the particle has moved to another cell. If so, load the new cell data
      if(loc%CellNumber .ne. this%TrackCell%CellData%CellNumber) then
          call this%FillCellBuffer(loc%CellNumber, this%TrackCell%CellData)
          ! RWPT
          if (this%TrackingOptions%RandomWalkParticleTracking) then
              call this%FillNeighborSubCellData( this%NeighborSubCellData )
              ! DEV
              call this%FillNeighborCellData( this%NeighborCellData )
          end if
              !! DEV
              !call this%FillNeighborCellData( this%NeighborCellData )
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
          call this%TrackCell%ExecuteRandomWalkParticleTracking(loc, stopTime, this%TrackCellResult, &
                                                       this%NeighborSubCellData, this%NeighborCellData)
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
         call WriteTraceData(traceModeUnit, this%TrackCell,                     &
           this%TrackCellResult, this%GetCurrentStressPeriod(),                 &
           this%GetCurrentTimeStep())
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
    ! TODO 
    !------------------------------------------------------------------
    !
    ! Specifications
    !------------------------------------------------------------------
    implicit none
    class(ParticleTrackingEngineType),target :: this
    type(ModpathCellDataType), dimension(18) :: neighborCellData
    type(ModpathCellDataType), dimension(8)  :: neighborCellBuffer
    type(ModpathCellDataType), target        :: auxCellDataBuffer
    type(ModpathCellDataType), target        :: parentCellDataBuffer
    integer :: n, m, id, cellCounter, firstNeighborFaceNumber
    integer, dimension(2)   :: orthogonalFaceNumbers, connectedSubCellIndexes
    integer, dimension(2,2) :: subCellIndexes
    doubleprecision, dimension(6) :: faceFlows
    !------------------------------------------------------------------


    ! Verify what is going on with the current TrackCell
    ! Is it refined ?
    print *, '*******************************************************************'
    if ( this%TrackCell%CellData%GetSubCellCount() .gt. 1 ) then 
        print *, '** ParticleTrackingEngine.FillNeighborCellData: TrackCell is refined ', &
        this%TrackCell%CellData%CellNumber, this%TrackCell%CellData%GetSubCellCount(), ' subcells.'
    else
        print *, '** ParticleTrackingEngine.FillNeighborCellData: TrackCell NOT REFINED ', &
        this%TrackCell%CellData%CellNumber, this%TrackCell%CellData%GetSubCellCount(), ' subcells.'

        ! Force refinement
        ! Force two rows and two columns 
        ! to compute internal flows 
        this%TrackCell%CellData%SubCellRowCount    = 2
        this%TrackCell%CellData%SubCellColumnCount = 2
        call this%TrackCell%CellData%ComputeSubCellFlows()
        ! Is this necessary ?

    end if 

    ! Loop through cell faces and fill the neighborhood
    cellCounter = 0
    do n = 1, 6

        ! Each cell face gives three neighbor cells. 
        ! Itself, and two other orthogonal connections
        cellCounter = ( n - 1 )*3 + 1

        ! From the cell connected at face n,
        ! create TrackCells from two connections 
        ! starting from neighborFaceNumber  
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

        ! Now, how each buffer is filled
        ! is handled internally by the
        ! function FillNeighborCellSubBuffer
        call pr_FillNeighborCellsSubBuffer( &
            this, this%TrackCell%CellData, n, firstNeighborFaceNumber, neighborCellData, cellCounter ) 


    end do ! End loop through cell faces


    ! Print a report
    print *, '********************************************************************************************'
    print *, '** FillNeighborCellData: neighbors information for TrackCell', this%TrackCell%CellData%CellNumber

    do n = 1, 18
        print *, '** FillNeighborCellData: neighbor cell ', n, ' cellNumber ', neighborCellData(n)%CellNumber
        if ( neighborCellData(n)%CellNumber .eq. -999 ) then
            print *, '** FillNeighborCellData: custom cell has sub cells '
            do m = 1,4
                print *, neighborCellData(n)%SubCellDataBuffer( m )%CellNumber
            end do 
        end if
    end do



end subroutine pr_FillNeighborCellData


! RWPT
subroutine pr_FillCellBufferFromSubCell( this, parentCellBuffer, subRow, subColumn, cellBuffer )
    ! 
    ! asd TODO
    ! 
    class(ParticleTrackingEngineType), target :: this
    ! input
    type(ModpathCellDataType), intent(in)  :: parentCellBuffer
    integer :: subRow, subColumn 
    ! output
    type(ModpathCellDataType), intent(out)  :: cellBuffer
    ! local
    doubleprecision, dimension(6) :: faceFlows
    !-----------------------------------------------------------------

    ! Fill faceFlows
    call parentCellBuffer%FillSubCellFaceFlowsBuffer( subRow, subColumn, faceFlows )

    ! Given faceFlows, fill cellBuffer faceFlows
    cellBuffer%Q1(1)       = faceFlows(1)
    cellBuffer%Q2(1)       = faceFlows(2)
    cellBuffer%Q3(1)       = faceFlows(3)
    cellBuffer%Q4(1)       = faceFlows(4)
    cellBuffer%Q5(1)       = faceFlows(5)
    cellBuffer%Q6(1)       = faceFlows(6)

    ! Fill sources/sinks/storage
    cellBuffer%SourceFlow  = parentCellBuffer%SourceFlow / 4d0
    cellBuffer%SinkFlow    = parentCellBuffer%SinkFlow / 4d0
    cellBuffer%StorageFlow = parentCellBuffer%StorageFlow / 4d0

    ! Force refinement 
    cellBuffer%SubCellRowCount    = 2
    cellBuffer%SubCellColumnCount = 2
    call cellBuffer%ComputeSubCellFlows()

    return 

end subroutine pr_FillCellBufferFromSubCell



recursive function pr_FillCellFromRefinedCells( this, currentCellData, faceNumber ) result(filledCellData)
    !
    ! TODO
    !
    ! 


    class(ParticleTrackingEngineType), target :: this
    ! input
    type(ModpathCellDataType), intent(in)  :: currentCellData
    integer, intent(in)                    :: faceNumber
    integer, dimension(2)                  :: subCellIndexes, newSubCellIndexes, connSubCellIndexes
    ! output
    type(ModpathCellDataType) :: filledCellData
    ! local 
    integer :: m, id , nid, cnid
    type(ModpathCellDataType) :: auxCellData
    type(ModpathCellDataType), dimension(:), pointer :: cellDataBuffer


    ! Assumptions
    !   - Smoothed unstructured grid
    !   - Horizontal faceNumber (?)


    ! Initialize filledCellData
    call filledCellData%Reset()

    ! Allocate celldatabuffer and assign pointer
    allocate( filledCellData%SubCellDataBuffer(4) )
    cellDataBuffer => filledCellData%SubCellDataBuffer
    ! Deallocation at some point ?


    ! Resolution if the cell is composed of
    ! a subcelldatabuffer, id -999
    if ( currentCellData%CellNumber .eq. -999 ) then

        print *, '**FillCellFromRefinedCells: initializing from a custom cell -999, faceNumber', faceNumber 

        do m = 1,4
            if( currentCellData%SubCellDataBuffer(m)%CellNumber .le. 0 ) then
                print *, '**FillCellFromRefinedCells: cell ', m, ' from buffer is not initialized, BEFORE'
            else
                print *, '**FillCellFromRefinedCells: cell ', m,&
                    ' from buffer is initialized ENTER', currentCellData%SubCellDataBuffer(m)%CellNumber
            end if
        end do


        ! subCellIndexes: location of subtrackcells in local buffer
        ! newSubCellIndexes: location of subtrackcells in buffer of cell to be created 
        ! connSubCellIndexes: location in new buffer of subtrackcell connected to direct connection
        select case (faceNumber)
            case (1)
                subCellIndexes(1)     = 1 
                subCellIndexes(2)     = 3
                newSubCellIndexes(1)  = 2 
                newSubCellIndexes(2)  = 4
                connSubCellIndexes(1) = 1 
                connSubCellIndexes(2) = 3
            case (2)
                subCellIndexes(1)     = 2 
                subCellIndexes(2)     = 4
                newSubCellIndexes(1)  = 1 
                newSubCellIndexes(2)  = 3
                connSubCellIndexes(1) = 2 
                connSubCellIndexes(2) = 4
            case (3)
                subCellIndexes(1)     = 3 
                subCellIndexes(2)     = 4
                newSubCellIndexes(1)  = 1 
                newSubCellIndexes(2)  = 2
                connSubCellIndexes(1) = 3 
                connSubCellIndexes(2) = 4
            case (4)
                subCellIndexes(1)     = 1 
                subCellIndexes(2)     = 2
                newSubCellIndexes(1)  = 3 
                newSubCellIndexes(2)  = 4
                connSubCellIndexes(1) = 1 
                connSubCellIndexes(2) = 2
        end select 


        ! It is possible for the custom cell
        ! to be connected to a bigger cell
        ! Verify if connections through faceNumber
        ! are equal, if cells are not custom
        if ( &
            ( currentCellData%SubCellDataBuffer( subCellIndexes(1) )%CellNumber .gt. 0 ) .and. &
            ( currentCellData%SubCellDataBuffer( subCellIndexes(2) )%CellNumber .gt. 0 )       &
        ) then
            ! In the context of a smoothed grid, if one of 
            ! the sub cells is custom, it means that 
            ! was built from smaller cells. In such case,
            ! the verification of being connected to a bigger cell
            ! makes no sense because it would
            ! imply the existence of a double jump in refinement
            ! between contiguous cells  

            ! If both subcells connected to the same cell,
            ! then it is a bigger one
            if ( &
                ( currentCellData%SubCellDataBuffer( subCellIndexes(1) )%GetFaceConnection( faceNumber, 1 ) ) .eq. &
                ( currentCellData%SubCellDataBuffer( subCellIndexes(2) )%GetFaceConnection( faceNumber, 1 ) )      &
            ) then
                print *, '**FillCellFromRefinedCells: sub cells related to faceNumber have the same connection, BIGGER'
                print *, '**FillCellFromRefinedCells: the connected cell number is cell: ', &
                    currentCellData%SubCellDataBuffer( subCellIndexes(2) )%GetFaceConnection( faceNumber, 1 )
                ! In this case, return filledCellData, initialized from
                ! real cell data, and refined.
                ! In this case, although this cell passed through
                ! the function of refined cells, is not created with 
                ! the method of filling the subcelldatabuffer and a custom
                ! cellNumber. 

                ! In any case this is a terminal case in building neighbors

                ! It is possible that this section resolution 
                ! is better placed in another function/place

                ! Fill the output buffer with the information of the bigger
                ! cell, refine it by computing subcell flows and leave
                call pr_FillCellFromRealCellData( this, &
                    currentCellData%SubCellDataBuffer( subCellIndexes(1) )%GetFaceConnection( faceNumber, 1 ), filledCellData )

                ! Done
                return

            end if

        end if


        ! Following section verify connections
        ! through subfaces and build buffer accordingly 

        ! Loop over horizontal sub faces. 
        ! Each loop should fill two entries 
        ! of cellDataBuffer
        do m = 1, 2

            ! Indexes
            id   = subCellIndexes( m )
            nid  = newSubCellIndexes( m ) 
            cnid = connSubCellIndexes( m )

            ! SubCell in buffer could also be a -999.
            if ( currentCellData%SubCellDataBuffer( id )%CellNumber .eq. -999 ) then
                print *, '**FillCellFromRefinedCells: one of the subcells is custom -999'

                ! Compare connections through faceNumber
                if ( &
                    currentCellData%SubCellDataBuffer( id )%SubCellDataBuffer(      &
                        subCellIndexes(1) )%GetFaceConnection( faceNumber, 1 ) .eq. &
                    currentCellData%SubCellDataBuffer( id )%SubCellDataBuffer(      &
                        subCellIndexes(2) )%GetFaceConnection( faceNumber, 1 ) ) then

                    ! Extract the cell and fill
                    ! direct connection
                    call this%FillCellBuffer( &
                        currentCellData%SubCellDataBuffer( id )%SubCellDataBuffer(  &
                            subCellIndexes(1) )%GetFaceConnection( faceNumber, 1 ), &
                                                              cellDataBuffer( nid ) )

                    ! Fill secondary connection
                    if ( cellDataBuffer( nid )%SubFaceCounts( faceNumber ) .gt. 1 ) then
                        print *, '*FillCellFromRefinedCells: going to recursion' 
                        ! recursion with direct assignment
                        cellDataBuffer( cnid ) = pr_FillCellFromRefinedCells( this, cellDataBuffer( nid ), faceNumber )
                    else if ( cellDataBuffer( nid )%SubFaceCounts( faceNumber ) .eq. 1 ) then 
                        print *, '*FillCellFromRefinedCells: just one connection' 
                        ! fill buffer and move on
                        ! Remember that filling of cellDataBuffer is not generalized yet.
                        ! It has to be independent from faceNumber orientation.
                        call this%FillCellBuffer( cellDataBuffer( nid )%GetFaceConnection( faceNumber, 1 ), cellDataBuffer( cnid ) )
                    else
                        ! Is this possible ?
                        print *, '*FillCellFromRefinedCells: something else, n subfaces', &
                            cellDataBuffer( id )%SubFaceCounts( faceNumber )
                    end if 

                else

                    print *, '**FillCellFromRefinedCells: going to recursion'

                    ! Recursion with direct assignment
                    cellDataBuffer( nid ) = pr_FillCellFromRefinedCells( this, currentCellData%SubCellDataBuffer( id ), faceNumber )

                    ! Just a verification in the meantime
                    if ( cellDataBuffer(nid)%CellNumber .gt. 0 ) then
                        print *, '**FillCellFromRefinedCells: the BUFFER is filled with REAL CELL DATA'
                    end if

                    ! Verify connection through faceNumber
                    if ( & 
                      cellDataBuffer( nid )%SubCellDataBuffer( subCellIndexes(1) )%GetFaceConnection( faceNumber, 1 ) .eq. &
                      cellDataBuffer( nid )%SubCellDataBuffer( subCellIndexes(2) )%GetFaceConnection( faceNumber, 1 ) ) then

                        ! Extract the cell and fill
                        call this%FillCellBuffer(                    &
                            cellDataBuffer( nid )%SubCellDataBuffer( &
                                 subCellIndexes(1) )%GetFaceConnection( faceNumber, 1 ), &
                                                                  cellDataBuffer( cnid ) )
                    
                    else

                        cellDataBuffer( cnid ) = pr_FillCellFromRefinedCells( this, cellDataBuffer( nid ), faceNumber )

                    end if

                end if
                   

            else if ( currentCellData%SubCellDataBuffer( id )%CellNumber .gt. 0 ) then

                print *, '**FillCellFromRefinedCells: sub cell is from real data'

                ! Initialize direct connection
                ! through faceNumber
                if ( currentCellData%SubCellDataBuffer( id )%SubFaceCounts( faceNumber ) .gt. 1 ) then
                    print *, '*FillCellFromRefinedCells: going to recursion' 
                    ! Recursion with direct assignment
                    cellDataBuffer( nid ) = pr_FillCellFromRefinedCells( this, currentCellData%SubCellDataBuffer( id ), faceNumber )
                else if ( currentCellData%SubCellDataBuffer( id )%SubFaceCounts( faceNumber ) .eq. 1 ) then 
                    print *, '*FillCellFromRefinedCells: just one connection'
                    print *, '*FillCellFromRefinedCellds: cell:', currentCellData%SubCellDataBuffer( id )%CellNumber, &
                           ' connected to ', currentCellData%SubCellDataBuffer( id )%GetFaceConnection( faceNumber, 1 ) 
                    ! Fill buffer and move on
                    call this%FillCellBuffer( currentCellData%SubCellDataBuffer( id )%GetFaceConnection( faceNumber, 1 ), &
                                                                                                    cellDataBuffer( nid ) )
                else
                    ! Is this possible ?
                    print *, '*FillCellFromRefinedCells: something else, n subfaces BLEHEHE'
                end if 

                ! Initialize secondary connection,
                ! cell connected to direct connection through
                ! faceNumber
                if ( cellDataBuffer( nid )%SubFaceCounts( faceNumber ) .gt. 1 ) then
                    print *, '*FillCellFromRefinedCells: going to recursion' 
                    ! Recursion with direct assignment
                    cellDataBuffer( cnid ) = pr_FillCellFromRefinedCells( this, cellDataBuffer( nid ), faceNumber )
                else if ( cellDataBuffer( nid )%SubFaceCounts( faceNumber ) .eq. 1 ) then 
                    print *, '*FillCellFromRefinedCells: just one connection' 
                    ! Fill buffer and move on
                    call this%FillCellBuffer( cellDataBuffer( nid )%GetFaceConnection( faceNumber, 1 ), cellDataBuffer( cnid ) )
                else
                    ! Is this possible ?
                    print *, '*FillCellFromRefinedCells: something else, n subfaces', &
                        cellDataBuffer( id )%SubFaceCounts( faceNumber )
                end if 

            else

                print *, '**FillCellFromRefinedCells: sub cell is empty, WHY ?'

            end if
           
            
        end do


    else

        ! If the center cell is a real one,
        ! request information of connected cells.
        ! This situation is the one ocurring in the very first 
        ! request process of neighbor cells, these are requested
        ! from a real center cell. Thus it is expected that
        ! that connections comply with smoothed grid.
        print *, '**FillCellFromRefinedCells: in the smooth unstructured case' 


        ! Loop over face connections and extract 
        ! their connected cell through same faceNumber
        do m = 1, currentCellData%SubFaceCounts( faceNumber )

            ! Compute cell id in buffer 
            id = (m-1)*currentCellData%SubFaceCounts( faceNumber ) + 1

            ! Fill buffer for connected cell
            ! Implied smoothed unstructured grid
            call this%FillCellBuffer( currentCellData%GetFaceConnection( faceNumber, m ), cellDataBuffer( id ) )

            ! Ask if the filled buffer is connected to multiple
            ! or single cell
            if ( cellDataBuffer( id )%SubFaceCounts( faceNumber ) .gt. 1 ) then
                print *, '**FillCellFromRefinedCells: going to recursion' 
                ! Recursion with direct assignment
                cellDataBuffer( id + 1 ) = pr_FillCellFromRefinedCells( this, cellDataBuffer( id ), faceNumber )
            else if ( cellDataBuffer( id )%SubFaceCounts( faceNumber ) .eq. 1 ) then 
                print *, '**FillCellFromRefinedCells: just one connection' 
                ! Fill buffer and move on
                call this%FillCellBuffer( cellDataBuffer( id )%GetFaceConnection( faceNumber, 1 ), cellDataBuffer( id + 1 ) )
            else
                ! Is this possible ?
                print *, '*FillCellFromRefinedCells: something else, n subfaces', cellDataBuffer( id )%SubFaceCounts( faceNumber )
            end if 

        end do


    end if    


    ! Something to recognize its a custom/linked cell
    filledCellData%CellNumber = -999

    ! Fill face flows based on buffer information
    ! Intercell 
    filledCellData%SubCellFlows(1) = cellDataBuffer(1)%GetAveragedFaceFlow(2)
    filledCellData%SubCellFlows(2) = cellDataBuffer(3)%GetAveragedFaceFlow(2)
    filledCellData%SubCellFlows(3) = cellDataBuffer(3)%GetAveragedFaceFlow(4)
    filledCellData%SubCellFlows(4) = cellDataBuffer(4)%GetAveragedFaceFlow(4) 

    ! Verify indexation
    ! Horizontal faces 
    filledCellData%Q1(1) = cellDataBuffer(1)%GetAveragedFaceFlow(1)
    filledCellData%Q1(2) = cellDataBuffer(3)%GetAveragedFaceFlow(1)
    filledCellData%Q2(1) = cellDataBuffer(2)%GetAveragedFaceFlow(2)
    filledCellData%Q2(2) = cellDataBuffer(4)%GetAveragedFaceFlow(2)
    filledCellData%Q3(1) = cellDataBuffer(3)%GetAveragedFaceFlow(3)
    filledCellData%Q3(2) = cellDataBuffer(4)%GetAveragedFaceFlow(3)
    filledCellData%Q4(1) = cellDataBuffer(1)%GetAveragedFaceFlow(4)
    filledCellData%Q4(2) = cellDataBuffer(2)%GetAveragedFaceFlow(4)

    ! Vertical faces
    filledCellData%Q5(1) = cellDataBuffer(1)%GetAveragedFaceFlow(5)
    filledCellData%Q5(2) = cellDataBuffer(2)%GetAveragedFaceFlow(5)
    filledCellData%Q5(3) = cellDataBuffer(3)%GetAveragedFaceFlow(5)
    filledCellData%Q5(4) = cellDataBuffer(4)%GetAveragedFaceFlow(5)
    filledCellData%Q6(1) = cellDataBuffer(1)%GetAveragedFaceFlow(6)
    filledCellData%Q6(2) = cellDataBuffer(2)%GetAveragedFaceFlow(6)
    filledCellData%Q6(3) = cellDataBuffer(3)%GetAveragedFaceFlow(6)
    filledCellData%Q6(4) = cellDataBuffer(4)%GetAveragedFaceFlow(6)

    ! Fill sources and sinks ?

    return 


end function pr_FillCellFromRefinedCells




subroutine pr_FillCellFromRealCellData( this, cellNumber, cellDataBuffer )
    !----------------------------------------------------------
    ! Forces computation of internal flows by redefining 
    ! the number of subcells.
    ! TODO
    !----------------------------------------------------------
    class(ParticleTrackingEngineType) :: this
    ! input
    integer, intent(in)                    :: cellNumber
    ! output
    type(ModpathCellDataType), intent(inout)  :: cellDataBuffer
    !----------------------------------------------------------

    ! If no connection, reset buffer and leave
    if ( cellNumber .le. 0 ) then
        print *, 'FillCellFromRealCellData: reset cell to zeroes' 
        call cellDataBuffer%Reset()
        return
    end if  

    ! Fill buffer from real (grid) cellNumber
    call this%FillCellBuffer( cellNumber, cellDataBuffer )

    ! Forces refinement and sub cell flows computation
    cellDataBuffer%SubCellRowCount    = 2
    cellDataBuffer%SubCellColumnCount = 2
    call cellDataBuffer%ComputeSubCellFlows()

    return

end subroutine pr_FillCellFromRealCellData




!RWPT
subroutine pr_FillNeighborCellsSubBuffer(       &
        this, centerCellDataBuffer, faceNumber, &
        firstNeighborFaceNumber, neighborCellsDataBuffer, directConnectionBufferIndex )
    !----------------------------------------------------------------------------------
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
    type(ModpathCellDataType), intent(in)                    :: centerCellDataBuffer
    type(ModpathCellDataType), dimension(18), intent(inout)  :: neighborCellsDataBuffer
    integer, intent(in) :: faceNumber, firstNeighborFaceNumber, directConnectionBufferIndex

    type(ModpathCellDataType), target        :: auxCellDataBuffer
    type(ModpathCellDataType), target        :: parentCellDataBuffer
    integer, dimension(2)         :: orthogonalFaceNumbers, connectedSubCellIndexes
    integer, dimension(2,2)       :: subCellIndexes
    doubleprecision, dimension(6) :: faceFlows

    integer :: directConnectionSubRow, indirectConnectionSubRow
    integer :: indirectConnectionBufferIndex
    integer :: lastConnectionBufferIndex    
    integer :: lastConnectionFaceNumber     
    integer :: n, m 
    !---------------------------------------------------------------------



    ! Detect orientation of faceNumber
    if ( ( faceNumber .eq. 5 ) .or. ( faceNumber .eq. 6 ) ) then 
        print *, '**FillNeighborCellsSubBuffer: will initialize from a vertical faceNumber', faceNumber
        ! If z-axis faceNumber
        ! Initialize direct connection buffer with connected cellNumber, 
        ! supposing no vertical refinement (single vertical connection)

        print *, '**FillNeighborCellsSubBuffer: connected to cell ', centerCellDataBuffer%GetFaceConnection( faceNumber, 1 )

        call pr_FillCellFromRealCellData( this, centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ), &
                                                  neighborCellsDataBuffer( directConnectionBufferIndex ) )

        ! There are two possible approaches here:
        !   - Assume/know remaining neighbors are horizontal and fill ! THIS 
        !   - Ask if remaining neighbors are horizontal and proceed 

        ! Process neighbor cells, using as center the already 
        ! filled directConnectionBuffer 
        !   - neighborCellsDataBuffer( directConnectionBufferIndex + 1 ) solving firstNeighborFaceNumber
        !   - neighborCellsDataBuffer( directConnectionBufferIndex + 2 ) solving firstNeighborFaceNumber + 1

        ! Both of previous buffers are solved with method of horizontal
        ! faceNumbers following below (else branch)

        ! Determine and fill only direct horizontal connections is a simplified 
        ! process of what is implemented downstream 

        ! Thus, it will be done manually now and then 
        ! it might be comparmentilized into a proper function 

        ! Take already filled buffer as the center buffer

        ! Remember that 
        ! firstNeighborFaceNumber linked to directConnectionBufferIndex + 1 
        ! firstNeighborFaceNumber + 1 linked to directConnectionBufferIndex + 2
        call pr_FillDirectBufferFromHorizontalFace( this,                                     & 
            neighborCellsDataBuffer( directConnectionBufferIndex ), firstNeighborFaceNumber,  &
                                   neighborCellsDataBuffer( directConnectionBufferIndex + 1 ) )

        call pr_FillDirectBufferFromHorizontalFace( this,                                         & 
            neighborCellsDataBuffer( directConnectionBufferIndex ), firstNeighborFaceNumber + 1,  &
                                       neighborCellsDataBuffer( directConnectionBufferIndex + 2 ) )

    else
        print *, '**FillNeighborCellsSubBuffer: will initialize from a horizontal faceNumber', faceNumber
        ! If horizontal faceNumber, 
        ! go through the process of connection 
        ! state detection.

        ! Connection state detection
        ! Connected to one or more cells ?
        if ( centerCellDataBuffer%SubFaceCounts( faceNumber ) .gt. 1 ) then
            print *, '**FillNeighborCellsSubBuffer: found more than one connection'

            ! If more than one connection
            ! Fill from refined cells
            neighborCellsDataBuffer( directConnectionBufferIndex ) = &
                pr_FillCellFromRefinedCells( this, centerCellDataBuffer, faceNumber )

            print *,  '** FillNeighborSubBuffer: custom buffer filled with indexes'
            do n = 1,4
                print *, neighborCellsDataBuffer( directConnectionBufferIndex )%SubCellDataBuffer(n)%CellNumber
            end do
            

            if ( firstNeighborFaceNumber .ge. 5 ) then  
                print *, '**FillNeighborCellsSubBuffer: will request vertical faces'
                ! If the others requested neighbor cells are through 
                ! vertical faces, then extract cellNumber of cells
                ! connected through vertical faces and FillFromRefinedCells
                ! taking center cell as the vertical connection and 
                ! fill requesting from faceNumber

                call pr_FillCellFromRealCellData( this,                                   &
                    centerCellDataBuffer%GetFaceConnection( firstNeighborFaceNumber, 1 ), &
                                                                        auxCellDataBuffer )
                if ( auxCellDataBuffer%CellNumber .le. 0 ) then
                    ! If no vertical connection, reset buffer
                    call neighborCellsDataBuffer( directConnectionBufferIndex + 1 )%Reset()
                else
                    neighborCellsDataBuffer( directConnectionBufferIndex + 1 ) = &
                        pr_FillCellFromRefinedCells( this, auxCellDataBuffer, faceNumber )
                end if 
                                                        
                call pr_FillCellFromRealCellData( this,                                       &
                    centerCellDataBuffer%GetFaceConnection( firstNeighborFaceNumber + 1, 1 ), &
                                                                            auxCellDataBuffer )
                if ( auxCellDataBuffer%CellNumber .le. 0 ) then
                    ! If no vertical connection, reset buffer
                    call neighborCellsDataBuffer( directConnectionBufferIndex + 2 )%Reset()
                else
                    neighborCellsDataBuffer( directConnectionBufferIndex + 2 ) = &
                        pr_FillCellFromRefinedCells( this, auxCellDataBuffer, faceNumber )
                end if 

                ! Done
                return

            else 
                print *, '**FillNeighborCellsSubBuffer: will request horizontal faces'
                ! If the others requested neighbor cells are through 
                ! horizontal faces, change centerCellDataBuffer to the 
                ! recently filled buffer and request info 
                neighborCellsDataBuffer( directConnectionBufferIndex + 1 ) = &
                    pr_FillCellFromRefinedCells( this, neighborCellsDataBuffer( directConnectionBufferIndex ), &
                                                                                       firstNeighborFaceNumber )

            print *,  '** FillNeighborSubBuffer: FIRST NEIGHBOR custom buffer filled with indexes'
            do n = 1,4
                print *, neighborCellsDataBuffer( directConnectionBufferIndex + 1 )%SubCellDataBuffer(n)%CellNumber
            end do

                neighborCellsDataBuffer( directConnectionBufferIndex + 2 ) = &
                    pr_FillCellFromRefinedCells( this, neighborCellsDataBuffer( directConnectionBufferIndex ), &
                                                                                 firstNeighborFaceNumber  +  1 )

            print *,  '** FillNeighborSubBuffer: SECOND NEIGHBOR custom buffer filled with indexes'
            do n = 1,4
                print *, neighborCellsDataBuffer( directConnectionBufferIndex + 2 )%SubCellDataBuffer(n)%CellNumber
            end do

            print *,  '** FillNeighborSubBuffer: VERIFYING DIRECT custom buffer filled with indexes'
            do n = 1,4
                print *, neighborCellsDataBuffer( directConnectionBufferIndex  )%SubCellDataBuffer(n)%CellNumber
            end do
                ! Done 
                return

            end if

        else if ( centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ) .gt. 0 ) then
            print *, '**FillNeighborCellsSubBuffer: there is at least one connection', &
                centerCellDataBuffer%GetFaceConnection(faceNumber, 1 )
            ! If one connection, verify if is an equal size 
            ! or bigger cell
            
            ! Information required for horizontal neighbor cell size 
            ! detection. This should be defined 
            ! in a different place maybe. 
            select case ( faceNumber )
                case (1)
                    orthogonalFaceNumbers(1)   = 3 
                    orthogonalFaceNumbers(2)   = 4
                    connectedSubCellIndexes(1) = 2
                    connectedSubCellIndexes(2) = 4
                    subCellIndexes(1,:)        = (1,2) 
                    subCellIndexes(2,:)        = (2,2) 
                case (2)
                    orthogonalFaceNumbers(1)   = 3 
                    orthogonalFaceNumbers(2)   = 4
                    connectedSubCellIndexes(1) = 1
                    connectedSubCellIndexes(2) = 3
                    subCellIndexes(1,:)        = (1,1) 
                    subCellIndexes(2,:)        = (2,1) 
                case (3)
                    orthogonalFaceNumbers(1)   = 1 
                    orthogonalFaceNumbers(2)   = 2
                    connectedSubCellIndexes(1) = 2
                    connectedSubCellIndexes(2) = 1
                    subCellIndexes(1,:)        = (1,2) 
                    subCellIndexes(2,:)        = (1,1) 
                case (4)
                    orthogonalFaceNumbers(1)   = 1 
                    orthogonalFaceNumbers(2)   = 2
                    connectedSubCellIndexes(1) = 4
                    connectedSubCellIndexes(2) = 3
                    subCellIndexes(1,:)        = (2,2) 
                    subCellIndexes(2,:)        = (2,1) 
            end select 

            
            ! Neighbor connection size detection 

            ! Analyze if current TrackCell is 
            ! connected to a bigger TrackCell. 
            ! Condition verifies if connections trough faceNumber
            ! of cells connected through orthogonal faces to current TrackCell, 
            ! are the same or different. If the same, then current TrackCell is 
            ! connected to a bigger cell.

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
                print *, '**FillNeighborCellsSubBuffer: connected to a bigger cell '
               
                ! Fill the parentCellDataBuffer
                call pr_FillCellFromRealCellData( this, centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ), &
                                                                                            parentCellDataBuffer )

                ! Location of current TrackCell relative to bigger cell
                ! defines indexes employed for filling buffers. 
                if ( &
                    ( centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ) .eq.                               &
                      this%Grid%GetFaceConnection(                                                               & 
                        centerCellDataBuffer%GetFaceConnection( orthogonalFaceNumbers(1), 1 ), faceNumber, 1 ) ) &
                ) then 
                    directConnectionSubRow        = 1
                    indirectConnectionSubRow      = 2
                    indirectConnectionBufferIndex = directConnectionBufferIndex + 1
                    lastConnectionBufferIndex     = directConnectionBufferIndex + 2
                    lastConnectionFaceNumber      = firstNeighborFaceNumber + 1 
                else
                    directConnectionSubRow        = 2
                    indirectConnectionSubRow      = 1
                    indirectConnectionBufferIndex = directConnectionBufferIndex + 2
                    lastConnectionBufferIndex     = directConnectionBufferIndex + 1
                    lastConnectionFaceNumber      = firstNeighborFaceNumber 
                end if 

                
                ! Direct Connection
                call pr_FillCellBufferFromSubCell( this,                               &
                    parentCellDataBuffer, subCellIndexes( directConnectionSubRow, 1 ), &
                    subCellIndexes( directConnectionSubRow, 2 ),                       &
                    neighborCellsDataBuffer( directConnectionBufferIndex ) )


                ! Final filling process
                if ( firstNeighborFaceNumber .ge. 5 ) then 
                    ! If remaining neighbor cells are vertical

                    ! Note that the case in which vertical connections
                    ! are requested as neighbor cells, can be solved
                    ! by taking advantage of the no vertical refinement 
                    ! assumption.

                    ! Basically, indexes determined in preliminary
                    ! if branch remain valid for indexation of
                    ! neighbor layers

                    ! So, this stage requires:
                    !   - Get vertical connection through firstNeighborFaceNumber
                    !     and firstNeighborFaceNumber + 1, from parentCellDataBuffer.
                    !   - Fill both buffers as parents buffers, maybe using aux buffer.
                    !   - Initialize cells from indexation solved for previous subcell,
                    !     using as source auxCellDataBuffer which acts as parent. 

                    ! Get cell connected through firstNeighborFaceNumber
                    call pr_FillCellFromRealCellData( this, &
                        parentCellDataBuffer%GetFaceConnection( firstNeighborFaceNumber, 1 ), &
                                                                            auxCellDataBuffer )

                    if ( auxCellDataBuffer%CellNumber .eq. 0 ) then
                        ! Empty buffer
                        call neighborCellsDataBuffer( directConnectionBufferIndex + 1 )%Reset()
                    else
                        ! Fill from subcell
                        call pr_FillCellBufferFromSubCell( this,                            &
                            auxCellDataBuffer, subCellIndexes( directConnectionSubRow, 1 ), &
                            subCellIndexes( directConnectionSubRow, 2 ),                    &
                            neighborCellsDataBuffer( directConnectionBufferIndex + 1 ) )
                    end if

                    ! And same for the last
                    call pr_FillCellFromRealCellData( this, &
                        parentCellDataBuffer%GetFaceConnection( firstNeighborFaceNumber + 1, 1 ), &
                                                                                auxCellDataBuffer )

                    if ( auxCellDataBuffer%CellNumber .eq. 0 ) then
                        ! Empty buffer
                        call neighborCellsDataBuffer( directConnectionBufferIndex + 2 )%Reset()
                    else
                        ! Fill from subcell
                        call pr_FillCellBufferFromSubCell( this,                            &
                            auxCellDataBuffer, subCellIndexes( directConnectionSubRow, 1 ), &
                            subCellIndexes( directConnectionSubRow, 2 ),                    &
                            neighborCellsDataBuffer( directConnectionBufferIndex + 2 ) )
                    end if


                else

                    ! If remaining neighbor cells are requested 
                    ! through horizontal faces

                    ! Indirect connection 
                    call pr_FillCellBufferFromSubCell( this,                                 &
                        parentCellDataBuffer, subCellIndexes( indirectConnectionSubRow, 1 ), &
                        subCellIndexes( indirectConnectionSubRow, 2 ),                       &
                        neighborCellsDataBuffer( indirectConnectionBufferIndex ) )

                    ! Populate last remaining neighbor by analyzing parent connections.
                    if ( parentCellDataBuffer%SubFaceCounts( lastConnectionFaceNumber ) .eq. 1 ) then
                        ! If parent connection through lastConnectionFaceNumber is single,
                        ! then it means big cell connection should undergo
                        ! double refinement.

                        ! Fill aux buffer with trackcell connected through 
                        ! lastConnectionFaceNumber to parentCellDataBuffer
                        call pr_FillCellFromRealCellData( this, &
                            parentCellDataBuffer%GetFaceConnection( lastConnectionFaceNumber, 1 ), &
                                                                                 auxCellDataBuffer )

                        ! Populate cell from subcell
                        call pr_FillCellBufferFromSubCell( this,                               &
                            auxCellDataBuffer, subCellIndexes( indirectConnectionSubRow, 1 ) , &
                            subCellIndexes( indirectConnectionSubRow, 2 ),                     &
                            neighborCellsDataBuffer( lastConnectionBufferIndex ) )


                    else if ( parentCellDataBuffer%SubFaceCounts( lastConnectionFaceNumber ) .gt. 1 ) then  
                        ! If the connection is multiple, then it should fill the 
                        ! buffer with the cell connected through 
                        ! lastConnectionFaceNumber, using column value obtained from 
                        ! subCellIndexes and directConnectionSubRow.

                        call pr_FillCellFromRealCellData( this, &
                           parentCellDataBuffer%GetFaceConnection( lastConnectionFaceNumber, & 
                                              subCellIndexes( directConnectionSubRow, 2 ) ), &
                                        neighborCellsDataBuffer( lastConnectionBufferIndex ) )

                    else
                        ! No connection through lastConnectionFaceNumber, empty buffer
                        call neighborCellsDataBuffer( lastConnectionBufferIndex )%Reset()

                    end if 


                end if 

                ! Done 
                return

            else
                print *, '**FillNeighborCellsSubBuffer: connected to an equal size cell '
                ! If equal size cell, fill buffers and leave,
                ! regardless of orientation/axis of requested faces

                ! Fill direct connection buffer
                call pr_FillCellFromRealCellData( this, centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ), &
                                                          neighborCellsDataBuffer( directConnectionBufferIndex ) )

                ! Fill remaining neighbor cells, taking as center the direct
                ! connection buffer
                call pr_FillCellFromRealCellData( this,                                                                     &
                    neighborCellsDataBuffer( directConnectionBufferIndex )%GetFaceConnection( firstNeighborFaceNumber, 1 ), &
                                                                 neighborCellsDataBuffer( directConnectionBufferIndex + 1 ) )
                call pr_FillCellFromRealCellData( this,                                                                         &
                    neighborCellsDataBuffer( directConnectionBufferIndex )%GetFaceConnection( firstNeighborFaceNumber + 1, 1 ), &
                                                                     neighborCellsDataBuffer( directConnectionBufferIndex + 2 ) )

                ! Done
                return

            end if 

        else

            ! No connection 
            ! Reset buffers
            call neighborCellsDataBuffer( directConnectionBufferIndex     )%Reset()
            call neighborCellsDataBuffer( directConnectionBufferIndex + 1 )%Reset()
            call neighborCellsDataBuffer( directConnectionBufferIndex + 2 )%Reset()

            ! Done
            return

        end if

    end if


end subroutine pr_FillNeighborCellsSubBuffer



subroutine pr_FillDirectBufferFromHorizontalFace( this, centerCellDataBuffer, faceNumber, outputCellBuffer )
    !--------------------------------------------------------------------------------------
    ! TODO
    !--------------------------------------------------------------------------------------
    class(ParticleTrackingEngineType) :: this
    ! input
    type(ModpathCellDataType), intent(in)    :: centerCellDataBuffer
    type(ModpathCellDataType), intent(inout) :: outputCellBuffer
    integer, intent(in) :: faceNumber
    type(ModpathCellDataType), target        :: auxCellDataBuffer
    type(ModpathCellDataType), target        :: parentCellDataBuffer
    integer, dimension(2)                    :: orthogonalFaceNumbers, connectedSubCellIndexes
    integer, dimension(2,2)                  :: subCellIndexes
    integer                                  :: directConnectionSubRow
    !--------------------------------------------------------------------------------------


    ! If invalid center cell, reset buffer and leave 
    if ( centerCellDataBuffer%CellNumber .le. 0 ) then
        print *, '*** FillDirectBufferFromHorizontalFace: center cell is empty, return empty cell'
        call outputCellBuffer%Reset()
        return
    end if 


    ! Verify horizontal connection state
    if ( centerCellDataBuffer%SubFaceCounts( faceNumber ) .gt. 1 ) then 
        ! More than one connection
        ! Fill from refined cells
        outputCellBuffer = pr_FillCellFromRefinedCells( this, centerCellDataBuffer, faceNumber )

    else if ( centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ) .gt. 0 ) then
        ! Only one connection, verify if 
        ! same size or bigger size cell 

        ! MOVE THIS CODE TO A PROPER LOCATION
        select case ( faceNumber )
            case (1)
                orthogonalFaceNumbers(1)   = 3 
                orthogonalFaceNumbers(2)   = 4
                connectedSubCellIndexes(1) = 2
                connectedSubCellIndexes(2) = 4
                subCellIndexes(1,:)        = (1,2) 
                subCellIndexes(2,:)        = (2,2) 
            case (2)
                orthogonalFaceNumbers(1)   = 3 
                orthogonalFaceNumbers(2)   = 4
                connectedSubCellIndexes(1) = 1
                connectedSubCellIndexes(2) = 3
                subCellIndexes(1,:)        = (1,1) 
                subCellIndexes(2,:)        = (2,1) 
            case (3)
                orthogonalFaceNumbers(1)   = 1 
                orthogonalFaceNumbers(2)   = 2
                connectedSubCellIndexes(1) = 2
                connectedSubCellIndexes(2) = 1
                subCellIndexes(1,:)        = (1,2) 
                subCellIndexes(2,:)        = (1,1) 
            case (4)
                orthogonalFaceNumbers(1)   = 1 
                orthogonalFaceNumbers(2)   = 2
                connectedSubCellIndexes(1) = 4
                connectedSubCellIndexes(2) = 3
                subCellIndexes(1,:)        = (2,2) 
                subCellIndexes(2,:)        = (2,1) 
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
            call pr_FillCellFromRealCellData( this, centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ), &
                                                                                        parentCellDataBuffer )

            ! Location of current TrackCell relative to bigger cell
            ! defines indexes employed for filling buffers. 
            if ( &
                ( centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ) .eq.                               &
                  this%Grid%GetFaceConnection(                                                               & 
                    centerCellDataBuffer%GetFaceConnection( orthogonalFaceNumbers(1), 1 ), faceNumber, 1 ) ) &
            ) then 
                directConnectionSubRow        = 1
            else
                directConnectionSubRow        = 2
            end if 
            
            ! Fill Direct Connection
            call pr_FillCellBufferFromSubCell( this,                               &
                parentCellDataBuffer, subCellIndexes( directConnectionSubRow, 1 ), &
                subCellIndexes( directConnectionSubRow, 2 ), outputCellBuffer )

        else
            ! If equal size cell, fill buffer and leave
            call pr_FillCellFromRealCellData( this,                      &
                centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ), &
                outputCellBuffer )

        end if 

    else
        ! No connection, empty buffer
        call outputCellBuffer%Reset()

    end if 

    ! Done
    return


end subroutine pr_FillDirectBufferFromHorizontalFace





end module ParticleTrackingEngineModule


!***************************************************************************************************************
!
!***************************************************************************************************************
!
! Specifications
!---------------------------------------------------------------------------------------------------------------
!  implicit none
!---------------------------------------------------------------------------------------------------------------

! TRASH


! 17/04/2021

!function pr_FillCellWithInternalFlows( this, currentCellData, faceNumber ) result(filledCellData)
!    !
!    ! TODO
!    !
!    ! 
!
!    class(ParticleTrackingEngineType), target :: this
!    ! input
!    type(ModpathCellDataType), intent(in)  :: currentCellData
!    integer, intent(in)                    :: faceNumber
!    ! output
!    type(ModpathCellDataType) :: filledCellData
!    ! local 
!    integer :: m, id
!    type(ModpathCellDataType), dimension(8) :: cellDataBuffer
!
!    
!    ! It is understood that there is only one connection
!    call this%FillCellBuffer( currentCellData%GetFaceConnection( faceNumber, 1 ), filledCellData )
!
!    ! Force two rows and two columns 
!    ! to compute internal flows 
!    filledCellData%SubCellRowCount    = 2
!    filledCellData%SubCellColumnCount = 2
!
!    ! Doit
!    call filledCellData%ComputeSubCellFlows()
!
!    return 
!
!end function pr_FillCellWithInternalFlows









! 15/04/2021

    ! Once the directConnectionBuffer has been filled, 
    ! it is required to fill the neighbors of directConnectionBuffer
    ! but this process will be determined by the orientation 
    ! of requested faces, and the kind of process undergone to 
    ! fill the directConnectionBuffer.




    ! FOLLOWING LINES COULD BE REMOVEDDDD



    ! DO THE MAGIC
    
    ! DO IT YOU PIECE OF SHIT

        !! Check if more than one connection at the face
        !if ( this%TrackCell%CellData%SubFaceCounts(n) .gt. 1 ) then 
        !    ! If one than more connection with face n, then
        !    ! TrackCell is connected to smaller cells, create TrackCell
        !    ! with flows filled from surrounding cells 
        !    neighborCellData( cellCounter )     = pr_FillCellFromRefinedCells( this, this%TrackCell%CellData, n )
        !    neighborCellData( cellCounter + 1 ) = pr_FillCellFromRefinedCells( this, neighborCellData( cellCounter ), &
        !                                                                                           neighborFaceNumber )
        !    neighborCellData( cellCounter + 2 ) = pr_FillCellFromRefinedCells( this, neighborCellData( cellCounter ), &
        !                                                                                       neighborFaceNumber + 1 )


        !else if ( this%TrackCell%CellData%GetFaceConnection(n,1) .gt. 0 ) then


        !    ! HORIZONTAL FACES !!!!!!

        !    

        !    ! If cell has only one connection through face n, 
        !    ! then verify if connection is bigger or 
        !    ! the same size than current TrackCell 
        !    
        !    ! For a given faceNumber...
        !    ! orthogonalFaceNumbers:   faces orthogonal to current faceNumber. First value always related 
        !    !                          to the smallest faceNumber. Orthogonal in the x-y plane.
        !    ! connectedSubCellIndexes: sub cell indexes to which the cell could be connected 
        !    !                          through faceNumber
        !    ! subCellIndexes:          define indexes in bigger trackCell of direct connection to subcell, depends 
        !    !                          on relative location of current TrackCell with respect to bigger cell.
        !    !                          If connection to bigger cell is shared with lowerFaceNumber track cell, then 
        !    !                          direct connection through faceNumber n is with subcell (1,:) of bigger trackcell.
        !    !                          Indexes in parentCellDataBuffer
        !    select case ( n )
        !        case (1)
        !            orthogonalFaceNumbers(1)   = 3 
        !            orthogonalFaceNumbers(2)   = 4
        !            connectedSubCellIndexes(1) = 2
        !            connectedSubCellIndexes(2) = 4
        !            subCellIndexes(1,:)        = (1,2) 
        !            subCellIndexes(2,:)        = (2,2) 
        !        case (2)
        !            orthogonalFaceNumbers(1)   = 3 
        !            orthogonalFaceNumbers(2)   = 4
        !            connectedSubCellIndexes(1) = 1
        !            connectedSubCellIndexes(2) = 3
        !            subCellIndexes(1,:)        = (1,1) 
        !            subCellIndexes(2,:)        = (2,1) 
        !        case (3)
        !            orthogonalFaceNumbers(1)   = 1 
        !            orthogonalFaceNumbers(2)   = 2
        !            connectedSubCellIndexes(1) = 2
        !            connectedSubCellIndexes(2) = 1
        !            subCellIndexes(1,:)        = (1,2) 
        !            subCellIndexes(2,:)        = (1,1) 
        !        case (4)
        !            orthogonalFaceNumbers(1)   = 1 
        !            orthogonalFaceNumbers(2)   = 2
        !            connectedSubCellIndexes(1) = 4
        !            connectedSubCellIndexes(2) = 3
        !            subCellIndexes(1,:)        = (2,2) 
        !            subCellIndexes(2,:)        = (2,1) 
        !    end select 


        !    ! Analyze if current TrackCell is 
        !    ! connected to a bigger TrackCell. 
        !    ! Condition verifies if connections trough faceNumber (n)
        !    ! of cells connected through orthogonal faces to current TrackCell, 
        !    ! are the same or different. If the same, then current TrackCell is 
        !    ! connected to a bigger cell. 
        !    ! Assumes smoothed usg grid.
        !    if ( &
        !        ( this%TrackCell%CellData%GetFaceConnection(n,1) .eq.                                    &
        !          this%Grid%GetFaceConnection(                                                           & 
        !              this%TrackCell%CellData%GetFaceConnection( orthogonalFaceNumbers(1), 1 ), n, 1 ) ) &
        !        .or.                                                                                     &
        !        ( this%TrackCell%CellData%GetFaceConnection(n,1) .eq.                                    &
        !          this%Grid%GetFaceConnection(                                                           & 
        !              this%TrackCell%CellData%GetFaceConnection( orthogonalFaceNumbers(2), 1 ), n, 1 ) ) &
        !    ) then 
        !        ! Current TrackCell is connected to a bigger cell

        !        ! Fill parentCellDataBuffer 
        !        ! And force two rows and two columns 
        !        ! to compute internal flows 
        !        call this%FillCellBuffer( this%TrackCell%CellData%GetFaceConnection(n,1), parentCellDataBuffer )
        !        parentCellDataBuffer%SubCellRowCount    = 2
        !        parentCellDataBuffer%SubCellColumnCount = 2
        !        call parentCellDataBuffer%ComputeSubCellFlows()

        !        ! Note: Following if branch is analog for both 
        !        ! scenarios, only changing the indexes over which 
        !        ! operates.

        !        ! Location of current TrackCell relative to bigger cell,
        !        ! determines location of new buffers on buffer array
        !        if (
        !            ( this%TrackCell%CellData%GetFaceConnection(n,1) .eq.                                    &
        !              this%Grid%GetFaceConnection(                                                           & 
        !                  this%TrackCell%CellData%GetFaceConnection( orthogonalFaceNumbers(1), 1 ), n, 1 ) ) &
        !        ) then 
        !            ! Both current TrackCell and cell connected through first orthogonal 
        !            ! faceNumber are connected to the same cell. Both 
        !            ! buffers are filled and assigned to the buffer array

        !            ! Populate cell from subcell
        !            ! Direct connection
        !            call pr_FillCellBufferFromSubCell( this,                                   &
        !                parentCellDataBuffer, subCellIndexes( 1, 1 ) , subCellIndexes( 1, 2 ), &
        !                neighborCellData( cellCounter ) )

        !            ! Assign ParentCellDataBuffer, required ?
        !            neighborCellData( cellCounter )%ParentCellDataBuffer => parentCellDataBuffer 

        !            ! Populate cell from subcell
        !            ! This buffer corresponds to neighborFaceNumber 
        !            call pr_FillCellBufferFromSubCell( this,                                   &
        !                parentCellDataBuffer, subCellIndexes( 2, 1 ) , subCellIndexes( 2, 2 ), &
        !                neighborCellData( cellCounter + 1 ) )

        !            ! Populate remaining neighbor by analyzing the parent connections.
        !            if ( parentCellDataBuffer%SubFaceCounts( neighborFaceNumber + 1 ) .eq. 1 ) then
        !                ! If parent connection through neighborFaceNumber + 1 is single,
        !                ! then it means big cell connection should undergo
        !                ! double refinement.

        !                ! Fill aux buffer with trackcell connected through
        !                ! neighborFaceNumber + 1, to parentCellDataBuffer  
        !                call this%FillCellBuffer( &
        !                    parentCellDataBuffer%GetFaceConnection( neighborFaceNumber + 1, 1 ), auxCellDataBuffer )
        !                auxCellDataBuffer%SubCellRowCount    = 2
        !                auxCellDataBuffer%SubCellColumnCount = 2
        !                call auxCellDataBuffer%ComputeSubCellFlows()

        !                ! Populate cell from subcell
        !                ! This buffer corresponds to neighborFaceNumber + 1 
        !                call pr_FillCellBufferFromSubCell( this,                                &
        !                    auxCellDataBuffer, subCellIndexes( 2, 1 ) , subCellIndexes( 2, 2 ), &
        !                    neighborCellData( cellCounter + 2 ) )

        !            else if ( parentCellDataBuffer%SubFaceCounts( neighborFaceNumber + 1 ) .gt. 1 ) then  
        !                ! If the connection is multiple, then it should fill the 
        !                ! buffer with the cell connected through 
        !                ! neighborFaceNumber + 1, using column value obtained from 
        !                ! subCellIndexes. Assumes smoothed usg grid. 
        !                call this%FillCellBuffer( &
        !                    parentCellDataBuffer%GetFaceConnection( neighborFaceNumber + 1, subCellIndexes( 1, 2 ) ), &
        !                                                                          neighborCellData( cellCounter + 2 ) )
        !                neighborCellData( cellCounter + 2 )%SubCellRowCount    = 2
        !                neighborCellData( cellCounter + 2 )%SubCellColumnCount = 2
        !                call neighborCellData( cellCounter + 2 )%ComputeSubCellFlows()

        !            else
        !                ! No connection through neighborFaceNumber + 1, empty buffer
        !                call neighborCellData( cellCounter + 2 )%Reset()

        !            end if 


        !        else
        !            ! Cell is connected to a bigger cell, shares connection 
        !            ! with trackCell connected through the second 
        !            ! orthogonal faceNumber

        !            ! Fill cell with data from subCell related 
        !            ! to second orthogonal faceNumber
        !            ! Direct connection
        !            call pr_FillCellBufferFromSubCell( this,                                   &
        !                parentCellDataBuffer, subCellIndexes( 2, 1 ) , subCellIndexes( 2, 2 ), &
        !                neighborCellData( cellCounter ) )

        !            ! Assign ParentCellDataBuffer
        !            neighborCellData( cellCounter )%ParentCellDataBuffer => parentCellDataBuffer 

        !            ! This corresponds to neighborFaceNumber + 1 
        !            call pr_FillCellBufferFromSubCell( this,                                   &
        !                parentCellDataBuffer, subCellIndexes( 1, 1 ) , subCellIndexes( 1, 2 ), &
        !                neighborCellData( cellCounter + 2 ) )
        !    

        !            ! Populate remaining neighbor by analyzing the parent connections
        !            if ( parentCellDataBuffer%SubFaceCounts( neighborFaceNumber ) .eq. 1 ) then

        !                ! Fill aux buffer with trackcell connected through 
        !                ! neighborFaceNumber, to parentCellDataBuffer
        !                call this%FillCellBuffer( &
        !                    parentCellDataBuffer%GetFaceConnection( neighborFaceNumber, 1 ), auxCellDataBuffer )
        !                auxCellDataBuffer%SubCellRowCount    = 2
        !                auxCellDataBuffer%SubCellColumnCount = 2
        !                call auxCellDataBuffer%ComputeSubCellFlows()

        !                ! Populate cell from subcell
        !                ! This buffer corresponds to neighborFaceNumber 
        !                call pr_FillCellBufferFromSubCell( this,                                &
        !                    auxCellDataBuffer, subCellIndexes( 1, 1 ) , subCellIndexes( 1, 2 ), &
        !                    neighborCellData( cellCounter + 1 ) )

        !            else if ( parentCellDataBuffer%SubFaceCounts( neighborFaceNumber ) .gt. 1 ) then  
        !                ! If the connection is multiple, then it should fill the 
        !                ! buffer with the cell connected through 
        !                ! neighborFaceNumber, using column value obtained from subCellIndexes 
        !                call this%FillCellBuffer( &
        !                    parentCellDataBuffer%GetFaceConnection( neighborFaceNumber, subCellIndexes( 2, 2 ) ), &
        !                                                                      neighborCellData( cellCounter + 1 ) )
        !                neighborCellData( cellCounter + 1 )%SubCellRowCount    = 2
        !                neighborCellData( cellCounter + 1 )%SubCellColumnCount = 2
        !                call neighborCellData( cellCounter + 1 )%ComputeSubCellFlows()

        !            else
        !                ! No connection through neighborFaceNumber, empty buffer
        !                call neighborCellData( cellCounter + 1 )%Reset()

        !            end if 

        !        end if 


        !    else
        !        ! Current TrackCell is connected to a cell of equal size
        !        ! then FillCellBuffer using connection and force subcellflows
        !        neighborCellData( cellCounter )     = pr_FillCellWithInternalFlows( this, this%TrackCell%CellData, n )

        !        ! It seems that these functions should go through the number of 
        !        ! connections detection process
        !        neighborCellData( cellCounter + 1 ) = pr_FillCellWithInternalFlows( this, neighborCellData( cellCounter ), &
        !                                                                                               neighborFaceNumber )
        !        neighborCellData( cellCounter + 2 ) = pr_FillCellWithInternalFlows( this, neighborCellData( cellCounter ), &
        !                                                                                           neighborFaceNumber + 1 )

        !    end if 


        !else

        !    ! No connection, 
        !    ! Reset TrackCell buffers
        !    call neighborCellData( cellCounter )%Reset()
        !    call neighborCellData( cellCounter + 1 )%Reset()
        !    call neighborCellData( cellCounter + 2 )%Reset()
        !    
        !end if 








! 14/04/2021


    ! CLEAN YOU MTF
    !! If requesting connections through vertical face 
    !! then, assuming no vertical refinement, implies that 
    !! cellCounter can be initialized immediatelly with 
    !! connected cell id. 
    !if ( ( n .eq. 5 ) .or. ( n .eq. 6 ) ) then 

    !    ! This buffer should be taken as the center for the 
    !    ! initialization of neighbors
    !    call this%FillCellBuffer( &
    !        this%TrackCell%CellData%GetFaceConnection(n,1), neighborCellData( cellCounter ) )

    !    ! Once the direct connection buffer is filled
    !    ! Fill neighbor cells connected through
    !    ! neighborFaceNumber and neighborFaceNumber + 1 

    !    ! Former could be horizontal faces, and in such case, 
    !    ! method should undergo process of horizontal verification 
    !    ! before filling respective buffers 
    !
    !end if 



    ! OR MORE LIKE
    ! FillCellFromFaceConnection( this, centerCellDataBuffer, faceNumber,  )



        ! Once the direct connection buffer is filled
        ! Fill neighbor cells connected through
        ! firstNeighborFaceNumber and firstNeighborFaceNumber + 1 

        ! Former could be horizontal faces, and in such case, 
        ! method should undergo process of horizontal verification 
        ! before filling respective buffers 


    ! Fill direct connection buffer 
    !neighborCellsDataBuffer( directConnectionBufferIndex )
    

    ! Fill indirect connections buffer  
    
    ! In each filling process, it is required to 
    ! identify the direction (axis) of the requested
    ! face. 


    




        ! Thus at this point it should be enough to do 


        ! PROTOTYPE
        !neighborCellData( cellCounter )     = pr_FillNeighborCellBuffer( this, this%TrackCell%CellData, n )
        !neighborCellData( cellCounter + 1 ) = pr_FillNeighborCellBuffer( this, neighborCellData( cellCounter ), &
        !                                                                       neighborFaceNumber )
        !neighborCellData( cellCounter + 2 ) = pr_FillNeighborCellBuffer( this, neighborCellData( cellCounter ), &
        !                                                                                 neighborFaceNumber + 1 )

        ! Although if the cell is connected to a bigger cell, 
        ! the function should be able to handler the situation 
        ! of filling the corresponding buffers immediatelly, 
        ! that is, the prototype function should be something 
        ! more like

        !print * , 'THIS IS A FAKE CELL'
        !print * , 'It is composed by other trackcells with ids'
        !print * , currentCellData%SubCellIds(1), currentCellData%SubCellIds(2), &
        !    currentCellData%SubCellIds(3), currentCellData%SubCellIds(4)


        !print *, 'TRY TO SEE ABOUT POINTER'
        !print *, currentCellData%SubCellDataBuffer(1)%CellNumber
        !print *, currentCellData%SubCellDataBuffer(2)%CellNumber
        !print *, currentCellData%SubCellDataBuffer(3)%CellNumber
        !print *, currentCellData%SubCellDataBuffer(4)%CellNumber
        

        !! When this happens FillCellBuffer wont do the work as 
        !! cell comes with a fake id

        !! Still it is required to fill cellDataBuffer

        !! Depending on faceNumber, orientation
        !! or indexes to extract cell ids for 
        !! creating sample cell will be 
        !! different

        !! e.g. - when querying faceNumber 2 from currentCellData, 
        !!      connected cells through that face would follow same
        !!      subindexes from new/sample cell.
        !!      - when querying faceNumber 4, first connected subcell
        !!      would occupy subindex 21 from the new/sample cell. 

        !
        !! There should be a loop through SubCellIds to 
        !! fill cellDataBuffer

        !! Subcells (actually cells, but subcells from the perspective of the sample cell)
        !! have connection information

        !! Buffer should be filled
        !! if faceNumber 4, then 41 is subcell 3 (21) and 42 is subcell 4 (22),
        !! e (1)face of the former subcells.




    !print *, ' Cell: ', currentCellData%CellNumber, ' faceNumber: ', faceNumber, 'GFC 1 ', &
    !    currentCellData%GetFaceConnection( faceNumber, 1 ), ' GFC 2', currentCellData%GetFaceConnection( faceNumber, 2 )

    ! Question is if filledCellData requires
    ! connection information

    ! SubFaceCounts and SubFaceConns are required

    ! Fill connections
    !
    ! Not generalized for vertical faces 
    ! 
    !!filledCellData%SubFaceConn1( 1 ) = cellDataBuffer( 1 )%GetFaceConnection( 1, 1 ) 
    !!filledCellData%SubFaceConn1( 2 ) = cellDataBuffer( 3 )%GetFaceConnection( 1, 1 )
    !!if ( filledCellData%SubFaceConn1( 1 ) .eq.  filledCellData%SubFaceConn1( 2 ) ) then
    !!    filledCellData%SubFaceCounts( 1 ) = 1
    !! else
    !!    filledCellData%SubFaceCounts( 1 ) = 2
    !!end if

    !!filledCellData%SubFaceConn2( 1 ) = cellDataBuffer( 2 )%CellNumber 
    !!filledCellData%SubFaceConn2( 2 ) = cellDataBuffer( 4 )%CellNumber
    !!if (                                    &
    !!  filledCellData%SubFaceConn2( 1 ) .eq. &
    !!  filledCellData%SubFaceConn2( 2 ) ) then
    !!    filledCellData%SubFaceCounts( 2 ) = 1
    !! else
    !!    filledCellData%SubFaceCounts( 2 ) = 2
    !!end if


    !!filledCellData%SubFaceConn3( 1 ) = cellDataBuffer( 3 )%CellNumber 
    !!filledCellData%SubFaceConn3( 2 ) = cellDataBuffer( 2 )%CellNumber
    !!if (                                    &
    !!  filledCellData%SubFaceConn3( 1 ) .eq. &
    !!  filledCellData%SubFaceConn3( 2 ) ) then
    !!    filledCellData%SubFaceCounts( 3 ) = 1
    !! else
    !!    filledCellData%SubFaceCounts( 3 ) = 2
    !!end if

    !!filledCellData%SubFaceConn4( 1 ) = 
    !!filledCellData%SubFaceConn4( 2 ) = 
    !!if (                                    &
    !!  filledCellData%SubFaceConn4( 1 ) .eq. &
    !!  filledCellData%SubFaceConn4( 2 ) ) then
    !!    filledCellData%SubFaceCounts( 4 ) = 1
    !! else
    !!    filledCellData%SubFaceCounts( 4 ) = 2
    !!end if

    !!filledCellData%SubFaceConn5( 1 ) = 
    !!filledCellData%SubFaceConn5( 2 ) = 
    !!filledCellData%SubFaceConn5( 3 ) = 
    !!filledCellData%SubFaceConn5( 4 ) = 

    !!filledCellData%SubFaceConn6( 1 ) = 
    !!filledCellData%SubFaceConn6( 2 ) = 
    !!filledCellData%SubFaceConn6( 3 ) = 
    !!filledCellData%SubFaceConn6( 4 ) = 


    !allocate( filledCellData%SubCellDataBuffer(4) )

    !filledCellData%SubCellDataBuffer(1) => cellDataBuffer(1)
    !filledCellData%SubCellDataBuffer(2) => cellDataBuffer(2)
    !filledCellData%SubCellDataBuffer(3) => cellDataBuffer(3)
    !filledCellData%SubCellDataBuffer(4) => cellDataBuffer(4)


    !print *, 'IM PRINTING INFO FROM THE POINTER'
    !print *, filledCellData%SubCellData%CellNumber

    !! The property way
    !filledCellData%SubCellDataBuffer(1) = cellDataBuffer(1)
    !filledCellData%SubCellDataBuffer(2) = cellDataBuffer(2)
    !filledCellData%SubCellDataBuffer(3) = cellDataBuffer(3)
    !filledCellData%SubCellDataBuffer(4) = cellDataBuffer(4)

    !! Fill Flows

    !! Intercell 
    !filledCellData%SubCellFlows(1) = cellDataBuffer(1)%GetAveragedFaceFlow(2)
    !filledCellData%SubCellFlows(2) = cellDataBuffer(3)%GetAveragedFaceFlow(2)
    !filledCellData%SubCellFlows(3) = cellDataBuffer(3)%GetAveragedFaceFlow(4)
    !filledCellData%SubCellFlows(4) = cellDataBuffer(4)%GetAveragedFaceFlow(4) 

    !! Horizontal faces 
    !filledCellData%Q1(1) = cellDataBuffer(1)%GetAveragedFaceFlow(1)
    !filledCellData%Q1(2) = cellDataBuffer(3)%GetAveragedFaceFlow(1)
    !filledCellData%Q2(1) = cellDataBuffer(2)%GetAveragedFaceFlow(2)
    !filledCellData%Q2(2) = cellDataBuffer(4)%GetAveragedFaceFlow(2)
    !filledCellData%Q3(1) = cellDataBuffer(3)%GetAveragedFaceFlow(3)
    !filledCellData%Q3(2) = cellDataBuffer(4)%GetAveragedFaceFlow(3)
    !filledCellData%Q4(1) = cellDataBuffer(1)%GetAveragedFaceFlow(4)
    !filledCellData%Q4(2) = cellDataBuffer(2)%GetAveragedFaceFlow(4)

    !! Vertical faces
    !filledCellData%Q5(1) = cellDataBuffer(1)%GetAveragedFaceFlow(5)
    !filledCellData%Q5(2) = cellDataBuffer(2)%GetAveragedFaceFlow(5)
    !filledCellData%Q5(3) = cellDataBuffer(3)%GetAveragedFaceFlow(5)
    !filledCellData%Q5(4) = cellDataBuffer(4)%GetAveragedFaceFlow(5)
    !filledCellData%Q6(1) = cellDataBuffer(1)%GetAveragedFaceFlow(6)
    !filledCellData%Q6(2) = cellDataBuffer(2)%GetAveragedFaceFlow(6)
    !filledCellData%Q6(3) = cellDataBuffer(3)%GetAveragedFaceFlow(6)
    !filledCellData%Q6(4) = cellDataBuffer(4)%GetAveragedFaceFlow(6)


    !            !! Sources/sinks and additional shit
    !            ! porosity, ibound, etc


! 13/04/2021


                        !! Populate faceFlows with complementary subcell indexes
                        !call auxCellDataBuffer%FillSubCellFaceFlowsBuffer( &
                        !    subCellIndexes( 1, 1 ) , subCellIndexes( 1, 2 ) , faceFlows )
                        !
                        !! Given faceFlows, fill trackcell
                        !! REMEMBER THAT THIS REQUIRES FILLING SOURCES AND SINKS
                        !neighborCellData( cellCounter + 1 )%Q1(1) = faceFlows(1)
                        !neighborCellData( cellCounter + 1 )%Q2(1) = faceFlows(2)
                        !neighborCellData( cellCounter + 1 )%Q3(1) = faceFlows(3)
                        !neighborCellData( cellCounter + 1 )%Q4(1) = faceFlows(4)
                        !neighborCellData( cellCounter + 1 )%Q5(1) = faceFlows(5)
                        !neighborCellData( cellCounter + 1 )%Q6(1) = faceFlows(6)
                        !! REMEMBER THAT THIS REQUIRES FILLING SOURCES AND SINKS
                        !neighborCellData( cellCounter + 1 )%SubCellRowCount    = 2
                        !neighborCellData( cellCounter + 1 )%SubCellColumnCount = 2
                        !call neighborCellData( cellCounter + 1 )%ComputeSubCellFlows()

                    !! Fill faceFlows with data from subcellindexes 
                    !! related to the second orthogonal faceNumber
                    !call parentCellDataBuffer%FillSubCellFaceFlowsBuffer( &
                    !    subCellIndexes( 2, 1 ) , subCellIndexes( 2, 2 ) , faceFlows )
                    !
                    !! Given the face flows, fill direct connection TrackCell
                    !! REMEMBER THAT THIS REQUIRES FILLING SOURCES AND SINKS
                    !neighborCellData( cellCounter )%Q1(1) = faceFlows(1)
                    !neighborCellData( cellCounter )%Q2(1) = faceFlows(2)
                    !neighborCellData( cellCounter )%Q3(1) = faceFlows(3)
                    !neighborCellData( cellCounter )%Q4(1) = faceFlows(4)
                    !neighborCellData( cellCounter )%Q5(1) = faceFlows(5)
                    !neighborCellData( cellCounter )%Q6(1) = faceFlows(6)
                    !! REMEMBER THAT THIS REQUIRES FILLING SOURCES AND SINKS
                    !neighborCellData( cellCounter )%SubCellRowCount    = 2
                    !neighborCellData( cellCounter )%SubCellColumnCount = 2
                    !call neighborCellData( cellCounter )%ComputeSubCellFlows()
                    !! Update faceFlows, using the other connected subcell
                    !call parentCellDataBuffer%FillSubCellFaceFlowsBuffer( &
                    !    subCellIndexes( 1, 1 ) , subCellIndexes( 1, 2 ) , faceFlows )
                    !! REMEMBER THAT THIS REQUIRES FILLING SOURCES AND SINKS
                    !neighborCellData( cellCounter + 2 )%Q1(1) = faceFlows(1)
                    !neighborCellData( cellCounter + 2 )%Q2(1) = faceFlows(2)
                    !neighborCellData( cellCounter + 2 )%Q3(1) = faceFlows(3)
                    !neighborCellData( cellCounter + 2 )%Q4(1) = faceFlows(4)
                    !neighborCellData( cellCounter + 2 )%Q5(1) = faceFlows(5)
                    !neighborCellData( cellCounter + 2 )%Q6(1) = faceFlows(6)
                    !! REMEMBER THAT THIS REQUIRES FILLING SOURCES AND SINKS
                    !neighborCellData( cellCounter + 2 )%SubCellRowCount    = 2
                    !neighborCellData( cellCounter + 2 )%SubCellColumnCount = 2
                    !call neighborCellData( cellCounter + 2 )%ComputeSubCellFlows()


                        !! Populate faceFlows with complementary subcell indexes
                        !call auxCellDataBuffer%FillSubCellFaceFlowsBuffer( &
                        !    subCellIndexes( 2, 1 ) , subCellIndexes( 2, 2 ) , faceFlows )
                        !
                        !! Given faceFlows, fill trackcell
                        !! REMEMBER THAT THIS REQUIRES FILLING SOURCES AND SINKS
                        !neighborCellData( cellCounter + 2 )%Q1(1) = faceFlows(1)
                        !neighborCellData( cellCounter + 2 )%Q2(1) = faceFlows(2)
                        !neighborCellData( cellCounter + 2 )%Q3(1) = faceFlows(3)
                        !neighborCellData( cellCounter + 2 )%Q4(1) = faceFlows(4)
                        !neighborCellData( cellCounter + 2 )%Q5(1) = faceFlows(5)
                        !neighborCellData( cellCounter + 2 )%Q6(1) = faceFlows(6)
                        !! REMEMBER THAT THIS REQUIRES FILLING SOURCES AND SINKS
                        !neighborCellData( cellCounter + 2 )%SubCellRowCount    = 2
                        !neighborCellData( cellCounter + 2 )%SubCellColumnCount = 2
                        !call neighborCellData( cellCounter + 2 )%ComputeSubCellFlows()


                    !! This buffer corresponds to neighborFaceNumber 
                    !call parentCellDataBuffer%FillSubCellFaceFlowsBuffer( &
                    !    subCellIndexes( 2, 1 ) , subCellIndexes( 2, 2 ) , faceFlows )
                    !! REMEMBER THAT THIS REQUIRES FILLING SOURCES AND SINKS
                    !neighborCellData( cellCounter + 1 )%Q1(1) = faceFlows(1)
                    !neighborCellData( cellCounter + 1 )%Q2(1) = faceFlows(2)
                    !neighborCellData( cellCounter + 1 )%Q3(1) = faceFlows(3)
                    !neighborCellData( cellCounter + 1 )%Q4(1) = faceFlows(4)
                    !neighborCellData( cellCounter + 1 )%Q5(1) = faceFlows(5)
                    !neighborCellData( cellCounter + 1 )%Q6(1) = faceFlows(6)
                    !! REMEMBER THAT THIS REQUIRES FILLING SOURCES AND SINKS
                    !neighborCellData( cellCounter + 1 )%SubCellRowCount    = 2
                    !neighborCellData( cellCounter + 1 )%SubCellColumnCount = 2
                    !call neighborCellData( cellCounter + 1 )%ComputeSubCellFlows()



                    !! Populate faceFlows, using subCell from bigger TrackCell
                    !! located at subCellIndexes(1,:)
                    !call parentCellDataBuffer%FillSubCellFaceFlowsBuffer( &
                    !    subCellIndexes( 1, 1 ) , subCellIndexes( 1, 2 ) , faceFlows )

                    !! Given the face flows, fill buffer faceflows
                    !! and sources/sinks/storage
                    !neighborCellData( cellCounter )%Q1(1)       = faceFlows(1)
                    !neighborCellData( cellCounter )%Q2(1)       = faceFlows(2)
                    !neighborCellData( cellCounter )%Q3(1)       = faceFlows(3)
                    !neighborCellData( cellCounter )%Q4(1)       = faceFlows(4)
                    !neighborCellData( cellCounter )%Q5(1)       = faceFlows(5)
                    !neighborCellData( cellCounter )%Q6(1)       = faceFlows(6)
                    !neighborCellData( cellCounter )%SourceFlow  = parentCellDataBuffer%SourceFlow / 4d0
                    !neighborCellData( cellCounter )%SinkFlow    = parentCellDataBuffer%SinkFlow / 4d0
                    !neighborCellData( cellCounter )%StorageFlow = parentCellDataBuffer%StorageFlow / 4d0
                    !neighborCellData( cellCounter )%SubCellRowCount    = 2
                    !neighborCellData( cellCounter )%SubCellColumnCount = 2
                    !call neighborCellData( cellCounter )%ComputeSubCellFlows()
                    !! NOTE: COMPUTESUBCELLFLOWS considers that sources/sinks 
                    !! at TrackCell are divided by 4, and then solves the balance
                    !! When filling a trackcell buffer, this means that 
                    !! it is required to just assign sources/sinks/storage and that 
                    !! should be enough for subcellflows computation.
                    !! The thing is though, that such sources/sinks/storage 
                    !! should actually be 1/4 of the total flow in the parent cell subcell




                    ! This conditional branch determines immediatelly 
                    ! two of the neighbor buffers 

                    ! If the direct connection is on the upper row, 
                    ! then it means that this process immediatelly returns 
                    ! the trackcell buffer connected through the  

                    ! FILL IT FROM SUBCELL


                    ! The current TrackCell is directly connected through face 
                    ! n with subcell, to be defined.
                    ! Idea is to connect the corresponding subcell

                    ! Function that populates TrackCell
                    ! Get flows related to a given subcell
                    !call this%TrackCell%CellData%FillSubCellFaceFlowsBuffer( subRow, subColumn, faceFlows )
                    ! Fill direct connection TrackCell with values from subcellindex 
                    ! obtained from connectedSubCellIndexes 1

            ! Is the cell connected to a bigger cell through face n ?
            ! Get the connection through face n.
            ! Compare it with connections through face n 
            ! of surrounding cells connected to center trackcell
            ! in direction perpendicular to requested face n.


        !else if ( this%TrackCell%CellData%SubFaceCounts(n) .eq. 1 ) then 
            ! Note: it is possible for cell data to have subFaceCounts == 1 and 
            ! GetFaceConnection == 0, thus it is a more restrictive 
            ! condition the one obtained from GetFaceConnection

    !else
    !    print *, '** ParticleTrackingEngine.FillNeighborCellData: found a normal single cell, it has ', &
    !    this%TrackCell%CellData%GetSubCellCount(), ' subcells.'
    !    print *, '** ParticleTrackingEngine.FillNeighborCellData: will try to modify internal values of column and row count'

    !    ! If it is single or multiple, in both 
    !    ! cases requires initialization of neighbor cells 
    !    ! connected to current TrackCell

    !    ! It seems that being a single cell only affects in how to 
    !    ! initialize the current track cell, and that is it. 
    !    ! This means that if the current TrackCell is not refined, 
    !    ! essentially the process of cells initialization would 
    !    ! be the same. 

    !    ! If the current cell is not refined, it does not have refined 
    !    ! direct neighbors but it can have diagonally refined cells 


    !    ! It could be possible to consider that although 
    !    ! the center cell is not refined, all of its neighbors 
    !    ! will be refined regardless. 

    !end if


! OLDER ThAN 13/04/2021

        !! OBS: setDataUnstructured fill values based on SubCellCounts, which is done 
        !! before even considering filling the values of counts.


        !! This is possible and yields intercell flows
        !neighborCellData(1)%SubCellRowCount    = 2
        !neighborCellData(1)%SubCellColumnCount = 2
        !call neighborCellData(1)%ComputeSubCellFlows()


        ! Thus when forced as before, flow rates through main 
        ! faces remain as the original values, but with intercell fluxes
        ! and it would be necessary to use GetSubCellBoundaryFlow in order to 
        ! get the proper divided value or divided it manually. 

        
        ! For this case it seems enough to just initialize 
        ! the cell and sub cell buffers
        ! If this occurs once, then all neighbor cell will 
        ! be the same, same refinement 


                ! Fill flows at the corresponding sampling track cell
                !! Which are the source cells ?
                !do m = 1, this%TrackCell%CellData%SubFaceCounts(n)

                !    ! Compute cell id in buffer and reset
                !    id = (m-1)*this%TrackCell%CellData%SubFaceCounts(n) + 1
                !    call neighborCellBuffer(   id   )%Reset()
                !    call neighborCellBuffer( id + 1 )%Reset()

                !    ! This should work always
                !    ! Fill buffer for connected cell
                !    call this%FillCellBuffer( this%TrackCell%CellData%GetFaceConnection(n,m), neighborCellBuffer( id ) )

                !    ! And then fill the buffer for the connection of the connected cell

                !    ! TODO
                !    ! To handle the case of connection with multiple cells
                !    ! cell its initialized with the same data as before and 
                !    ! relevant flows are updated with averaged values
                !    ! It seems better to apply a derefinement, given a set of cells
                !    ! but treat it like a specific case
                !    !call this%FillCellBuffer( this%TrackCell%CellData%GetFaceConnection(n,m), neighborCellBuffer( id + 1 ) )

                !    ! For the simple unstructured this should do it 
                !    call this%FillCellBuffer( neighborCellBuffer( id )%GetFaceConnection(n,1), neighborCellBuffer( id + 1 ) )
                !    
                !    print *, 'Buffer id has cell :', neighborCellBuffer( id )%CellNumber
                !    print *, 'Buffer id+1 has cell :', neighborCellBuffer( id + 1 )%CellNumber

                !end do


                !! Once information is known, fill the sampling track cell

                !! Reset the sampling cell holder 
                !call neighborCellData( cellCounter )%Reset()
                !

                !! Buffer cells (fake subcells) transfer information to 
                !! sampling track cell

                !! Do I need connections information ?
                !! It may be useful to handle creation of 
                !! diagonally located neighbor cells 
                !! This immediatelly allows detection if 
                !! in the diagonally located cell needs 
                !! DeRefinement

                !! Fill connections
                !!neighborCellData( cellCounter )%SubFaceConn1(1) = 
                !!neighborCellData( cellCounter )%SubFaceConn1(2) = 
                !!neighborCellData( cellCounter )%SubFaceConn2(1) = 
                !!neighborCellData( cellCounter )%SubFaceConn2(2) = 
                !!neighborCellData( cellCounter )%SubFaceConn3(1) = 
                !!neighborCellData( cellCounter )%SubFaceConn3(2) = 
                !!neighborCellData( cellCounter )%SubFaceConn4(1) = 
                !!neighborCellData( cellCounter )%SubFaceConn4(2) = 
                !!neighborCellData( cellCounter )%SubFaceConn5(1) = 
                !!neighborCellData( cellCounter )%SubFaceConn5(2) = 
                !!neighborCellData( cellCounter )%SubFaceConn5(3) = 
                !!neighborCellData( cellCounter )%SubFaceConn5(4) = 
                !!neighborCellData( cellCounter )%SubFaceConn6(1) = 
                !!neighborCellData( cellCounter )%SubFaceConn6(2) = 
                !!neighborCellData( cellCounter )%SubFaceConn6(3) = 
                !!neighborCellData( cellCounter )%SubFaceConn6(4) = 

                ! 

                !! Fill Flows

                !! Intercell 
                !neighborCellData( cellCounter )%SubCellFlows(1) = neighborCellBuffer(1)%GetAveragedFaceFlow(2)
                !neighborCellData( cellCounter )%SubCellFlows(2) = neighborCellBuffer(3)%GetAveragedFaceFlow(2)
                !neighborCellData( cellCounter )%SubCellFlows(3) = neighborCellBuffer(3)%GetAveragedFaceFlow(4)
                !neighborCellData( cellCounter )%SubCellFlows(4) = neighborCellBuffer(4)%GetAveragedFaceFlow(4) 

                !! Horizontal faces 
                !neighborCellData( cellCounter )%Q1(1) = neighborCellBuffer(1)%GetAveragedFaceFlow(1)
                !neighborCellData( cellCounter )%Q1(2) = neighborCellBuffer(3)%GetAveragedFaceFlow(1)
                !neighborCellData( cellCounter )%Q2(1) = neighborCellBuffer(2)%GetAveragedFaceFlow(2)
                !neighborCellData( cellCounter )%Q2(2) = neighborCellBuffer(4)%GetAveragedFaceFlow(2)
                !neighborCellData( cellCounter )%Q3(1) = neighborCellBuffer(3)%GetAveragedFaceFlow(3)
                !neighborCellData( cellCounter )%Q3(2) = neighborCellBuffer(4)%GetAveragedFaceFlow(3)
                !neighborCellData( cellCounter )%Q4(1) = neighborCellBuffer(1)%GetAveragedFaceFlow(4)
                !neighborCellData( cellCounter )%Q4(2) = neighborCellBuffer(2)%GetAveragedFaceFlow(4)

                !! Vertical faces
                !neighborCellData( cellCounter )%Q5(1) = neighborCellBuffer(1)%GetAveragedFaceFlow(5)
                !neighborCellData( cellCounter )%Q5(2) = neighborCellBuffer(2)%GetAveragedFaceFlow(5)
                !neighborCellData( cellCounter )%Q5(3) = neighborCellBuffer(3)%GetAveragedFaceFlow(5)
                !neighborCellData( cellCounter )%Q5(4) = neighborCellBuffer(4)%GetAveragedFaceFlow(5)
                !neighborCellData( cellCounter )%Q6(1) = neighborCellBuffer(1)%GetAveragedFaceFlow(6)
                !neighborCellData( cellCounter )%Q6(2) = neighborCellBuffer(2)%GetAveragedFaceFlow(6)
                !neighborCellData( cellCounter )%Q6(3) = neighborCellBuffer(3)%GetAveragedFaceFlow(6)
                !neighborCellData( cellCounter )%Q6(4) = neighborCellBuffer(4)%GetAveragedFaceFlow(6)



!!!!! 08/04/2021

            ! The thing is, now that I have filled information 
            ! for sampling cells connected through trackcell faces
            ! what to do about diagonal cells ?

            ! Depending on the number of connections
            ! of a connected cell through the face, 
            ! it will determine if the diagonal 
            ! should be filled as a normal cell buffer 
            ! or a DeRefinement procedure.        

            ! Lets take an specific case 

            
            ! For example face 2
            ! In previous implementation it was required for this connection 
            ! to also return information about its connections 3 and 4



    ! SPECIFICATION FOR FILL CELL BUFFER
    !pr_FillCellBuffer(this, cellNumber, cellBuffer)
    !if(cellBuffer%GetSubCellCount() .gt. 1) then
    !    call cellBuffer%ComputeSubCellFlows()
    !end if
    ! call this%FillCellBuffer(loc%CellNumber,  this%TrackCell%CellData)

    ! If unstructured, and the current TrackCell has subCells, then 
    ! compute subCell internal flows for all neighbor cells.

    ! So, something is needed to loop properly 
    ! over neighbor cells and fill corresponding buffers 
    ! The issue seems to be how to force the computation of 
    ! internal flows.

    ! The routine FillCellBuffer invokes the mpcelldata function 
    ! ComputeSubCellFlows. 

    ! Such function performs a check to detect if GetSubCellCount() == 1
    ! If that is the case, then the function leaves
    ! Otherwise, compute sub cell flows

    ! GetSubCellCount returns the product between subCell columns 
    ! and subCell rows

    ! Function mpcelldata%SetDataUnstructured is the one that 
    ! assigns values for the number of subcell columns and sub
    ! cell rows.


    ! The question is, is it possible to compute internal fluxes 
    ! even in the case in which the TrackCell has no subcells from 
    ! the default data loading from mp7.


    ! Note that in mpcelldata%SetDataUnstructured if any of the 
    ! connected faces has SubFaceCount higher than 1, 
    ! then it is assumed/understood that the cell has subcells in the 
    ! format nrows and ncols = 2.


    ! At ModpathCellDataType 
    !this%SubCellRowCount = 2
    !this%SubCellColumnCount = 2
    
    ! Computation of subcellflows extracts data using 
    ! GetSubCellBoundaryFlow. This last function 
    ! computes the flow through the boundary. 
    ! If the face has only one cell connection 
    ! in a given face, then the function GetSubCellBoundaryFlow
    ! computes that the flow through such face is equal to half 
    ! of the total flow.

    ! In the case that a neigbor cell is already multi sub cell, 
    ! then this function will access to the corresponding boundary flows. 
    ! The latter is verified by looking into the array SubFaceCounts 
    
    ! So, forcing computation of internal flows should be done only for those 
    ! neighbors that are not already multi sub cell.

    ! Internal flows are saved into the array SubCellFlows, 4 elements
    ! at mpcelldata
  
    ! So, when looping to fill neighborcelldata
    ! it will be necessary to verify that a given neighbor cell 
    ! has or not subcells. If it has, then no problem and move on to populate. 
    ! If it has not subcells, then force them by modifying the parameters 
    ! subcellrowcount and subcellcolumncount. Product between these two will 
    ! give the getsubcellcount gt 1, then subcell flows will be computed 
    ! at setdataunstructured.


    ! The thing with unstructured grid is how to extract the corresponding cell 
    ! number when next to a multi cell. That is, there will be several cells connected 
    ! through a face. Smaller cells, but cells at the end. 


    ! Should we be build a fake neighbor cell that replicates neighbor multi cells 
    ! as a single trackcell ? 

    ! Loop to fill neighbor cell data
    ! Remember that a single cell will 
    ! also provide information about surrounding 
    ! corner cells.  


        ! This is the case when there are sub cells.
        ! This configuration sets the required refinement
        ! and neighbor cells.

        ! Number of subcells at TrackCell, define
        ! which refinement level should have 
        ! neighbor cells.

        ! ... and also which cells are considered neighbors.

        ! In modpath docs (Fig 4b, p. 8), its seen a case in which a 
        ! refined trackcell is in contact to the left with
        ! a coarser cell (smaller refinement) and immediately
        ! to the right with a finer cell (higher refinement).

        ! This means that for interpolation purposes
        ! it will be necessary to apply a "double" refinement
        ! to the coarse cell.

        ! First refinement return 4 subcells and  
        ! subcells 12, 22 are refined again to reach the
        ! refinement level of the current cell.

        ! TRYING TO SEE WHAT WE GET  

        !this%TrackCell%CellData%SubCellRowCount    = 2
        !this%TrackCell%CellData%SubCellColumnCount = 2
        !print *, '** ParticleTrackingEngine.FillNeighborCellData: after modifying internal values... ', &
        !this%TrackCell%CellData%GetSubCellCount(), ' subcells.'

        !!call this%FillCellBuffer( this%TrackCell%CellData%CellNumber, this%TrackCell%CellData ) 
        !call this%FillCellBuffer( this%TrackCell%CellData%CellNumber , neighborCellData(1) ) 
