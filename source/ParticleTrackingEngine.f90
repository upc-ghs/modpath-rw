module ParticleTrackingEngineModule
  use TrackPathResultModule,only : TrackPathResultType
  use ParticleLocationModule,only : ParticleLocationType
  use ParticleLocationListModule,only : ParticleLocationListType
  use ParticleCoordinateModule,only : ParticleCoordinateType
  use TrackCellModule,only : TrackCellType
  use TrackCellResultModule,only : TrackCellResultType
  use ModflowRectangularGridModule,only : ModflowRectangularGridType
  use ModpathCellDataModule,only : ModpathCellDataType
  use ModpathSubCellDataModule,only : ModpathSubCellDataType
  use ParticleTrackingOptionsModule,only : ParticleTrackingOptionsType
  use FlowModelDataModule,only : FlowModelDataType
  implicit none
  
  ! Set default access status to private
  private

  type,public :: ParticleTrackingEngineType

    type(ParticleTrackingOptionsType) :: TrackingOptions
    type(ModpathCellDataType)         :: CellDataBuffer
    type(ParticleLocationListType)    :: LocBuffP
    type(ParticleLocationListType)    :: LocBuffTS
    logical :: Initialized = .false.

    ! Derived type pointers
    type(FlowModelDataType), pointer :: FlowModelData => null()
    class(ModflowRectangularGridType),pointer :: Grid => null()

    ! Private variables
    type(TrackCellType),private :: TrackCell
    type(TrackCellResultType),private :: TrackCellResult

    ! RWPT
    !type(ModpathSubCellDataType), dimension(18) :: NeighborSubCellData
    !type(ModpathCellDataType), dimension(18) :: NeighborCellData
    procedure(FillNeighborCells), pass, pointer :: FillNeighborCellData=>null()
    procedure(FillCellData)     , pass, pointer :: FillCellBuffer=>null()

  contains

    procedure :: Initialize=>pr_Initialize
    procedure :: Reset=>pr_Reset
    procedure :: TrackPath=>pr_TrackPath
    !procedure :: FillCellBuffer=>pr_FillCellBuffer
    procedure :: FillFaceFlowsBuffer=>pr_FillFaceFlowsBuffer
    procedure :: FillCellFlowsBuffer=>pr_FillCellFlowsBuffer
    procedure :: GetVolumetricBalanceSummary=>pr_GetVolumetricBalanceSummary
    procedure :: WriteCellBuffer=>pr_WriteCellBuffer
    procedure :: GetTopMostActiveCell=>pr_GetTopMostActiveCell

    ! RWPT
    !procedure :: FillNeighborCellData=>pr_FillNeighborCellData
    ! DEPRECATION WARNING
    !procedure :: FillNeighborSubCellData=>pr_FillNeighborSubCellData

  end type

  ! RWPT
  ! Interfaces
  abstract interface
  
      ! FillNeighborCells
      subroutine FillNeighborCells( this, neighborCellData )
          import ParticleTrackingEngineType
          import ModpathCellDataType
          !----------------------------------------
          class(ParticleTrackingEngineType) :: this
          type(ModpathCellDataType), dimension(2,18), intent(inout) :: neighborCellData
      end subroutine FillNeighborCells

      ! FillCellData 
      subroutine FillCellData( this, cellNumber, cellBuffer )
          import ParticleTrackingEngineType
          import ModpathCellDataType
          !----------------------------------------
          class(ParticleTrackingEngineType) :: this
          integer,intent(in) :: cellNumber
          type(ModpathCellDataType),intent(inout) :: cellBuffer
      end subroutine FillCellData

  end interface


contains


    function pr_GetTopMostActiveCell(this, cellNumber) result(n)
    !***************************************************************************************************************
    ! Description goes here
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(ParticleTrackingEngineType) :: this
    integer,intent(in) :: cellNumber
    integer :: n
    !---------------------------------------------------------------------------------------------------------------
   

        n = cellNumber
        do while(.true.)
            if(n .eq. 0) return
            if(this%FlowModelData%IBoundTS(n) .ne. 0) return
            n = this%Grid%GetFaceConnection(n, 5, 1)
        end do
    

    end function pr_GetTopMostActiveCell


    subroutine pr_WriteCellBuffer(this, unit, cellData, backwardTracking)
    !***************************************************************************************************************
    ! Description goes here
    !***************************************************************************************************************
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
    !---------------------------------------------------------------------------------------------------------------
      
      
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
                    if(m .eq. intervalCount) balanceCounts(intervalCount + 1) = &
                      balanceCounts(intervalCount + 1) + 1
                end do
            end if
        end do
   

    end subroutine pr_GetVolumetricBalanceSummary
    
    
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


    subroutine pr_Initialize(this, grid, trackingOptions, flowModelData)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(ParticleTrackingEngineType) :: this
    class(ModflowRectangularGridType),intent(inout),pointer :: grid
    type(ParticleTrackingOptionsType),intent(in) :: trackingOptions
    type(FlowModelDataType), intent(in), target :: flowModelData
    !---------------------------------------------------------------------------------------------------------------
   

        this%Initialized = .false.
        
        ! Call Reset to make sure that all arrays are initially unallocated
        call this%Reset()
        
        ! Initialize pointers and tracking options
        this%FlowModelData => flowModelData
        this%Grid => grid
        this%TrackingOptions = trackingOptions
        
        
        this%Initialized = .true.
       
        ! RWPT
        ! Assign interface for neighbor cell data generation
        select case( this%Grid%GridType ) 
            case (1)
                ! Set for structured grid
                this%FillCellBuffer=>pr_FillCellBufferStructured
                this%FillNeighborCellData=> pr_FillNeighborCellDataStructured
            case (2)
                ! Set for MODFLOW-USG unstructured grid
                this%FillCellBuffer=>pr_FillCellBufferUnstructured
                this%FillNeighborCellData=> pr_FillNeighborCellDataUnstructured
            case (3)
                ! Set for MODFLOW-6 structured grid (DIS)
                this%FillCellBuffer=>pr_FillCellBufferUnstructured
                this%FillNeighborCellData=> pr_FillNeighborCellDataStructured
            case (4)
                ! Set for MODFLOW-6 unstructured grid (DISV)
                this%FillCellBuffer=>pr_FillCellBufferUnstructured
                this%FillNeighborCellData=> pr_FillNeighborCellDataUnstructured
        end select

        
    end subroutine pr_Initialize


    subroutine  pr_FillFaceFlowsBuffer(this,buffer,bufferSize,count)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
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
        
        count = size(this%FlowModelData%FlowsJA)
        do n = 1, count
            buffer(n) = this%FlowModelData%FlowsJA(n)
        end do
    

    end subroutine pr_FillFaceFlowsBuffer
 

    subroutine pr_FillCellFlowsBuffer(this,cellNumber,buffer,bufferSize,count)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
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
            buffer(n) = this%FlowModelData%FlowsJA(offset + n)
        end do
      
    end subroutine pr_FillCellFlowsBuffer
 
  
    !subroutine pr_FillCellBuffer(this, cellNumber, cellBuffer)
    !!***************************************************************************************************************
    !!
    !!***************************************************************************************************************
    !! Specifications
    !!---------------------------------------------------------------------------------------------------------------
    !implicit none
    !class(ParticleTrackingEngineType) :: this
    !integer,intent(in) :: cellNumber
    !type(ModpathCellDataType),intent(inout) :: cellBuffer
    !doubleprecision,dimension(6) :: boundaryFlows
    !integer :: n, layer, boundaryFlowsOffset, gridType, cellType
    !!---------------------------------------------------------------------------------------------------------------
    !

    !    boundaryFlowsOffset = 6 * (cellNumber - 1)
    !    do n = 1, 6
    !        boundaryFlows(n) = this%FlowModelData%BoundaryFlows(boundaryFlowsOffset + n)
    !    end do
    !    
    !    layer = this%Grid%GetLayer(cellNumber)
    !    
    !    gridType = this%Grid%GridType
    !    cellType = this%Grid%CellType(cellNumber)
    !    select case (gridType)
    !        case (1)
    !            ! Set cell buffer data for a structured grid
    !            call cellBuffer%SetDataStructured(cellNumber,this%Grid%CellCount,     &
    !              this%Grid,this%FlowModelData%IBound,this%FlowModelData%IBoundTS,                        &
    !              this%FlowModelData%Porosity(cellNumber),this%FlowModelData%Retardation(cellNumber),     &
    !              this%FlowModelData%StorageFlows(cellNumber),this%FlowModelData%SourceFlows(cellNumber), &
    !              this%FlowModelData%SinkFlows(cellNumber), this%FlowModelData%FlowsRightFace,            &
    !              this%FlowModelData%FlowsFrontFace, this%FlowModelData%FlowsLowerFace, boundaryFlows,    &
    !              this%FlowModelData%Heads(cellNumber), cellType,                                   &
    !              this%FlowModelData%Zones(cellNumber))
    !        case (2)
    !            ! Set cell buffer data for a MODFLOW-USG unstructured grid
    !            call cellBuffer%SetDataUnstructured(cellNumber,this%Grid%CellCount,   &
    !              this%Grid%JaCount,this%Grid,                                        &
    !              this%FlowModelData%IBound,this%FlowModelData%IBoundTS,                                  &
    !              this%FlowModelData%Porosity(cellNumber),this%FlowModelData%Retardation(cellNumber),     &
    !              this%FlowModelData%StorageFlows(cellNumber),this%FlowModelData%SourceFlows(cellNumber), &
    !              this%FlowModelData%SinkFlows(cellNumber), this%FlowModelData%FlowsJA, boundaryFlows,    &
    !              this%FlowModelData%Heads(cellNumber), cellType,                                   &
    !              this%FlowModelData%Zones(cellNumber))
    !            ! Compute internal sub-cell face flows for cells with multiple sub-cell
    !            if(cellBuffer%GetSubCellCount() .gt. 1) then
    !                call cellBuffer%ComputeSubCellFlows()
    !            end if
    !        case (3)
    !            ! Set cell buffer data for a MODFLOW-6 structured grid (DIS)
    !            call cellBuffer%SetDataUnstructured(cellNumber,this%Grid%CellCount,   &
    !              this%Grid%JaCount,this%Grid,                                        &
    !              this%FlowModelData%IBound,this%FlowModelData%IBoundTS,                                  &
    !              this%FlowModelData%Porosity(cellNumber),this%FlowModelData%Retardation(cellNumber),     &
    !              this%FlowModelData%StorageFlows(cellNumber),this%FlowModelData%SourceFlows(cellNumber), &
    !              this%FlowModelData%SinkFlows(cellNumber), this%FlowModelData%FlowsJA, boundaryFlows,    &
    !              this%FlowModelData%Heads(cellNumber), cellType,                                   &
    !              this%FlowModelData%Zones(cellNumber))
    !        case (4)
    !            ! Set cell buffer data for a MODFLOW-6 unstructured grid (DISV)
    !            call cellBuffer%SetDataUnstructured(cellNumber,this%Grid%CellCount,   &
    !              this%Grid%JaCount,this%Grid,                                        &
    !              this%FlowModelData%IBound,this%FlowModelData%IBoundTS,                                  &
    !              this%FlowModelData%Porosity(cellNumber),this%FlowModelData%Retardation(cellNumber),     &
    !              this%FlowModelData%StorageFlows(cellNumber),this%FlowModelData%SourceFlows(cellNumber), &
    !              this%FlowModelData%SinkFlows(cellNumber), this%FlowModelData%FlowsJA, boundaryFlows,    &
    !              this%FlowModelData%Heads(cellNumber), cellType,                                   &
    !              this%FlowModelData%Zones(cellNumber))
    !             ! Compute internal sub-cell face flows for cells with multiple sub-cells
    !            if(cellBuffer%GetSubCellCount() .gt. 1) then
    !                call cellBuffer%ComputeSubCellFlows()
    !            end if
    !           
    !        case default
    !        ! Write error message and stop
    !    end select


    !end subroutine pr_FillCellBuffer


    subroutine pr_FillCellBufferUnstructured(this, cellNumber, cellBuffer)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
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
            boundaryFlows(n) = this%FlowModelData%BoundaryFlows(boundaryFlowsOffset + n)
        end do
        
        cellType = this%Grid%CellType(cellNumber)

        ! Set cell buffer data for a MODFLOW-USG unstructured grid
        ! Set cell buffer data for a MODFLOW-6 structured grid (DIS)
        ! Set cell buffer data for a MODFLOW-6 unstructured grid (DISV)
        call cellBuffer%SetDataUnstructured(cellNumber,this%Grid%CellCount,   &
          this%Grid%JaCount,this%Grid,                                        &
          this%FlowModelData%IBound,this%FlowModelData%IBoundTS,                                  &
          this%FlowModelData%Porosity(cellNumber),this%FlowModelData%Retardation(cellNumber),     &
          this%FlowModelData%StorageFlows(cellNumber),this%FlowModelData%SourceFlows(cellNumber), &
          this%FlowModelData%SinkFlows(cellNumber), this%FlowModelData%FlowsJA, boundaryFlows,    &
          this%FlowModelData%Heads(cellNumber), cellType,                                   &
          this%FlowModelData%Zones(cellNumber))


        return

    end subroutine pr_FillCellBufferUnstructured


    subroutine pr_FillCellBufferStructured(this, cellNumber, cellBuffer)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
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
            boundaryFlows(n) = this%FlowModelData%BoundaryFlows(boundaryFlowsOffset + n)
        end do
        
        cellType = this%Grid%CellType(cellNumber)

        ! Set cell buffer data for a structured grid
        call cellBuffer%SetDataStructured(cellNumber,this%Grid%CellCount,     &
          this%Grid,this%FlowModelData%IBound,this%FlowModelData%IBoundTS,                        &
          this%FlowModelData%Porosity(cellNumber),this%FlowModelData%Retardation(cellNumber),     &
          this%FlowModelData%StorageFlows(cellNumber),this%FlowModelData%SourceFlows(cellNumber), &
          this%FlowModelData%SinkFlows(cellNumber), this%FlowModelData%FlowsRightFace,            &
          this%FlowModelData%FlowsFrontFace, this%FlowModelData%FlowsLowerFace, boundaryFlows,    &
          this%FlowModelData%Heads(cellNumber), cellType,                                   &
          this%FlowModelData%Zones(cellNumber))

        return

    end subroutine pr_FillCellBufferStructured



    subroutine pr_Reset(this)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
    !
    ! Specifications
    !---------------------------------------------------------------------------------------------------------------
    implicit none
    class(ParticleTrackingEngineType) :: this
    !---------------------------------------------------------------------------------------------------------------
   
        this%Grid => null()
        this%FlowModelData => null()
        
    end subroutine pr_Reset


    subroutine pr_TrackPath(this, trackPathResult, traceModeOn, traceModeUnit,      &
      group, particleID, seqNumber, location, maximumTrackingTime, timeseriesPoints,&
      timeseriesPointCount)
    !***************************************************************************************************************
    !
    !***************************************************************************************************************
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

    ! RWPT
    type(ModpathCellDataType), dimension( 2, 18 ) :: neighborCellData

    ! OBS
    integer   :: idObservationCell
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
        this%TrackCell%SteadyState = this%FlowModelData%SteadyState
        this%TrackCell%TrackingOptions = this%TrackingOptions
        call this%FillCellBuffer(loc%CellNumber,  this%TrackCell%CellData)


        ! RWPT
        if (this%TrackingOptions%RandomWalkParticleTracking) then
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
                call this%TrackCell%ExecuteRandomWalkParticleTracking(loc, stopTime, this%TrackCellResult, &
                                                                                          neighborCellData )
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
               
                !! SOMEWHERE AROUND THIS BLOCK SHOULD SOLVE THE RWPT PARTICLE REBOUND

                ! If NextCellNumber is > 0, it means the particle has moved to another cell. 
                ! If so, convert loc from the current cell coordinates to the equivalent location in the new cell.
                nextCell = this%TrackCellResult%NextCellNumber
                if(nextCell .gt. 0) then
                    if(this%FlowModelData%IBoundTS(nextCell) .ne. 0) then
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
                        print *, 'INACTIVEEEE', nextCell
                        ! If next cell is inactive, it implies that a boundary face has been reached. 
                        ! Set status and return.
                        continueLoop = .false.
                        trackPathResult%Status = trackPathResult%Status_ReachedBoundaryFace()        
                    end if
                else
                    print *, 'CELLNUMBER ', nextCell 
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
               call WriteTraceData(traceModeUnit, this%TrackCell,                   &
                 this%TrackCellResult, this%FlowModelData%GetCurrentStressPeriod(), &
                 this%FlowModelData%GetCurrentTimeStep())
               !$omp end critical (tracedata)
            end if
        
            ! If there are observation cells
            if ( this%TrackingOptions%observationSimulation ) then
                ! Determine if the current TrackCell is on 
                ! array of observations cells and extract
                ! the corresponding unit

                ! Get id of observation cell in observation cell list
                ! Default is -999
                idObservationCell = this%TrackingOptions%IdObservationCell( this%TrackCell%CellData%CellNumber ) 
                if ( idObservationCell .ge. 0 ) then
                    call WriteObservationCellRecord( this, group, particleID,      &
                         this%TrackCell,                                           &
                         this%TrackingOptions%observationUnits( idObservationCell ))
                end if
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
            call this%Grid%ConvertToModelXYZ(pCoord%CellNumber, pCoord%LocalX, &
              pCoord%LocalY, pCoord%LocalZ, pCoord%GlobalX, pCoord%GlobalY,    &
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
            call this%Grid%ConvertToModelXYZ(pCoord%CellNumber, pCoord%LocalX, &
              pCoord%LocalY, pCoord%LocalZ, pCoord%GlobalX, pCoord%GlobalY,    &
              pCoord%GlobalZ)
            call trackPathResult%ParticlePath%Timeseries%AddItem(pCoord)
        end do
    

    end subroutine pr_TrackPath


! RWPT
subroutine pr_FillNeighborCellDataUnstructured( this, neighborCellData )
    !------------------------------------------------------------------
    !------------------------------------------------------------------
    ! Specifications
    !------------------------------------------------------------------
    implicit none
    class(ParticleTrackingEngineType) :: this
    type(ModpathCellDataType), dimension(2, 18), intent(inout) :: neighborCellData
    integer :: n, m, id, cellCounter, firstNeighborFaceNumber
    logical :: forceCellRefinement
    !------------------------------------------------------------------
    ! Reset all buffers
    do m = 1, 18 
        call neighborCellData( 1, m )%Reset() 
        call neighborCellData( 2, m )%Reset() 
    end do 

    ! Verify if current TrackCell is refined 
    if ( this%TrackCell%CellData%GetSubCellCount() .gt. 1 ) then 
        forceCellRefinement = .true.
    else
        forceCellRefinement = .false.
    end if 

    ! Loop through cell faces and fill neighborhood
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

        ! Thus this loop fill buffers corresponding 
        ! to: 
        !   - cellCounter     : direct connection 
        !   - cellCounter + 1 : connected neighbor 1, (firstNeighborFaceNumber)
        !   - cellCounter + 2 : connected neighbor 2, (firstNeighborFaceNumber+1)
        call pr_FillNeighborCellsSubBuffer( this, this%TrackCell%CellData, n, & 
            firstNeighborFaceNumber, neighborCellData, cellCounter, forceCellRefinement ) 

    end do


end subroutine pr_FillNeighborCellDataUnstructured


!RWPT
subroutine pr_FillNeighborCellsSubBuffer( this, &
        centerCellDataBuffer, faceNumber, firstNeighborFaceNumber, & 
        neighborCellsDataBuffer, directConnectionBufferIndex, forceCellRefinement )
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
    implicit none
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
        ! Consider skipping this if simulation is 2D
        
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


    ! Horizontal axis faceNumber

    ! Connected to one or more cells ?
    if ( centerCellDataBuffer%SubFaceCounts( faceNumber ) .gt. 1 ) then
        ! Connected to more than one cell

        ! Fill direct connection buffer 
        call pr_FillNeighborCellsConnectionFromHorizontalFace( this, centerCellDataBuffer, &
                    faceNumber, neighborCellsDataBuffer( :, directConnectionBufferIndex ), &
                                                                       forceCellRefinement )

        ! Fill indirect connection buffers

        ! Vertical face request
        if ( firstNeighborFaceNumber .ge. 5 ) then 
            ! Consider skipping this if simulation is 2D

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
            ! If requested neighbor cells are through 
            ! horizontal faces, change centerCellDataBuffer to the 
            ! recently filled buffer and request info. 
            ! Do not forceRefinement 

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
                                                                                          .false. )

            call pr_FillNeighborCellsConnectionFromHorizontalFace( this,                              &
                    neighborCellsDataBuffer( faceConnectionIndexes(2), directConnectionBufferIndex ), &
                                                                         firstNeighborFaceNumber + 1, &
                                       neighborCellsDataBuffer( :, directConnectionBufferIndex + 2 ), &
                                                                                              .false. )

            ! Done 
            return


        end if


    else if ( centerCellDataBuffer%GetFaceConnection( faceNumber, 1 ) .gt. 0 ) then
        ! Connected to one cell

        ! Fill direct connection buffer 
        call pr_FillNeighborCellsConnectionFromHorizontalFace( this, centerCellDataBuffer, &
                    faceNumber, neighborCellsDataBuffer( :, directConnectionBufferIndex ), &
                                                                       forceCellRefinement )

        ! If the buffer was filled with one of the subcells
        ! from a bigger cell
        if ( neighborCellsDataBuffer( 1, directConnectionBufferIndex )%fromSubCell ) then 

            ! Fill indirect connection buffers

            ! Vertical face request
            if ( firstNeighborFaceNumber .ge. 5 ) then 

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

            ! Horizontal face request
            else 

                ! Requires  a parent cellDataBuffer
                call pr_FillCellFromRealData( this, &
                    neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentCellNumber, &
                                                                      auxCellDataBuffer, .true. )
                auxCellDataBuffer%isParentCell    = .true. 
                auxCellDataBuffer%parentSubRow    = neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentSubRow
                auxCellDataBuffer%parentSubColumn = neighborCellsDataBuffer( 1, directConnectionBufferIndex )%parentSubColumn

                ! Verify direction of requested faces
                if ( firstNeighborFaceNumber .ge. 3 ) then
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
                    ! This case does not happens under the current scheme. 
                    ! That is because x axis faces are requested after 
                    ! a vertical faceNumber. As the grid cell is assumed
                    ! to have same refinement features for all layers, 
                    ! the vertical connection is never filled from a 
                    ! a bigger cell, only from cells same size as the 
                    ! current/center cell.
                    directionId = 1
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

        ! Not filled from subcell of bigger cell
        else

            ! Fill indirect connection buffers
            ! taking as center the already filled buffer 

            ! Vertical face request
            if ( firstNeighborFaceNumber .ge. 5 ) then 
                ! Consider skipping this if simulation is 2D
        
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
                call pr_FillCellFromRealData( this,    &
                    neighborCellsDataBuffer(           &
                        1, directConnectionBufferIndex &
                    )%GetFaceConnection( firstNeighborFaceNumber + 1, 1 ), &
                    neighborCellsDataBuffer( 1, directConnectionBufferIndex + 2 ), &
                                                               forceCellRefinement )
                neighborCellsDataBuffer( 1, directConnectionBufferIndex + 2 )%requestedFromDirection = 3

            ! Horizontal face request
            else

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
    implicit none
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


    ! Verify how many horizontal connections
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
                print *, 'ParticleTrackingEngine:FillNeighborCellsConnectionFromHorizontalFace: inconsistency'
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
                call pr_FillCellBufferFromSubCell( this, auxCellDataBuffer,    &
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
    implicit none 
    class(ParticleTrackingEngineType) :: this
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
    implicit none
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

    ! Done
    return


end subroutine pr_FillCellFromRealData



! RWPT
subroutine pr_FillNeighborCellDataStructured( this, neighborCellData )
    !------------------------------------------------------------------
    ! 
    !------------------------------------------------------------------
    ! Specifications
    !------------------------------------------------------------------
    implicit none
    class(ParticleTrackingEngineType) :: this
    type(ModpathCellDataType), dimension(2, 18), intent(inout) :: neighborCellData
    integer :: n, m, cellCounter, subCount
    !------------------------------------------------------------------


    cellCounter = 0
    do n = 1, 6
        cellCounter = cellCounter + 1
  
        ! If connection exists
        if ( this%TrackCell%CellData%GetFaceConnection(n,1) .gt. 0) then
  
            ! Fill the data buffer
            call this%FillCellBuffer( this%TrackCell%CellData%GetFaceConnection(n,1) , neighborCellData( 1, cellCounter ) )
  
            ! Extract surrounding cells from connection
            ! to complete information required for 
            ! computation of cell corner quantities 
            select case (n)
                case (1)
                    subCount = 0
                    do m = 3,4
                        subCount = subCount + 1 
                        if ( neighborCellData( 1, cellCounter )%GetFaceConnection(m,1) .gt. 0 ) then
                            call this%FillCellBuffer( neighborCellData( 1, cellCounter )%GetFaceConnection( m, 1 ), & 
                                                                      neighborCellData( 1, cellCounter + subCount ) )
                        else
                           call neighborCellData( 1, cellCounter + subCount )%Reset() 
                        end if
                    end do
                    cellCounter = cellCounter + subCount
                case (2)
                    subCount = 0
                    do m = 3,4
                        subCount = subCount + 1 
                        if ( neighborCellData( 1, cellCounter )%GetFaceConnection(m,1) .gt. 0 ) then
                            call this%FillCellBuffer( neighborCellData( 1, cellCounter )%GetFaceConnection( m, 1 ), & 
                                                                      neighborCellData( 1, cellCounter + subCount ) )
                        else
                           call neighborCellData( 1, cellCounter + subCount )%Reset() 
                        end if
                    end do
                    cellCounter = cellCounter + subCount
                case (3)
                    subCount = 0
                    do m = 5,6
                        subCount = subCount + 1 
                        if ( neighborCellData( 1, cellCounter )%GetFaceConnection(m,1) .gt. 0 ) then
                            call this%FillCellBuffer( neighborCellData( 1, cellCounter )%GetFaceConnection( m, 1 ), & 
                                                                      neighborCellData( 1, cellCounter + subCount ) )
                        else
                           call neighborCellData( 1, cellCounter + subCount )%Reset() 
                        end if
                    end do
                    cellCounter = cellCounter + subCount
                case (4)
                    subCount = 0
                    do m = 5,6
                        subCount = subCount + 1 
                        if ( neighborCellData( 1, cellCounter )%GetFaceConnection(m,1) .gt. 0 ) then
                            call this%FillCellBuffer( neighborCellData( 1, cellCounter )%GetFaceConnection( m, 1 ), & 
                                                                      neighborCellData( 1, cellCounter + subCount ) )
                        else
                           call neighborCellData( 1, cellCounter + subCount )%Reset() 
                        end if
                    end do
                    cellCounter = cellCounter + subCount
                case (5)
                    subCount = 0
                    do m = 1,2
                        subCount = subCount + 1 
                        if ( neighborCellData( 1, cellCounter )%GetFaceConnection(m,1) .gt. 0 ) then
                            call this%FillCellBuffer( neighborCellData( 1, cellCounter )%GetFaceConnection( m, 1 ), & 
                                                                      neighborCellData( 1, cellCounter + subCount ) )
                        else
                           call neighborCellData( 1, cellCounter + subCount )%Reset() 
                        end if
                    end do
                    cellCounter = cellCounter + subCount
                case (6)
                    subCount = 0
                    do m = 1,2
                        subCount = subCount + 1 
                        if ( neighborCellData( 1, cellCounter )%GetFaceConnection(m,1) .gt. 0 ) then
                            call this%FillCellBuffer( neighborCellData( 1, cellCounter )%GetFaceConnection( m, 1 ), & 
                                                                      neighborCellData( 1, cellCounter + subCount ) )
                        else
                           call neighborCellData( 1, cellCounter + subCount )%Reset() 
                        end if
                    end do
                    cellCounter = cellCounter + subCount
            end select
        else
            ! If no connection through cell face n, 
            ! reset buffers
            call neighborCellData( 1, cellCounter )%Reset()
            call neighborCellData( 1, cellCounter + 1 )%Reset()
            call neighborCellData( 1, cellCounter + 2 )%Reset()
            cellCounter = cellCounter + 2
        endif
    end do


end subroutine pr_FillNeighborCellDataStructured


! OBS
subroutine WriteObservationCellRecord( this, groupIndex, particleID, trackCell, outUnit)
    !--------
    ! Write observation cell record
    ! Doc me
    !----------
    ! Specifications
    !-----------------
    implicit none
    class(ParticleTrackingEngineType) :: this
    ! input
    type(TrackCellType), intent(in) :: trackCell
    integer, intent(in)             :: outUnit, groupIndex, particleID
    ! local
    doubleprecision     :: initialGlobalX, initialGlobalY, initialGlobalZ
    doubleprecision     :: finalGlobalX, finalGlobalY, finalGlobalZ
    doubleprecision     :: initialTime, finalTime
    !----------------------------------------------

    initialTime = trackCell%TrackSubCell%TrackSubCellResult%InitialLocation%TrackingTime
    call this%Grid%ConvertToModelXYZ( trackCell%CellData%CellNumber,      &
        trackCell%TrackSubCell%TrackSubCellResult%InitialLocation%LocalX, & 
        trackCell%TrackSubCell%TrackSubCellResult%InitialLocation%LocalY, &
        trackCell%TrackSubCell%TrackSubCellResult%InitialLocation%LocalZ, &
        initialGlobalX, initialGlobalY, initialGlobalZ )

    finalTime = trackCell%TrackSubCell%TrackSubCellResult%FinalLocation%TrackingTime
    call this%Grid%ConvertToModelXYZ( trackCell%CellData%CellNumber,    &
        trackCell%TrackSubCell%TrackSubCellResult%FinalLocation%LocalX, & 
        trackCell%TrackSubCell%TrackSubCellResult%FinalLocation%LocalY, &
        trackCell%TrackSubCell%TrackSubCellResult%FinalLocation%LocalZ, &
        finalGlobalX, finalGlobalY, finalGlobalZ )

    write(outUnit, '(2I8,es18.9e3,8es18.9e3)')                      &
      groupIndex, particleID,                                       & 
      initialTime, initialGlobalX, initialGlobalY, initialGlobalZ,  &
      finalTime, finalGlobalX, finalGlobalY, finalGlobalZ                        

end subroutine WriteObservationCellRecord




end module ParticleTrackingEngineModule
