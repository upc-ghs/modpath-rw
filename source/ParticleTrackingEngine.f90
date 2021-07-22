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

  
  contains

    procedure :: Initialize=>pr_Initialize
    procedure :: Reset=>pr_Reset
    procedure :: TrackPath=>pr_TrackPath
    procedure :: FillCellBuffer=>pr_FillCellBuffer
    procedure :: FillFaceFlowsBuffer=>pr_FillFaceFlowsBuffer
    procedure :: FillCellFlowsBuffer=>pr_FillCellFlowsBuffer
    procedure :: GetVolumetricBalanceSummary=>pr_GetVolumetricBalanceSummary
    procedure :: WriteCellBuffer=>pr_WriteCellBuffer
    procedure :: GetTopMostActiveCell=>pr_GetTopMostActiveCell

  end type
 

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
 
  
    subroutine pr_FillCellBuffer(this, cellNumber, cellBuffer)
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
        
        layer = this%Grid%GetLayer(cellNumber)
        
        gridType = this%Grid%GridType
        cellType = this%Grid%CellType(cellNumber)
        select case (gridType)
            case (1)
                ! Set cell buffer data for a structured grid
                call cellBuffer%SetDataStructured(cellNumber,this%Grid%CellCount,     &
                  this%Grid,this%FlowModelData%IBound,this%FlowModelData%IBoundTS,                        &
                  this%FlowModelData%Porosity(cellNumber),this%FlowModelData%Retardation(cellNumber),     &
                  this%FlowModelData%StorageFlows(cellNumber),this%FlowModelData%SourceFlows(cellNumber), &
                  this%FlowModelData%SinkFlows(cellNumber), this%FlowModelData%FlowsRightFace,            &
                  this%FlowModelData%FlowsFrontFace, this%FlowModelData%FlowsLowerFace, boundaryFlows,    &
                  this%FlowModelData%Heads(cellNumber), cellType,                                   &
                  this%FlowModelData%Zones(cellNumber))
            case (2)
                ! Set cell buffer data for a MODFLOW-USG unstructured grid
                call cellBuffer%SetDataUnstructured(cellNumber,this%Grid%CellCount,   &
                  this%Grid%JaCount,this%Grid,                                        &
                  this%FlowModelData%IBound,this%FlowModelData%IBoundTS,                                  &
                  this%FlowModelData%Porosity(cellNumber),this%FlowModelData%Retardation(cellNumber),     &
                  this%FlowModelData%StorageFlows(cellNumber),this%FlowModelData%SourceFlows(cellNumber), &
                  this%FlowModelData%SinkFlows(cellNumber), this%FlowModelData%FlowsJA, boundaryFlows,    &
                  this%FlowModelData%Heads(cellNumber), cellType,                                   &
                  this%FlowModelData%Zones(cellNumber))
                ! Compute internal sub-cell face flows for cells with multiple sub-cell
                if(cellBuffer%GetSubCellCount() .gt. 1) then
                    call cellBuffer%ComputeSubCellFlows()
                end if
            case (3)
                ! Set cell buffer data for a MODFLOW-6 structured grid (DIS)
                call cellBuffer%SetDataUnstructured(cellNumber,this%Grid%CellCount,   &
                  this%Grid%JaCount,this%Grid,                                        &
                  this%FlowModelData%IBound,this%FlowModelData%IBoundTS,                                  &
                  this%FlowModelData%Porosity(cellNumber),this%FlowModelData%Retardation(cellNumber),     &
                  this%FlowModelData%StorageFlows(cellNumber),this%FlowModelData%SourceFlows(cellNumber), &
                  this%FlowModelData%SinkFlows(cellNumber), this%FlowModelData%FlowsJA, boundaryFlows,    &
                  this%FlowModelData%Heads(cellNumber), cellType,                                   &
                  this%FlowModelData%Zones(cellNumber))
            case (4)
                ! Set cell buffer data for a MODFLOW-6 unstructured grid (DISV)
                call cellBuffer%SetDataUnstructured(cellNumber,this%Grid%CellCount,   &
                  this%Grid%JaCount,this%Grid,                                        &
                  this%FlowModelData%IBound,this%FlowModelData%IBoundTS,                                  &
                  this%FlowModelData%Porosity(cellNumber),this%FlowModelData%Retardation(cellNumber),     &
                  this%FlowModelData%StorageFlows(cellNumber),this%FlowModelData%SourceFlows(cellNumber), &
                  this%FlowModelData%SinkFlows(cellNumber), this%FlowModelData%FlowsJA, boundaryFlows,    &
                  this%FlowModelData%Heads(cellNumber), cellType,                                   &
                  this%FlowModelData%Zones(cellNumber))
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
        
        continueLoop = .true.
        isTimeSeriesPoint = .false.
        isMaximumTime = .false.
        
        do while(continueLoop)
            ! Check to see if the particle has moved to another cell. If so, load the new cell data
            if(loc%CellNumber .ne. this%TrackCell%CellData%CellNumber) then
                call this%FillCellBuffer(loc%CellNumber, this%TrackCell%CellData)
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
            
            ! Start with the particle loacion loc and track it through the cell until it reaches
            ! an exit face or the tracking time reaches the value specified by stopTime
            call this%TrackCell%ExecuteTracking(loc, stopTime, this%TrackCellResult)
            
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
               call WriteTraceData(traceModeUnit, this%TrackCell,                   &
                 this%TrackCellResult, this%FlowModelData%GetCurrentStressPeriod(), &
                 this%FlowModelData%GetCurrentTimeStep())
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
 

end module ParticleTrackingEngineModule
