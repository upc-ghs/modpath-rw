module ModpathCellDataModule
  use ModflowRectangularGridModule,only : ModflowRectangularGridType
  use ParticleLocationModule,only : ParticleLocationType
  use ModpathSubCellDataModule,only : ModpathSubCellDataType
  implicit none
  
! Set default access status to private
  private

! Private data type declarations

  
! Public derived data type definitions
!--------------------------------------
! type: 
!--------------------------------------
  type,public :: ModpathCellDataType
    ! public data
    integer :: CellNumber, Layer, Ibound, IboundTS, Zone, LayerType
    doubleprecision :: DX, DY, MinX, MinY, Bottom, Top, Head, Porosity,Retardation,SourceFlow,SinkFlow,StorageFlow
    ! RWPT: more info for mass transport models could be considered 
    doubleprecision :: ICBound, ICBoundTS 

    ! private data
    ! RWPT: remove privates 
    integer :: SubCellRowCount,SubCellColumnCount,ReducedConnectionCount
    integer,dimension(6) :: SubFaceCounts,PotentialConnectionsCount,SubFaceBoundaryCounts
    integer,dimension(2) :: SubFaceConn1,SubFaceConn2,SubFaceConn3,SubFaceConn4
    integer,dimension(4) :: SubFaceConn5,SubFaceConn6
    integer,dimension(2) :: MassBoundarySubFace1,MassBoundarySubFace2,MassBoundarySubFace3,MassBoundarySubFace4
    integer,dimension(4) :: MassBoundarySubFace5,MassBoundarySubFace6
    doubleprecision,dimension(4) :: SubCellFlows,Q5,Q6
    doubleprecision,dimension(2) :: Q1,Q2,Q3,Q4
    integer :: ArraySizeMode = 1

    ! RWPT-USG NEW PROPERTIES
    logical :: fromSubCell    = .false.
    logical :: fromSubSubCell = .false.
    logical :: isParentCell   = .false.
    integer :: parentCellNumber
    integer :: parentSubRow, parentSubColumn
    integer :: requestedFromDirection

    ! RWPT transport parameters
    doubleprecision, public :: alphaL, alphaT ! to be deprecated
    doubleprecision, public :: alphaLH, alphaLV, alphaTH, alphaTV
    doubleprecision, public :: dMEff

    ! RWPT convertible cells parameters
    logical :: dry
    logical :: partiallyDry

    ! RWPT Initialized flag
    logical :: initialized

  contains
    procedure :: GetDZ=>pr_GetDZ
    procedure :: GetDZRW=>pr_GetDZRW
    procedure :: GetArraySizeMode=>pr_GetArraySizeMode
    procedure :: SetArraySizeMode=>pr_SetArraySizeMode
    procedure :: SetFlowAndPropertyData=>pr_SetFlowAndPropertyData
    procedure :: Reset=>pr_Reset
    procedure :: ResetArrays=>pr_ResetArrays
    procedure :: HasExitFace=>pr_HasExitFace
    procedure :: GetSubCellRowCount=>pr_GetSubCellRowCount
    procedure :: GetSubCellColumnCount=>pr_GetSubCellColumnCount
    procedure :: GetSubCellCount=>pr_GetSubCellCount
    procedure :: GetReducedConnectionCount=>pr_GetReducedConnectionCount
    procedure :: GetSubFaceCount=>pr_GetSubFaceCount
    procedure :: GetFaceFlow=>pr_GetFaceFlow
    procedure :: GetFaceConnection=>pr_GetFaceConnection
    procedure :: SetSubCellFlows=>pr_SetSubCellFlows
    procedure :: GetSubCellFlow=>pr_GetSubCellFlow
    procedure :: GetAveragedFaceFlow=>pr_GetAveragedFaceFlow
    procedure :: FillMassSubCellDataBuffer=>pr_FillMassSubCellDataBuffer
    procedure :: FillSubCellDataBuffer=>pr_FillSubCellDataBuffer
    procedure :: GetSubCellData=>pr_GetSubCellData
    procedure :: FillSubCellFaceFlowsBuffer=>pr_FillSubCellFaceFlowsBuffer
    procedure :: AssignAveragedFaceFlowArray=>pr_AssignAveragedFaceFlowArray
    procedure :: GetVolumetricBalance=>pr_GetVolumetricBalance
    procedure :: GetVolumetricBalanceComponents=>pr_GetVolumetricBalanceComponents
    procedure :: ComputeSubCellFlows=>pr_ComputeSubCellFlows
    !procedure,private :: GetSubFaceBoundaryCount=>pr_GetSubFaceBoundaryCount
    procedure,private :: GetSubCellBoundaryFlow=>pr_GetSubCellBoundaryFlow
    procedure,private :: pr_SetDataUnstructured
    procedure,private :: pr_SetDataStructured
    procedure :: SetDataUnstructured=>pr_SetDataUnstructured
    procedure :: SetMassTransportDataUnstructured=>pr_SetMassTransportDataUnstructured
    procedure :: SetDataStructured=>pr_SetDataStructured
    procedure :: SetMassTransportDataStructured=>pr_SetMassTransportDataStructured
    generic :: SetData=>pr_SetDataUnstructured

    ! RWPT-USG
    procedure :: GetNeighborSubCellIndexes => pr_GetNeighborSubCellIndexes
    procedure :: FillSubCellFaceAreas      => pr_FillSubCellFaceAreas
    procedure :: FillSubCellFaceFlows      => pr_FillSubCellFaceFlows
    procedure :: GetVolume                 => pr_GetVolume
    procedure :: GetWaterVolume            => pr_GetWaterVolume

    ! RWPT: could be replaced by getdzrw
    procedure :: VerifyDryCell             => pr_VerifyDryCell

  end type


contains


  function pr_GetDZ(this) result(dz)
  implicit none
  class(ModpathCellDataType) :: this
  doubleprecision :: dz
  
  dz = this%Top - this%Bottom

  ! If the layer is convertible, set dz = Head - Bottom if Head < Top
  if(this%LayerType .eq. 1) then
      if(this%Head .lt. this%Top) dz = this%Head - this%Bottom
      ! If dz < 0, set dz to an arbitrary, small positive value
      if(dz .lt. 0.0d0) dz = 1.0d-4
  end if
  
  end function pr_GetDZ




  function pr_GetVolumetricBalance(this) result(balance)
  implicit none
  class(ModpathCellDataType) :: this
  doubleprecision :: balance,inflow,outflow
  integer :: n
  
  balance = 0.0d0
  inflow = 0.0d0
  outflow = 0.0d0
  
  inflow = this%SourceFlow
  outflow = -this%SinkFlow
  if(this%StorageFlow .ge. 0.0d0) then
      inflow = inflow + this%StorageFlow
  else
      outflow = outflow - this%StorageFlow
  end if
  
  ! Face 1
  do n = 1, this%SubFaceCounts(1)
      if(this%Q1(n) .ge. 0.0) then
          inflow = inflow + this%Q1(n)
      else
          outflow = outflow - this%Q1(n)
      end if
  end do
  
  ! Face 2
  do n = 1, this%SubFaceCounts(2)
      if(this%Q2(n) .le. 0.0) then
          inflow = inflow - this%Q2(n)
      else
          outflow = outflow + this%Q2(n)
      end if
  end do
  
  ! Face 3
  do n = 1, this%SubFaceCounts(3)
      if(this%Q3(n) .ge. 0.0) then
          inflow = inflow + this%Q3(n)
      else
          outflow = outflow - this%Q3(n)
      end if
  end do
  
  ! Face 4
  do n = 1, this%SubFaceCounts(4)
      if(this%Q4(n) .le. 0.0) then
          inflow = inflow - this%Q4(n)
      else
          outflow = outflow + this%Q4(n)
      end if
  end do
  
  ! Face 5
  do n = 1, this%SubFaceCounts(5)
      if(this%Q5(n) .ge. 0.0) then
          inflow = inflow + this%Q5(n)
      else
          outflow = outflow - this%Q5(n)
      end if
  end do
  
  ! Face 6
  do n = 1, this%SubFaceCounts(6)
      if(this%Q6(n) .le. 0.0) then
          inflow = inflow - this%Q6(n)
      else
          outflow = outflow + this%Q6(n)
      end if
  end do
  
  if((inflow .eq. 0.0d0) .and. (outflow .eq. 0.0d0)) then
      balance = 0.0d0
  else
      balance = 100.0d0 * (inflow - outflow) / ((inflow + outflow) / 2.0d0)
  end if
  
  end function pr_GetVolumetricBalance
  

  subroutine pr_GetVolumetricBalanceComponents(this, totalFaceInflow,           &
    totalFaceOutflow, sourceFlow, sinkFlow, storageInflow, storageOutflow,      &
    balance)
  implicit none
  class(ModpathCellDataType) :: this
  doubleprecision :: inflow,outflow
  doubleprecision,intent(inout) :: totalFaceInflow, totalFaceOutflow,           &
    sourceFlow, sinkFlow, storageInflow, storageOutflow, balance
  integer :: n
  
  balance = 0.0d0
  totalFaceInflow = 0.0d0
  totalFaceOutflow = 0.0d0
  sourceFlow = this%SourceFlow
  sinkFlow = -this%SinkFlow
  storageInflow = 0.0d0
  storageOutflow = 0.0d0
  if(this%StorageFlow .gt. 0.0d0) then
      storageInflow = this%StorageFlow
  else if(this%StorageFlow .lt. 0.0d0) then
      storageOutflow = -this%StorageFlow
  end if
  
  ! Face 1
  do n = 1, this%SubFaceCounts(1)
      if(this%Q1(n) .ge. 0.0) then
          totalFaceInflow = totalFaceInflow + this%Q1(n)
      else
          totalFaceOutflow = totalFaceOutflow - this%Q1(n)
      end if
  end do
  
  ! Face 2
  do n = 1, this%SubFaceCounts(2)
      if(this%Q2(n) .le. 0.0) then
          totalFaceInflow = totalFaceInflow - this%Q2(n)
      else
          totalFaceOutflow = totalFaceOutflow + this%Q2(n)
      end if
  end do
  
  ! Face 3
  do n = 1, this%SubFaceCounts(3)
      if(this%Q3(n) .ge. 0.0) then
          totalFaceInflow = totalFaceInflow + this%Q3(n)
      else
          totalFaceOutflow = totalFaceOutflow - this%Q3(n)
      end if
  end do
  
  ! Face 4
  do n = 1, this%SubFaceCounts(4)
      if(this%Q4(n) .le. 0.0) then
          totalFaceInflow = totalFaceInflow - this%Q4(n)
      else
          totalFaceOutflow = totalFaceOutflow + this%Q4(n)
      end if
  end do
  
  ! Face 5
  do n = 1, this%SubFaceCounts(5)
      if(this%Q5(n) .ge. 0.0) then
          totalFaceInflow = totalFaceInflow + this%Q5(n)
      else
          totalFaceOutflow = totalFaceOutflow - this%Q5(n)
      end if
  end do
  
  ! Face 6
  do n = 1, this%SubFaceCounts(6)
      if(this%Q6(n) .le. 0.0) then
          totalFaceInflow = totalFaceInflow - this%Q6(n)
      else
          totalFaceOutflow = totalFaceOutflow + this%Q6(n)
      end if
  end do
  
  ! Add in the internal sources and sinks to the appropriate inflow or outflow component
  inflow = totalFaceInflow + sourceFlow + storageInflow
  outflow = totalFaceOutflow + sinkFlow + storageOutflow
  
  if((inflow .eq. 0.0d0) .and. (outflow .eq. 0.0d0)) then
      balance = 0.0d0
  else
      balance = 100.0d0 * (inflow - outflow) / ((inflow + outflow) / 2.0d0)
  end if
  
  end subroutine pr_GetVolumetricBalanceComponents
  
!------------------------------------------
  function pr_GetArraySizeMode(this) result(mode)
  implicit none
  class(ModpathCellDataType) :: this
  integer :: mode
  
  mode = this%ArraySizeMode
  
  end function pr_GetArraySizeMode

!------------------------------------------
  subroutine pr_SetArraySizeMode(this,mode)
  implicit none
  class(ModpathCellDataType) :: this
  integer,intent(in) :: mode
  
  this%ArraySizeMode = mode
  call this%Reset()
  
  end subroutine pr_SetArraySizeMode
  
!------------------------------------------
  function pr_HasExitFace(this,backwardTracking) result(hasExit)
  implicit none
  class(ModpathCellDataType) :: this
  doubleprecision :: sign
  logical,intent(in) :: backwardTracking
  logical :: hasExit
  integer :: n
  doubleprecision :: q
  
  sign = 1.0d0
  if(backwardTracking) sign = - sign
  
  hasExit = .false.
  
  ! Check face 1
  do n = 1, this%subFaceCounts(1)
    q = sign*this%Q1(n)
    if(q .lt. 0.0d0) then
      hasExit = .true.
      return
    end if
  end do
  
  ! Check face 2
  do n = 1, this%subFaceCounts(2)
    q = sign*this%Q2(n)
    if(q .gt. 0.0d0) then
      hasExit = .true.
      return
    end if
  end do
  
  ! Check face 3
  do n = 1, this%subFaceCounts(3)
    q = sign*this%Q3(n)
    if(q .lt. 0.0d0) then
      hasExit = .true.
      return
    end if
  end do
  
  ! Check face 4
  do n = 1, this%subFaceCounts(4)
    q = sign*this%Q4(n)
    if(q .gt. 0.0d0) then
      hasExit = .true.
      return
    end if
  end do
  
  ! Check face 5
  do n = 1, this%subFaceCounts(5)
    q = sign*this%Q5(n)
    if(q .lt. 0.0d0) then
      hasExit = .true.
      return
    end if
  end do
  
  ! Check face 6
  do n = 1, this%subFaceCounts(6)
    q = sign*this%Q6(n)
    if(q .gt. 0.0d0) then
      hasExit = .true.
      return
    end if
  end do
  
  end function pr_HasExitFace
  
!------------------------------------------
  subroutine pr_Reset(this)
  implicit none
  class(ModpathCellDataType) :: this
  
  this%CellNumber = 0
  this%DX = 0
  this%DY = 0
  this%MinX = 0.0d0
  this%MinY = 0.0d0
  this%Bottom = 0.0d0
  this%Top = 0.0d0
  this%Ibound = 0
  this%Zone = 0
  this%Porosity = 0d0
  this%Retardation = 0d0
  this%SourceFlow = 0d0
  this%SinkFlow =  0d0
  this%StorageFlow = 0d0
  this%ReducedConnectionCount = 0
  this%SubCellRowCount = 1
  this%SubCellColumnCount = 1
  
  call this%ResetArrays()

  ! RWPT-USG 
  this%fromSubCell     = .false.
  this%fromSubSubCell  = .false.
  this%isParentCell    = .false.
  this%parentCellNumber = 0
  this%parentSubRow     = 0
  this%parentSubColumn  = 0
  this%requestedFromDirection = 0

  this%dry          = .false.
  this%partiallyDry = .false.

  this%initialized = .false.

  end subroutine pr_Reset


  !------------------------------------------
  subroutine pr_ResetArrays(this)
  implicit none
  class(ModpathCellDataType) :: this
  !------------------------------------------
  
  this%SubFaceCounts = 1
  this%SubFaceBoundaryCounts = 1
  this%Q1 = 0.0d0
  this%Q2 = 0.0d0
  this%Q3 = 0.0d0
  this%Q4 = 0.0d0
  this%Q5 = 0.0d0
  this%Q6 = 0.0d0
  this%SubFaceConn1 = 0
  this%SubFaceConn2 = 0
  this%SubFaceConn3 = 0
  this%SubFaceConn4 = 0
  this%SubFaceConn5 = 0
  this%SubFaceConn6 = 0
  this%MassBoundarySubFace1 = 0 
  this%MassBoundarySubFace2 = 0 
  this%MassBoundarySubFace3 = 0 
  this%MassBoundarySubFace4 = 0 
  this%MassBoundarySubFace5 = 0 
  this%MassBoundarySubFace6 = 0 
  this%SubCellFlows = 0.0d0


  end subroutine pr_ResetArrays


!------------------------------------------
  subroutine pr_SetFlowAndPropertyData(this,ibound,porosity,retardation,        &
    arraySize,faceFlowsCount,faceFlows,connectionList,storageFlow,sourceFlow,   &
    sinkFlow,boundaryFlows)
  implicit none
  class(ModpathCellDataType) :: this
  integer :: n,index,cellNumber,count
  integer,intent(in) :: ibound,faceFlowsCount,arraySize
  integer,intent(in),dimension(arraySize) :: connectionList
  doubleprecision :: flow
  doubleprecision,intent(in) :: porosity,retardation,storageFlow,sourceFlow,sinkFlow
  doubleprecision,intent(in),dimension(6) :: boundaryFlows
  doubleprecision,intent(in),dimension(arraySize) :: faceFlows
  
  ! Assign property data
  this%Ibound = ibound
  this%Porosity = porosity
  this%Retardation = retardation
  this%SourceFlow = sourceFlow
  this%SinkFlow = sinkFlow
  this%StorageFlow = storageFlow

  ! Process face data
  
  ! Face 1
  count = this%SubFaceCounts(1)
  do n = 1, count
    cellNumber = this%SubFaceConn1(n)
    if(cellNumber .gt. 0) then
      index = pr_FindConnectionNumberIndex(cellNumber,connectionList,arraySize,faceFlowsCount)
      this%Q1(n) = faceFlows(index)
    else
      this%Q1(n) = 0.0d0
    end if
  end do
  
  ! Face 2
  count = this%SubFaceCounts(2)
  do n = 1, count
    cellNumber = this%SubFaceConn2(n)
    if(cellNumber .gt. 0) then
      index = pr_FindConnectionNumberIndex(cellNumber,connectionList,arraySize,faceFlowsCount)
      this%Q2(n) = faceFlows(index)
    else
      this%Q2(n) = 0.0d0
    end if
  end do
  
  ! Face 3
  count = this%SubFaceCounts(3)
  do n = 1, count
    cellNumber = this%SubFaceConn3(n)
    if(cellNumber .gt. 0) then
      index = pr_FindConnectionNumberIndex(cellNumber,connectionList,arraySize,faceFlowsCount)
      this%Q3(n) = faceFlows(index)
    else
      this%Q3(n) = 0.0d0
    end if
  end do
  
  ! Face 4
  count = this%SubFaceCounts(4)
  do n = 1, count
    cellNumber = this%SubFaceConn4(n)
    if(cellNumber .gt. 0) then
      index = pr_FindConnectionNumberIndex(cellNumber,connectionList,arraySize,faceFlowsCount)
      this%Q4(n) = faceFlows(index)
    else
      this%Q4(n) = 0.0d0
    end if
  end do
  
  ! Face 5
  count = this%SubFaceCounts(5)
  do n = 1, count
    cellNumber = this%SubFaceConn5(n)
    if(cellNumber .gt. 0) then
      index = pr_FindConnectionNumberIndex(cellNumber,connectionList,arraySize,faceFlowsCount)
      this%Q5(n) = faceFlows(index)
    else
      this%Q5(n) = 0.0d0
    end if
  end do
  
  ! Face 6
  count = this%SubFaceCounts(6)
  do n = 1, count
    cellNumber = this%SubFaceConn1(n)
    if(cellNumber .gt. 0) then
      index = pr_FindConnectionNumberIndex(cellNumber,connectionList,arraySize,faceFlowsCount)
      this%Q6(n) = faceFlows(index)
    else
      this%Q6(n) = 0.0d0
    end if
  end do
  
  ! Process boundary flow data
  
  ! Face 1
  if(boundaryFlows(1) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(1) .gt. 0) then
        flow = boundaryFlows(1)/this%SubFaceBoundaryCounts(1)
        do n = 1, this%SubFaceCounts(1)
          if(this%SubFaceConn1(n) .eq. 0) then
            this%Q1(n) = flow
          end if
        end do
      else
        if(boundaryFlows(1) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(1)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(1)
        end if
      end if
  end if

  ! Face 2
  if(boundaryFlows(2) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(2) .gt. 0) then
        flow = boundaryFlows(2)/this%SubFaceBoundaryCounts(2)
        do n = 1, this%SubFaceCounts(2)
          if(this%SubFaceConn2(n) .eq. 0) then
            this%Q2(n) = -flow
          end if
        end do
      else
        if(boundaryFlows(2) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(2)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(2)
        end if
      end if
  end if
  
  ! Face 3
  if(boundaryFlows(3) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(3) .gt. 0) then
        flow = boundaryFlows(3)/this%SubFaceBoundaryCounts(3)
        do n = 1, this%SubFaceCounts(3)
          if(this%SubFaceConn3(n) .eq. 0) then
            this%Q3(n) = flow
          end if
        end do
      else
        if(boundaryFlows(3) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(3)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(3)
        end if
      end if
  end if

  ! Face 4
  if(boundaryFlows(4) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(4) .gt. 0) then
        flow = boundaryFlows(4)/this%SubFaceBoundaryCounts(4)
        do n = 1, this%SubFaceCounts(4)
          if(this%SubFaceConn4(n) .eq. 0) then
            this%Q4(n) = -flow
          end if
        end do
      else
        if(boundaryFlows(4) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(4)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(4)
        end if
      end if
  end if
  
  ! Face 5
  if(boundaryFlows(5) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(5) .gt. 0) then
        flow = boundaryFlows(5)/this%SubFaceBoundaryCounts(5)
        do n = 1, this%SubFaceCounts(5)
          if(this%SubFaceConn5(n) .eq. 0) then
            this%Q5(n) = flow
          end if
        end do
      else
        if(boundaryFlows(5) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(5)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(5)
        end if
      end if
  end if

  ! Face 6
  if(boundaryFlows(6) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(6) .gt. 0) then
        flow = boundaryFlows(6)/this%SubFaceBoundaryCounts(6)
        do n = 1, this%SubFaceCounts(6)
          if(this%SubFaceConn6(n) .eq. 0) then
            this%Q6(n) = -flow
          end if
        end do
      else
        if(boundaryFlows(6) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(6)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(6)
        end if
      end if
  end if
  
  end subroutine pr_SetFlowAndPropertyData

  ! RWPT
  !---------------------------------------------------------------------------------------
  subroutine pr_SetMassTransportDataUnstructured(this,cellNumber,cellCount, &
          reducedConnectionCount,grid, ibound,iboundTS,porosity,retardation,& 
               storageFlow,sourceFlow,sinkFlow,faceFlows,boundaryFlows,head,& 
                                 layerType, zone, icboundTS, defaultICBound )
    !---------------------------------------------------------------------------------------
    implicit none
    class(ModpathCellDataType) :: this
    class(ModflowRectangularGridType),intent(in) :: grid
    integer,intent(in) :: cellNumber,cellCount
    integer,intent(in) :: reducedConnectionCount
    integer,intent(in) :: layerType, zone
    integer,intent(in),dimension(cellCount) :: ibound, iboundTS
    integer,intent(in),dimension(cellCount) :: icboundTS
    integer,intent(in)                      :: defaultICBound
    doubleprecision,intent(in) :: porosity, retardation, storageFlow, sourceFlow, sinkFlow
    doubleprecision,intent(in),dimension(6) :: boundaryFlows
    doubleprecision,intent(in),dimension(reducedConnectionCount) :: faceFlows
    doubleprecision,intent(in) :: head
    integer :: n,index,count,i,conn
    doubleprecision :: flow
    !---------------------------------------------------------------------------------------
  
    ! Is this reset necessary ?
    ! All scalar variables are reassigned in the following lines after reset
    !call this%Reset()
    call this%ResetArrays()

    this%CellNumber = cellNumber
    this%Layer = grid%GetLayer(cellNumber)
    this%DX = grid%DelX(cellNumber)
    this%DY = grid%DelY(cellNumber)
    this%MinX = grid%GetLeft(cellNumber)
    this%MinY = grid%GetFront(cellNumber)
    this%Bottom = grid%Bottom(cellNumber)
    this%Top = grid%Top(cellNumber)
    this%ReducedConnectionCount = grid%GetJaCellConnectionsCount(cellNumber)
    !this%ReducedConnectionCount = grid%GetReducedCellConnectionCount(cellNumber)
    
    ! Assign property data
    this%Zone = zone
    this%LayerType = layerType
    this%Head = head
    this%Ibound = ibound(cellNumber)
    this%IboundTS = iboundTS(cellNumber)
    this%Porosity = porosity
    this%Retardation = retardation
    this%SourceFlow = sourceFlow
    this%SinkFlow = sinkFlow
    this%StorageFlow = storageFlow
 

    ! RWPT 
    this%ICBoundTS = icboundTS(cellNumber)


    ! Process face flow data
    ! Face 1
    count = grid%GetPotentialFaceConnectionCount(cellNumber, 1)
    this%PotentialConnectionsCount(1) = count
    if(count .gt. 0) then
        this%SubFaceCounts(1) = count
        ! i: is a counter of boundary faces
        i = 0
        do n = 1, count
          conn = grid%GetFaceConnection(cellNumber,1,n)
          this%SubFaceConn1(n) = conn

          ! RWPT: it also includes processing of mass boundaries
          ! To avoid repeating if block
          if(conn .eq. 0) then
              i = i + 1
              ! Set default value when no connection 
              this%MassBoundarySubFace1(n) = defaultICBound
          else
              ! Flow boundary, increase i
              if(iboundTS(conn) .eq. 0) i = i + 1
              ! Mass boundary
              this%MassBoundarySubFace1(n) = icboundTS(conn)
          end if

        end do
        this%SubFaceBoundaryCounts(1) = i
    else
        ! If no connections, is boundary, defaultICBound 
        this%MassBoundarySubFace1 = defaultICBound
    end if
    
    
    ! Face 2
    count = grid%GetPotentialFaceConnectionCount(cellNumber, 2)
    this%PotentialConnectionsCount(2) = count
    if(count .gt. 0) then
        this%SubFaceCounts(2) = count
        i = 0
        do n = 1, count
          conn = grid%GetFaceConnection(cellNumber,2,n)
          this%SubFaceConn2(n) = conn

          ! RWPT: it also includes processing of mass boundaries
          ! To avoid repeating if block
          if(conn .eq. 0) then
              i = i + 1
              ! Set default value when no connection 
              this%MassBoundarySubFace2(n) = defaultICBound
          else
              ! Flow boundary, increase i
              if(iboundTS(conn) .eq. 0) i = i + 1
              ! Mass boundary
              this%MassBoundarySubFace2(n) = icboundTS(conn)
          end if
        end do
        this%SubFaceBoundaryCounts(2) = i
    else
        ! If no connections, is boundary, defaultICBound 
        this%MassBoundarySubFace2 = defaultICBound
    end if
    
    ! Face 3
    count = grid%GetPotentialFaceConnectionCount(cellNumber, 3)
    this%PotentialConnectionsCount(3) = count
    if(count .gt. 0) then
       this%SubFaceCounts(3) = count
        i = 0
        do n = 1, count
          conn = grid%GetFaceConnection(cellNumber,3,n)
          this%SubFaceConn3(n) = conn
          if(conn .eq. 0) then
              i = i + 1
              ! Set default value when no connection 
              this%MassBoundarySubFace3(n) = defaultICBound
          else
              if(iboundTS(conn) .eq. 0) i = i + 1
              ! Mass boundary
              this%MassBoundarySubFace3(n) = icboundTS(conn)
          end if
        end do
        this%SubFaceBoundaryCounts(3) = i
    else
        ! If no connections, is boundary, defaultICBound 
        this%MassBoundarySubFace3 = defaultICBound
    end if
    
    ! Face 4
    count = grid%GetPotentialFaceConnectionCount(cellNumber, 4)
    this%PotentialConnectionsCount(4) = count
    if(count .gt. 0) then
        this%SubFaceCounts(4) = count
        i = 0
        do n = 1, count
          conn = grid%GetFaceConnection(cellNumber,4,n)
          this%SubFaceConn4(n) = conn
          if(conn .eq. 0) then
              i = i + 1
              ! Set default value when no connection 
              this%MassBoundarySubFace4(n) = defaultICBound 
          else
              if(iboundTS(conn) .eq. 0) i = i + 1
              ! Mass boundary
              this%MassBoundarySubFace4(n) = icboundTS(conn)
          end if
        end do
        this%SubFaceBoundaryCounts(4) = i
    else
        ! If no connections, is boundary, defaultICBound 
        this%MassBoundarySubFace4 = defaultICBound 
    end if
    
    ! Face 5
    count = grid%GetPotentialFaceConnectionCount(cellNumber, 5)
    this%PotentialConnectionsCount(5) = count
    if(count .gt. 0) then
        this%SubFaceCounts(5) = count
        i = 0
        do n = 1, count
          conn = grid%GetFaceConnection(cellNumber,5,n)
          this%SubFaceConn5(n) = conn
          if(conn .eq. 0) then
              i = i + 1
              ! Set default value when no connection 
              this%MassBoundarySubFace5(n) =  defaultICBound
          else
              if(iboundTS(conn) .eq. 0) i = i + 1
              ! Mass boundary
              this%MassBoundarySubFace5(n) = icboundTS(conn)
          end if
        end do
        this%SubFaceBoundaryCounts(5) = i
    else
        ! If no connections, is boundary, defaultICBound 
        this%MassBoundarySubFace5 = defaultICBound
    end if
    
    ! Face 6
    count = grid%GetPotentialFaceConnectionCount(cellNumber, 6)
    this%PotentialConnectionsCount(6) = count
    if(count .gt. 0) then
        this%SubFaceCounts(6) = count
        i = 0
        do n = 1, count
          conn = grid%GetFaceConnection(cellNumber,6,n)
          this%SubFaceConn6(n) = conn
          if(conn .eq. 0) then
              i = i + 1
              ! Set default value when no connection 
              this%MassBoundarySubFace6(n) = defaultICBound
          else
              if(iboundTS(conn) .eq. 0) i = i + 1
              ! Mass boundary
              this%MassBoundarySubFace6(n) = icboundTS(conn)
          end if
        end do
        this%SubFaceBoundaryCounts(6) = i
    else
        ! If no connections, is boundary, defaultICBound 
        this%MassBoundarySubFace6 = defaultICBound
    end if
    

    ! Process face data
    ! Face 1
    count = this%SubFaceCounts(1)
    do n = 1, count
      conn = this%SubFaceConn1(n)
      if(conn .gt. 0) then
        index = grid%FindConnectionIndex(cellNumber, conn)
        this%Q1(n) = faceFlows(index)
      end if
    end do
    
    ! Face 2
    count = this%SubFaceCounts(2)
    do n = 1, count
      conn = this%SubFaceConn2(n)
      if(conn .gt. 0) then
        index = grid%FindConnectionIndex(cellNumber, conn)
        this%Q2(n) = -faceFlows(index)
      end if
    end do
    
    ! Face 3
    count = this%SubFaceCounts(3)
    do n = 1, count
      conn = this%SubFaceConn3(n)
      if(conn .gt. 0) then
        index = grid%FindConnectionIndex(cellNumber, conn)
        this%Q3(n) = faceFlows(index)
      end if
    end do
    
    ! Face 4
    count = this%SubFaceCounts(4)
    do n = 1, count
      conn = this%SubFaceConn4(n)
      if(conn .gt. 0) then
        index = grid%FindConnectionIndex(cellNumber, conn)
        this%Q4(n) = -faceFlows(index)
      end if
    end do
    
    ! Face 5
    count = this%SubFaceCounts(5)
    do n = 1, count
      conn = this%SubFaceConn5(n)
      if(conn .gt. 0) then
        index = grid%FindConnectionIndex(cellNumber, conn)
        this%Q5(n) = faceFlows(index)
      end if
    end do
    
    ! Face 6
    count = this%SubFaceCounts(6)
    do n = 1, count
      conn = this%SubFaceConn6(n)
      if(conn .gt. 0) then
        index = grid%FindConnectionIndex(cellNumber, conn)
        this%Q6(n) = -faceFlows(index)
      end if
    end do
    
    ! Process boundary flow data
    ! Face 1
    if(boundaryFlows(1) .ne. 0.0d0) then
        if(this%SubFaceBoundaryCounts(1) .gt. 0) then
          flow = boundaryFlows(1)/this%SubFaceBoundaryCounts(1)
          do n = 1, this%SubFaceCounts(1)
            conn = this%SubFaceConn1(n)
            if (conn .eq. 0) then
              this%Q1(n) = flow
            else if (iboundTS(conn) .eq. 0) then
              this%Q1(n) = flow
            end if
          end do
        else
          if(boundaryFlows(1) .gt. 0d0) then
            this%SourceFlow = this%SourceFlow + boundaryFlows(1)
          else
            this%SinkFlow = this%SinkFlow + boundaryFlows(1)
          end if
        end if
    end if

    ! Face 2
    if(boundaryFlows(2) .ne. 0.0d0) then
        if(this%SubFaceBoundaryCounts(2) .gt. 0) then
          flow = boundaryFlows(2)/this%SubFaceBoundaryCounts(2)
          do n = 1, this%SubFaceCounts(2)
            conn = this%SubFaceConn2(n)
            if (conn .eq. 0) then
              this%Q2(n) = -flow
            else if (iboundTS(conn) .eq. 0) then
              this%Q2(n) = -flow
            end if
          end do
        else
          if(boundaryFlows(2) .gt. 0d0) then
            this%SourceFlow = this%SourceFlow + boundaryFlows(2)
          else
            this%SinkFlow = this%SinkFlow + boundaryFlows(2)
          end if
        end if
    end if
    
    ! Face 3
    if(boundaryFlows(3) .ne. 0.0d0) then
        if(this%SubFaceBoundaryCounts(3) .gt. 0) then
          flow = boundaryFlows(3)/this%SubFaceBoundaryCounts(3)
          do n = 1, this%SubFaceCounts(3)
            conn = this%SubFaceConn3(n)
            if (conn .eq. 0) then
              this%Q3(n) = flow
            else if (iboundTS(conn) .eq. 0) then
              this%Q3(n) = flow
            end if
          end do
        else
          if(boundaryFlows(3) .gt. 0d0) then
            this%SourceFlow = this%SourceFlow + boundaryFlows(3)
          else
            this%SinkFlow = this%SinkFlow + boundaryFlows(3)
          end if
        end if
    end if

    ! Face 4
    if(boundaryFlows(4) .ne. 0.0d0) then
        if(this%SubFaceBoundaryCounts(4) .gt. 0) then
          flow = boundaryFlows(4)/this%SubFaceBoundaryCounts(4)
          do n = 1, this%SubFaceCounts(4)
            conn = this%SubFaceConn4(n)
            if (conn .eq. 0) then
              this%Q4(n) = -flow
            else if (iboundTS(conn) .eq. 0) then
              this%Q4(n) = -flow
            end if
          end do
        else
          if(boundaryFlows(4) .gt. 0d0) then
            this%SourceFlow = this%SourceFlow + boundaryFlows(4)
          else
            this%SinkFlow = this%SinkFlow + boundaryFlows(4)
          end if
        end if
    end if
    
    ! Face 5
    if(boundaryFlows(5) .ne. 0.0d0) then
        if(this%SubFaceBoundaryCounts(5) .gt. 0) then
          flow = boundaryFlows(5)/this%SubFaceBoundaryCounts(5)
          do n = 1, this%SubFaceCounts(5)
            conn = this%SubFaceConn5(n)
            if (conn .eq. 0) then
              this%Q5(n) = flow
            else if (iboundTS(conn) .eq. 0) then
              this%Q5(n) = flow
            end if
          end do
        else
          if(boundaryFlows(5) .gt. 0d0) then
            this%SourceFlow = this%SourceFlow + boundaryFlows(5)
          else
            this%SinkFlow = this%SinkFlow + boundaryFlows(5)
          end if
        end if
    end if

    ! Face 6
    if(boundaryFlows(6) .ne. 0.0d0) then
        if(this%SubFaceBoundaryCounts(6) .gt. 0) then
          flow = boundaryFlows(6)/this%SubFaceBoundaryCounts(6)
          do n = 1, this%SubFaceCounts(6)
            conn = this%SubFaceConn6(n)
            if (conn .eq. 0) then
              this%Q6(n) = -flow
            else if (iboundTS(conn) .eq. 0) then
              this%Q6(n) = -flow
            end if
          end do
        else
          if(boundaryFlows(6) .gt. 0d0) then
            this%SourceFlow = this%SourceFlow + boundaryFlows(6)
          else
            this%SinkFlow = this%SinkFlow + boundaryFlows(6)
          end if
        end if
    end if
    
    ! Set sub-cell row and column count
    this%SubCellRowCount = 1
    this%SubCellColumnCount = 1
    do n = 1, 6
        if(this%SubFaceCounts(n) .gt. 1) then
            this%SubCellRowCount = 2
            this%SubCellColumnCount = 2
            ! Compute internal sub-cell face flows for cells with multiple sub-cell
            call this%ComputeSubCellFlows()
            return
        end if
    end do

    ! Mark the intialized flag 
    this%initialized = .true.

  end subroutine pr_SetMassTransportDataUnstructured


  subroutine pr_SetMassTransportDataStructured(this,cellNumber,cellCount,grid,ibound, &
    iboundTS,porosity,retardation,storageFlow,sourceFlow,sinkFlow,                    &
    flowsRightFace,flowsFrontFace,flowsLowerFace,boundaryFlows, head,                 &
    layerType, zone, icboundTS, defaultICBound )
    !-----------------------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------------------
    ! Specifications
    !-----------------------------------------------------------------------------------
    implicit none
    class(ModpathCellDataType) :: this
    class(ModflowRectangularGridType),intent(in) :: grid
    integer,intent(in) :: cellNumber
    integer,intent(in) :: cellCount
    integer,intent(in) :: layerType, zone
    integer,intent(in),dimension(cellCount) :: ibound, iboundTS
    integer,intent(in),dimension(cellCount) :: icboundTS
    integer,intent(in)                      :: defaultICBound
    doubleprecision,intent(in) :: porosity, retardation, storageFlow, sourceFlow, sinkFlow
    doubleprecision,intent(in),dimension(6) :: boundaryFlows
    doubleprecision,intent(in),dimension(cellCount) :: flowsRightFace
    doubleprecision,intent(in),dimension(cellCount) :: flowsFrontFace
    doubleprecision,intent(in),dimension(cellCount) :: flowsLowerFace
    doubleprecision,intent(in) :: head
    integer :: n,count,i,conn
    doubleprecision :: flow
    !-----------------------------------------------------------------------------------
    
    ! Is this reset necessary ?
    ! All scalar variables are reassigned in the following lines after reset
    !call this%Reset()
    call this%ResetArrays()

    this%CellNumber = cellNumber
    this%Layer = grid%GetLayer(cellNumber)
    this%DX = grid%DelX(cellNumber)
    this%DY = grid%DelY(cellNumber)
    this%MinX = grid%GetLeft(cellNumber)
    this%MinY = grid%GetFront(cellNumber)
    this%Bottom = grid%Bottom(cellNumber)
    this%Top = grid%Top(cellNumber)
    this%ReducedConnectionCount = grid%GetJaCellConnectionsCount(cellNumber)
    !this%ReducedConnectionCount = grid%GetReducedCellConnectionCount(cellNumber)
    
    ! Assign property data
    this%Zone = zone
    this%LayerType = layerType
    this%Head = head
    this%Ibound = ibound(cellNumber)
    this%IboundTS = iboundTS(cellNumber)
    this%Porosity = porosity
    this%Retardation = retardation
    this%SourceFlow = sourceFlow
    this%SinkFlow = sinkFlow
    this%StorageFlow = storageFlow
 

    ! RWPT 
    this%ICBoundTS = icboundTS(cellNumber)


    ! Process face flow data
    ! Face 1
    count = grid%GetPotentialFaceConnectionCount(cellNumber, 1)
    this%PotentialConnectionsCount(1) = count
    if(count .gt. 0) then
        this%SubFaceCounts(1) = count
        i = 0
        do n = 1, count
          conn = grid%GetFaceConnection(cellNumber,1,n)
          this%SubFaceConn1(n) = conn
          if(conn .eq. 0) then
              i = i + 1
              ! Set default value when no connection 
              this%MassBoundarySubFace1(n) = defaultICBound
          else
               if(iboundTS(conn) .eq. 0) i = i + 1
               ! Mass boundary
               this%MassBoundarySubFace1(n) = icboundTS(conn)
          end if
        end do
        this%SubFaceBoundaryCounts(1) = i
    else
        ! If no connections, is boundary, defaultICBound 
        this%MassBoundarySubFace1 = defaultICBound
    end if
    
    
    ! Face 2
    count = grid%GetPotentialFaceConnectionCount(cellNumber, 2)
    this%PotentialConnectionsCount(2) = count
    if(count .gt. 0) then
        this%SubFaceCounts(2) = count
        i = 0
        do n = 1, count
          conn = grid%GetFaceConnection(cellNumber,2,n)
          this%SubFaceConn2(n) = conn
          if(conn .eq. 0) then
              i = i + 1
              ! Set default value when no connection 
              this%MassBoundarySubFace2(n) = defaultICBound
          else
              if(iboundTS(conn) .eq. 0) i = i + 1
              ! Mass boundary
              this%MassBoundarySubFace2(n) = icboundTS(conn)
          end if
        end do
        this%SubFaceBoundaryCounts(2) = i
    else
        ! If no connections, is boundary, defaultICBound 
        this%MassBoundarySubFace2 = defaultICBound
    end if
    
    ! Face 3
    count = grid%GetPotentialFaceConnectionCount(cellNumber, 3)
    this%PotentialConnectionsCount(3) = count
    if(count .gt. 0) then
       this%SubFaceCounts(3) = count
        i = 0
        do n = 1, count
          conn = grid%GetFaceConnection(cellNumber,3,n)
          this%SubFaceConn3(n) = conn
          if(conn .eq. 0) then
              i = i + 1
              ! Set default value when no connection 
              this%MassBoundarySubFace3(n) = defaultICBound
          else
              if(iboundTS(conn) .eq. 0) i = i + 1
              ! Mass boundary
              this%MassBoundarySubFace3(n) = icboundTS(conn)
          end if
        end do
        this%SubFaceBoundaryCounts(3) = i
    else
        ! If no connections, is boundary, defaultICBound 
        this%MassBoundarySubFace3 = defaultICBound
    end if
    
    ! Face 4
    count = grid%GetPotentialFaceConnectionCount(cellNumber, 4)
    this%PotentialConnectionsCount(4) = count
    if(count .gt. 0) then
        this%SubFaceCounts(4) = count
        i = 0
        do n = 1, count
          conn = grid%GetFaceConnection(cellNumber,4,n)
          this%SubFaceConn4(n) = conn
          if(conn .eq. 0) then
              i = i + 1
              ! Set default value when no connection 
              this%MassBoundarySubFace4(n) = defaultICBound 
          else
              if(iboundTS(conn) .eq. 0) i = i + 1
              ! Mass boundary
              this%MassBoundarySubFace4(n) = icboundTS(conn)
          end if
        end do
        this%SubFaceBoundaryCounts(4) = i
    else
        ! If no connections, is boundary, defaultICBound 
        this%MassBoundarySubFace4 = defaultICBound 
    end if
    
    ! Face 5
    count = grid%GetPotentialFaceConnectionCount(cellNumber, 5)
    this%PotentialConnectionsCount(5) = count
    if(count .gt. 0) then
        this%SubFaceCounts(5) = count
        i = 0
        do n = 1, count
          conn = grid%GetFaceConnection(cellNumber,5,n)
          this%SubFaceConn5(n) = conn
          if(conn .eq. 0) then
              i = i + 1
              ! Set default value when no connection 
              this%MassBoundarySubFace5(n) =  defaultICBound
          else
              if(iboundTS(conn) .eq. 0) i = i + 1
              ! Mass boundary
              this%MassBoundarySubFace5(n) = icboundTS(conn)
          end if
        end do
        this%SubFaceBoundaryCounts(5) = i
    else
        ! If no connections, is boundary, defaultICBound 
        this%MassBoundarySubFace5 = defaultICBound
    end if
    
    ! Face 6
    count = grid%GetPotentialFaceConnectionCount(cellNumber, 6)
    this%PotentialConnectionsCount(6) = count
    if(count .gt. 0) then
        this%SubFaceCounts(6) = count
        i = 0
        do n = 1, count
          conn = grid%GetFaceConnection(cellNumber,6,n)
          this%SubFaceConn6(n) = conn
          if(conn .eq. 0) then
              i = i + 1
              ! Set default value when no connection 
              this%MassBoundarySubFace6(n) =  defaultICBound
          else
              if(iboundTS(conn) .eq. 0) i = i + 1
              ! Mass boundary
              this%MassBoundarySubFace6(n) = icboundTS(conn)
          end if
        end do
        this%SubFaceBoundaryCounts(6) = i
    else
        ! If no connections, is boundary, defaultICBound 
        this%MassBoundarySubFace6 = defaultICBound
    end if
    

    ! Process face data
    ! Face 1
    if(this%SubFaceCounts(1) .eq. 1) then
      conn = this%SubFaceConn1(1)
      if(conn .gt. 0) then
        this%Q1(1) = flowsRightFace(conn)
      end if
    end if
    
    ! Face 2
    if(this%SubFaceCounts(2) .eq. 1) then
      this%Q2(1) = flowsRightFace(cellNumber)
    end if
    
    ! Face 3
    if(this%SubFaceCounts(3) .eq. 1) then
      this%Q3(1) = -flowsFrontFace(cellNumber)
    end if
    
    ! Face 4
    if(this%SubFaceCounts(4) .eq. 1) then
      conn = this%SubFaceConn4(1)
      if(conn .gt. 0) then
        this%Q4(1) = -flowsFrontFace(conn)
      end if
    end if
    
    ! Face 5
    if(this%SubFaceCounts(5) .eq. 1) then
      this%Q5(1) = -flowsLowerFace(cellNumber)
    end if
    
    ! Face 6
    if(this%SubFaceCounts(6) .eq. 1) then
      conn = this%SubFaceConn6(1)
      if(conn .gt. 0) then
        this%Q6(1) = -flowsLowerFace(conn)
      end if
    end if
    
    ! Process boundary flow data
    ! Face 1
    if(boundaryFlows(1) .ne. 0.0d0) then
        if(this%SubFaceBoundaryCounts(1) .gt. 0) then
          flow = boundaryFlows(1)/this%SubFaceBoundaryCounts(1)
          do n = 1, this%SubFaceCounts(1)
            conn = this%SubFaceConn1(n)
            if (conn .eq. 0) then
              this%Q1(n) = flow
            else if (iboundTS(conn) .eq. 0) then
              this%Q1(n) = flow
            end if
          end do
        else
          if(boundaryFlows(1) .gt. 0d0) then
            this%SourceFlow = this%SourceFlow + boundaryFlows(1)
          else
            this%SinkFlow = this%SinkFlow + boundaryFlows(1)
          end if
        end if
    end if

    ! Face 2
    if(boundaryFlows(2) .ne. 0.0d0) then
        if(this%SubFaceBoundaryCounts(2) .gt. 0) then
          flow = boundaryFlows(2)/this%SubFaceBoundaryCounts(2)
          do n = 1, this%SubFaceCounts(2)
            conn = this%SubFaceConn2(n)
            if (conn .eq. 0) then
              this%Q2(n) = -flow
            else if (iboundTS(conn) .eq. 0) then
              this%Q2(n) = -flow
            end if
          end do
        else
          if(boundaryFlows(2) .gt. 0d0) then
            this%SourceFlow = this%SourceFlow + boundaryFlows(2)
          else
            this%SinkFlow = this%SinkFlow + boundaryFlows(2)
          end if
        end if
    end if
    
    ! Face 3
    if(boundaryFlows(3) .ne. 0.0d0) then
        if(this%SubFaceBoundaryCounts(3) .gt. 0) then
          flow = boundaryFlows(3)/this%SubFaceBoundaryCounts(3)
          do n = 1, this%SubFaceCounts(3)
            conn = this%SubFaceConn3(n)
            if (conn .eq. 0) then
              this%Q3(n) = flow
            else if (iboundTS(conn) .eq. 0) then
              this%Q3(n) = flow
            end if
          end do
        else
          if(boundaryFlows(3) .gt. 0d0) then
            this%SourceFlow = this%SourceFlow + boundaryFlows(3)
          else
            this%SinkFlow = this%SinkFlow + boundaryFlows(3)
          end if
        end if
    end if

    ! Face 4
    if(boundaryFlows(4) .ne. 0.0d0) then
        if(this%SubFaceBoundaryCounts(4) .gt. 0) then
          flow = boundaryFlows(4)/this%SubFaceBoundaryCounts(4)
          do n = 1, this%SubFaceCounts(4)
            conn = this%SubFaceConn4(n)
            if (conn .eq. 0) then
              this%Q4(n) = -flow
            else if (iboundTS(conn) .eq. 0) then
              this%Q4(n) = -flow
            end if
          end do
        else
          if(boundaryFlows(4) .gt. 0d0) then
            this%SourceFlow = this%SourceFlow + boundaryFlows(4)
          else
            this%SinkFlow = this%SinkFlow + boundaryFlows(4)
          end if
        end if
    end if
    
    ! Face 5
    if(boundaryFlows(5) .ne. 0.0d0) then
        if(this%SubFaceBoundaryCounts(5) .gt. 0) then
          flow = boundaryFlows(5)/this%SubFaceBoundaryCounts(5)
          do n = 1, this%SubFaceCounts(5)
            conn = this%SubFaceConn5(n)
            if (conn .eq. 0) then
              this%Q5(n) = flow
            else if (iboundTS(conn) .eq. 0) then
              this%Q5(n) = flow
            end if
          end do
        else
          if(boundaryFlows(5) .gt. 0d0) then
            this%SourceFlow = this%SourceFlow + boundaryFlows(5)
          else
            this%SinkFlow = this%SinkFlow + boundaryFlows(5)
          end if
        end if
    end if

    ! Face 6
    if(boundaryFlows(6) .ne. 0.0d0) then
        if(this%SubFaceBoundaryCounts(6) .gt. 0) then
          flow = boundaryFlows(6)/this%SubFaceBoundaryCounts(6)
          do n = 1, this%SubFaceCounts(6)
            conn = this%SubFaceConn6(n)
            if (conn .eq. 0) then
              this%Q6(n) = -flow
            else if (iboundTS(conn) .eq. 0) then
              this%Q6(n) = -flow
            end if
          end do
        else
          if(boundaryFlows(6) .gt. 0d0) then
            this%SourceFlow = this%SourceFlow + boundaryFlows(6)
          else
            this%SinkFlow = this%SinkFlow + boundaryFlows(6)
          end if
        end if
    end if
    
    ! Set sub-cell row and column count. For structured grids, all cells have only 1 sub-cell.
    this%SubCellRowCount = 1
    this%SubCellColumnCount = 1

    ! Mark the intialized flag 
    this%initialized = .true.

  end subroutine pr_SetMassTransportDataStructured



! ORIGINAL
!------------------------------------------
  subroutine pr_SetDataUnstructured(this,cellNumber,cellCount,reducedConnectionCount,grid,&
    ibound,iboundTS,porosity,retardation,storageFlow,sourceFlow,sinkFlow,       &
    faceFlows,boundaryFlows, head, layerType, zone)
  implicit none
  class(ModpathCellDataType) :: this
  class(ModflowRectangularGridType),intent(in) :: grid
  integer,intent(in) :: cellNumber,cellCount
  integer,intent(in) :: reducedConnectionCount
  integer,intent(in) :: layerType, zone
  integer,intent(in),dimension(cellCount) :: ibound, iboundTS
  doubleprecision,intent(in) :: porosity, retardation, storageFlow, sourceFlow, sinkFlow
  doubleprecision,intent(in),dimension(6) :: boundaryFlows
  doubleprecision,intent(in),dimension(reducedConnectionCount) :: faceFlows
  doubleprecision,intent(in) :: head
  integer :: n,index,count,i,conn
  doubleprecision :: flow
  
  call this%Reset()

  this%CellNumber = cellNumber
  this%Layer = grid%GetLayer(cellNumber)
  this%DX = grid%DelX(cellNumber)
  this%DY = grid%DelY(cellNumber)
  this%MinX = grid%GetLeft(cellNumber)
  this%MinY = grid%GetFront(cellNumber)
  this%Bottom = grid%Bottom(cellNumber)
  this%Top = grid%Top(cellNumber)
  this%ReducedConnectionCount = grid%GetJaCellConnectionsCount(cellNumber)
  !this%ReducedConnectionCount = grid%GetReducedCellConnectionCount(cellNumber)
  
  ! Assign property data
  this%Zone = zone
  this%LayerType = layerType
  this%Head = head
  this%Ibound = ibound(cellNumber)
  this%IboundTS = iboundTS(cellNumber)
  this%Porosity = porosity
  this%Retardation = retardation
  this%SourceFlow = sourceFlow
  this%SinkFlow = sinkFlow
  this%StorageFlow = storageFlow
  
  ! Process face flow data
  ! Face 1
  count = grid%GetPotentialFaceConnectionCount(cellNumber, 1)
  this%PotentialConnectionsCount(1) = count
  if(count .gt. 0) then
      this%SubFaceCounts(1) = count
      i = 0
      do n = 1, count
        conn = grid%GetFaceConnection(cellNumber,1,n)
        this%SubFaceConn1(n) = conn
        if(conn .eq. 0) then
            i = i + 1
        else
            if(iboundTS(conn) .eq. 0) i = i + 1
        end if
      end do
      this%SubFaceBoundaryCounts(1) = i
  end if
  
  
  ! Face 2
  count = grid%GetPotentialFaceConnectionCount(cellNumber, 2)
  this%PotentialConnectionsCount(2) = count
  if(count .gt. 0) then
      this%SubFaceCounts(2) = count
      i = 0
      do n = 1, count
        conn = grid%GetFaceConnection(cellNumber,2,n)
        this%SubFaceConn2(n) = conn
        if(conn .eq. 0) then
            i = i + 1
        else
            if(iboundTS(conn) .eq. 0) i = i + 1
        end if
      end do
      this%SubFaceBoundaryCounts(2) = i
  end if
  
  ! Face 3
  count = grid%GetPotentialFaceConnectionCount(cellNumber, 3)
  this%PotentialConnectionsCount(3) = count
  if(count .gt. 0) then
     this%SubFaceCounts(3) = count
      i = 0
      do n = 1, count
        conn = grid%GetFaceConnection(cellNumber,3,n)
        this%SubFaceConn3(n) = conn
        if(conn .eq. 0) then
            i = i + 1
        else
            if(iboundTS(conn) .eq. 0) i = i + 1
        end if
      end do
      this%SubFaceBoundaryCounts(3) = i
  end if
  
  ! Face 4
  count = grid%GetPotentialFaceConnectionCount(cellNumber, 4)
  this%PotentialConnectionsCount(4) = count
  if(count .gt. 0) then
      this%SubFaceCounts(4) = count
      i = 0
      do n = 1, count
        conn = grid%GetFaceConnection(cellNumber,4,n)
        this%SubFaceConn4(n) = conn
        if(conn .eq. 0) then
            i = i + 1
        else
            if(iboundTS(conn) .eq. 0) i = i + 1
        end if
      end do
      this%SubFaceBoundaryCounts(4) = i
  end if
  
  ! Face 5
  count = grid%GetPotentialFaceConnectionCount(cellNumber, 5)
  this%PotentialConnectionsCount(5) = count
  if(count .gt. 0) then
      this%SubFaceCounts(5) = count
      i = 0
      do n = 1, count
        conn = grid%GetFaceConnection(cellNumber,5,n)
        this%SubFaceConn5(n) = conn
        if(conn .eq. 0) then
            i = i + 1
        else
             if(iboundTS(conn) .eq. 0) i = i + 1
        end if
      end do
      this%SubFaceBoundaryCounts(5) = i
  end if
  
  ! Face 6
  count = grid%GetPotentialFaceConnectionCount(cellNumber, 6)
  this%PotentialConnectionsCount(6) = count
  if(count .gt. 0) then
      this%SubFaceCounts(6) = count
      i = 0
      do n = 1, count
        conn = grid%GetFaceConnection(cellNumber,6,n)
        this%SubFaceConn6(n) = conn
        if(conn .eq. 0) then
            i = i + 1
        else
            if(iboundTS(conn) .eq. 0) i = i + 1
        end if
      end do
      this%SubFaceBoundaryCounts(6) = i
  end if
  

  ! Process face data
  ! Face 1
  count = this%SubFaceCounts(1)
  do n = 1, count
    conn = this%SubFaceConn1(n)
    if(conn .gt. 0) then
      index = grid%FindConnectionIndex(cellNumber, conn)
      this%Q1(n) = faceFlows(index)
    end if
  end do
  
  ! Face 2
  count = this%SubFaceCounts(2)
  do n = 1, count
    conn = this%SubFaceConn2(n)
    if(conn .gt. 0) then
      index = grid%FindConnectionIndex(cellNumber, conn)
      this%Q2(n) = -faceFlows(index)
    end if
  end do
  
  ! Face 3
  count = this%SubFaceCounts(3)
  do n = 1, count
    conn = this%SubFaceConn3(n)
    if(conn .gt. 0) then
      index = grid%FindConnectionIndex(cellNumber, conn)
      this%Q3(n) = faceFlows(index)
    end if
  end do
  
  ! Face 4
  count = this%SubFaceCounts(4)
  do n = 1, count
    conn = this%SubFaceConn4(n)
    if(conn .gt. 0) then
      index = grid%FindConnectionIndex(cellNumber, conn)
      this%Q4(n) = -faceFlows(index)
    end if
  end do
  
  ! Face 5
  count = this%SubFaceCounts(5)
  do n = 1, count
    conn = this%SubFaceConn5(n)
    if(conn .gt. 0) then
      index = grid%FindConnectionIndex(cellNumber, conn)
      this%Q5(n) = faceFlows(index)
    end if
  end do
  
  ! Face 6
  count = this%SubFaceCounts(6)
  do n = 1, count
    conn = this%SubFaceConn6(n)
    if(conn .gt. 0) then
      index = grid%FindConnectionIndex(cellNumber, conn)
      this%Q6(n) = -faceFlows(index)
    end if
  end do
  
  ! Process boundary flow data
  ! Face 1
  if(boundaryFlows(1) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(1) .gt. 0) then
        flow = boundaryFlows(1)/this%SubFaceBoundaryCounts(1)
        do n = 1, this%SubFaceCounts(1)
          conn = this%SubFaceConn1(n)
          if (conn .eq. 0) then
            this%Q1(n) = flow
          else if (iboundTS(conn) .eq. 0) then
            this%Q1(n) = flow
          end if
        end do
      else
        if(boundaryFlows(1) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(1)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(1)
        end if
      end if
  end if

  ! Face 2
  if(boundaryFlows(2) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(2) .gt. 0) then
        flow = boundaryFlows(2)/this%SubFaceBoundaryCounts(2)
        do n = 1, this%SubFaceCounts(2)
          conn = this%SubFaceConn2(n)
          if (conn .eq. 0) then
            this%Q2(n) = -flow
          else if (iboundTS(conn) .eq. 0) then
            this%Q2(n) = -flow
          end if
        end do
      else
        if(boundaryFlows(2) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(2)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(2)
        end if
      end if
  end if
  
  ! Face 3
  if(boundaryFlows(3) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(3) .gt. 0) then
        flow = boundaryFlows(3)/this%SubFaceBoundaryCounts(3)
        do n = 1, this%SubFaceCounts(3)
          conn = this%SubFaceConn3(n)
          if (conn .eq. 0) then
            this%Q3(n) = flow
          else if (iboundTS(conn) .eq. 0) then
            this%Q3(n) = flow
          end if
        end do
      else
        if(boundaryFlows(3) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(3)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(3)
        end if
      end if
  end if

  ! Face 4
  if(boundaryFlows(4) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(4) .gt. 0) then
        flow = boundaryFlows(4)/this%SubFaceBoundaryCounts(4)
        do n = 1, this%SubFaceCounts(4)
          conn = this%SubFaceConn4(n)
          if (conn .eq. 0) then
            this%Q4(n) = -flow
          else if (iboundTS(conn) .eq. 0) then
            this%Q4(n) = -flow
          end if
        end do
      else
        if(boundaryFlows(4) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(4)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(4)
        end if
      end if
  end if
  
  ! Face 5
  if(boundaryFlows(5) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(5) .gt. 0) then
        flow = boundaryFlows(5)/this%SubFaceBoundaryCounts(5)
        do n = 1, this%SubFaceCounts(5)
          conn = this%SubFaceConn5(n)
          if (conn .eq. 0) then
            this%Q5(n) = flow
          else if (iboundTS(conn) .eq. 0) then
            this%Q5(n) = flow
          end if
        end do
      else
        if(boundaryFlows(5) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(5)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(5)
        end if
      end if
  end if

  ! Face 6
  if(boundaryFlows(6) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(6) .gt. 0) then
        flow = boundaryFlows(6)/this%SubFaceBoundaryCounts(6)
        do n = 1, this%SubFaceCounts(6)
          conn = this%SubFaceConn6(n)
          if (conn .eq. 0) then
            this%Q6(n) = -flow
          else if (iboundTS(conn) .eq. 0) then
            this%Q6(n) = -flow
          end if
        end do
      else
        if(boundaryFlows(6) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(6)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(6)
        end if
      end if
  end if
  
  ! Set sub-cell row and column count
  this%SubCellRowCount = 1
  this%SubCellColumnCount = 1
  do n = 1, 6
      if(this%SubFaceCounts(n) .gt. 1) then
          this%SubCellRowCount = 2
          this%SubCellColumnCount = 2
          ! Compute internal sub-cell face flows for cells with multiple sub-cell
          call this%ComputeSubCellFlows()
          return
      end if
  end do

  end subroutine pr_SetDataUnstructured

!------------------------------------------
  subroutine pr_SetDataStructured(this,cellNumber,cellCount,grid,ibound,        &
    iboundTS,porosity,retardation,storageFlow,sourceFlow,sinkFlow,              &
    flowsRightFace,flowsFrontFace,flowsLowerFace,boundaryFlows, head,           &
    layerType, zone)
  implicit none
  class(ModpathCellDataType) :: this
  class(ModflowRectangularGridType),intent(in) :: grid
  integer,intent(in) :: cellNumber
  integer,intent(in) :: cellCount
  integer,intent(in) :: layerType, zone
  integer,intent(in),dimension(cellCount) :: ibound, iboundTS
  doubleprecision,intent(in) :: porosity, retardation, storageFlow, sourceFlow, sinkFlow
  doubleprecision,intent(in),dimension(6) :: boundaryFlows
  doubleprecision,intent(in),dimension(cellCount) :: flowsRightFace
  doubleprecision,intent(in),dimension(cellCount) :: flowsFrontFace
  doubleprecision,intent(in),dimension(cellCount) :: flowsLowerFace
  doubleprecision,intent(in) :: head
  integer :: n,count,i,conn
  doubleprecision :: flow
  
  call this%Reset()

  this%CellNumber = cellNumber
  this%Layer = grid%GetLayer(cellNumber)
  this%DX = grid%DelX(cellNumber)
  this%DY = grid%DelY(cellNumber)
  this%MinX = grid%GetLeft(cellNumber)
  this%MinY = grid%GetFront(cellNumber)
  this%Bottom = grid%Bottom(cellNumber)
  this%Top = grid%Top(cellNumber)
  this%ReducedConnectionCount = grid%GetJaCellConnectionsCount(cellNumber)
  !this%ReducedConnectionCount = grid%GetReducedCellConnectionCount(cellNumber)
  
  ! Assign property data
  this%Zone = zone
  this%LayerType = layerType
  this%Head = head
  this%Ibound = ibound(cellNumber)
  this%IboundTS = iboundTS(cellNumber)
  this%Porosity = porosity
  this%Retardation = retardation
  this%SourceFlow = sourceFlow
  this%SinkFlow = sinkFlow
  this%StorageFlow = storageFlow
  
  ! Process face flow data
  ! Face 1
  count = grid%GetPotentialFaceConnectionCount(cellNumber, 1)
  this%PotentialConnectionsCount(1) = count
  if(count .gt. 0) then
      this%SubFaceCounts(1) = count
      i = 0
      do n = 1, count
        conn = grid%GetFaceConnection(cellNumber,1,n)
        this%SubFaceConn1(n) = conn
        if(conn .eq. 0) then
            i = i + 1
        else
             if(iboundTS(conn) .eq. 0) i = i + 1
        end if
      end do
      this%SubFaceBoundaryCounts(1) = i
  end if
  
  
  ! Face 2
  count = grid%GetPotentialFaceConnectionCount(cellNumber, 2)
  this%PotentialConnectionsCount(2) = count
  if(count .gt. 0) then
      this%SubFaceCounts(2) = count
      i = 0
      do n = 1, count
        conn = grid%GetFaceConnection(cellNumber,2,n)
        this%SubFaceConn2(n) = conn
        if(conn .eq. 0) then
            i = i + 1
        else
            if(iboundTS(conn) .eq. 0) i = i + 1
        end if
      end do
      this%SubFaceBoundaryCounts(2) = i
  end if
  
  ! Face 3
  count = grid%GetPotentialFaceConnectionCount(cellNumber, 3)
  this%PotentialConnectionsCount(3) = count
  if(count .gt. 0) then
     this%SubFaceCounts(3) = count
      i = 0
      do n = 1, count
        conn = grid%GetFaceConnection(cellNumber,3,n)
        this%SubFaceConn3(n) = conn
        if(conn .eq. 0) then
            i = i + 1
        else
            if(iboundTS(conn) .eq. 0) i = i + 1
        end if
      end do
      this%SubFaceBoundaryCounts(3) = i
  end if
  
  ! Face 4
  count = grid%GetPotentialFaceConnectionCount(cellNumber, 4)
  this%PotentialConnectionsCount(4) = count
  if(count .gt. 0) then
      this%SubFaceCounts(4) = count
      i = 0
      do n = 1, count
        conn = grid%GetFaceConnection(cellNumber,4,n)
        this%SubFaceConn4(n) = conn
        if(conn .eq. 0) then
            i = i + 1
        else
            if(iboundTS(conn) .eq. 0) i = i + 1
        end if
      end do
      this%SubFaceBoundaryCounts(4) = i
  end if
  
  ! Face 5
  count = grid%GetPotentialFaceConnectionCount(cellNumber, 5)
  this%PotentialConnectionsCount(5) = count
  if(count .gt. 0) then
      this%SubFaceCounts(5) = count
      i = 0
      do n = 1, count
        conn = grid%GetFaceConnection(cellNumber,5,n)
        this%SubFaceConn5(n) = conn
        if(conn .eq. 0) then
            i = i + 1
        else
            if(iboundTS(conn) .eq. 0) i = i + 1
        end if
      end do
      this%SubFaceBoundaryCounts(5) = i
  end if
  
  ! Face 6
  count = grid%GetPotentialFaceConnectionCount(cellNumber, 6)
  this%PotentialConnectionsCount(6) = count
  if(count .gt. 0) then
      this%SubFaceCounts(6) = count
      i = 0
      do n = 1, count
        conn = grid%GetFaceConnection(cellNumber,6,n)
        this%SubFaceConn6(n) = conn
        if(conn .eq. 0) then
            i = i + 1
        else
            if(iboundTS(conn) .eq. 0) i = i + 1
        end if
      end do
      this%SubFaceBoundaryCounts(6) = i
  end if
  

  ! Process face data
  ! Face 1
  if(this%SubFaceCounts(1) .eq. 1) then
    conn = this%SubFaceConn1(1)
    if(conn .gt. 0) then
      this%Q1(1) = flowsRightFace(conn)
    end if
  end if
  
  ! Face 2
  if(this%SubFaceCounts(2) .eq. 1) then
    this%Q2(1) = flowsRightFace(cellNumber)
  end if
  
  ! Face 3
  if(this%SubFaceCounts(3) .eq. 1) then
    this%Q3(1) = -flowsFrontFace(cellNumber)
  end if
  
  ! Face 4
  if(this%SubFaceCounts(4) .eq. 1) then
    conn = this%SubFaceConn4(1)
    if(conn .gt. 0) then
      this%Q4(1) = -flowsFrontFace(conn)
    end if
  end if
  
  ! Face 5
  if(this%SubFaceCounts(5) .eq. 1) then
    this%Q5(1) = -flowsLowerFace(cellNumber)
  end if
  
  ! Face 6
  if(this%SubFaceCounts(6) .eq. 1) then
    conn = this%SubFaceConn6(1)
    if(conn .gt. 0) then
      this%Q6(1) = -flowsLowerFace(conn)
    end if
  end if
  
  ! Process boundary flow data
  ! Face 1
  if(boundaryFlows(1) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(1) .gt. 0) then
        flow = boundaryFlows(1)/this%SubFaceBoundaryCounts(1)
        do n = 1, this%SubFaceCounts(1)
          conn = this%SubFaceConn1(n)
          if (conn .eq. 0) then
            this%Q1(n) = flow
          else if (iboundTS(conn) .eq. 0) then
            this%Q1(n) = flow
          end if
        end do
      else
        if(boundaryFlows(1) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(1)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(1)
        end if
      end if
  end if

  ! Face 2
  if(boundaryFlows(2) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(2) .gt. 0) then
        flow = boundaryFlows(2)/this%SubFaceBoundaryCounts(2)
        do n = 1, this%SubFaceCounts(2)
          conn = this%SubFaceConn2(n)
          if (conn .eq. 0) then
            this%Q2(n) = -flow
          else if (iboundTS(conn) .eq. 0) then
            this%Q2(n) = -flow
          end if
        end do
      else
        if(boundaryFlows(2) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(2)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(2)
        end if
      end if
  end if
  
  ! Face 3
  if(boundaryFlows(3) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(3) .gt. 0) then
        flow = boundaryFlows(3)/this%SubFaceBoundaryCounts(3)
        do n = 1, this%SubFaceCounts(3)
          conn = this%SubFaceConn3(n)
          if (conn .eq. 0) then
            this%Q3(n) = flow
          else if (iboundTS(conn) .eq. 0) then
            this%Q3(n) = flow
          end if
        end do
      else
        if(boundaryFlows(3) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(3)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(3)
        end if
      end if
  end if

  ! Face 4
  if(boundaryFlows(4) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(4) .gt. 0) then
        flow = boundaryFlows(4)/this%SubFaceBoundaryCounts(4)
        do n = 1, this%SubFaceCounts(4)
          conn = this%SubFaceConn4(n)
          if (conn .eq. 0) then
            this%Q4(n) = -flow
          else if (iboundTS(conn) .eq. 0) then
            this%Q4(n) = -flow
          end if
        end do
      else
        if(boundaryFlows(4) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(4)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(4)
        end if
      end if
  end if
  
  ! Face 5
  if(boundaryFlows(5) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(5) .gt. 0) then
        flow = boundaryFlows(5)/this%SubFaceBoundaryCounts(5)
        do n = 1, this%SubFaceCounts(5)
          conn = this%SubFaceConn5(n)
          if (conn .eq. 0) then
            this%Q5(n) = flow
          else if (iboundTS(conn) .eq. 0) then
            this%Q5(n) = flow
          end if
        end do
      else
        if(boundaryFlows(5) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(5)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(5)
        end if
      end if
  end if

  ! Face 6
  if(boundaryFlows(6) .ne. 0.0d0) then
      if(this%SubFaceBoundaryCounts(6) .gt. 0) then
        flow = boundaryFlows(6)/this%SubFaceBoundaryCounts(6)
        do n = 1, this%SubFaceCounts(6)
          conn = this%SubFaceConn6(n)
          if (conn .eq. 0) then
            this%Q6(n) = -flow
          else if (iboundTS(conn) .eq. 0) then
            this%Q6(n) = -flow
          end if
        end do
      else
        if(boundaryFlows(6) .gt. 0d0) then
          this%SourceFlow = this%SourceFlow + boundaryFlows(6)
        else
          this%SinkFlow = this%SinkFlow + boundaryFlows(6)
        end if
      end if
  end if
  
  ! Set sub-cell row and column count. For structured grids, all cells have only 1 sub-cell.
  this%SubCellRowCount = 1
  this%SubCellColumnCount = 1

  end subroutine pr_SetDataStructured

!------------------------------------------
!  function pr_GetSubFaceBoundaryCount(this,faceNumber) result(count)
!  implicit none
!  class(ModpathCellDataType) :: this
!  integer,intent(in) :: faceNumber
!  integer :: count
!  
!  end function

!------------------------------------------
  function pr_FindConnectionNumberIndex(cellNumber,connectionList,arraySize,connectionCount) result(index)
  implicit none
  integer,intent(in) :: cellNumber,connectionCount,arraySize
  integer,intent(in),dimension(connectionCount) :: connectionList
  integer :: index,n
  
  index = 0
  do n = 1,arraySize
    if(connectionList(n) .eq. cellNumber) then
      index = n
      exit
    end if
  end do
  
  end function pr_FindConnectionNumberIndex

!------------------------------------------
  function pr_GetReducedConnectionCount(this) result(fval)
  implicit none
  class(ModpathCellDataType) :: this
  integer :: fval
  
  fval = this%ReducedConnectionCount

  end function pr_GetReducedConnectionCount

!------------------------------------------
  function pr_GetSubFaceCount(this,faceNumber) result(count)
  implicit none
  class(ModpathCellDataType) :: this
  integer,intent(in) :: faceNumber
  integer :: count
  
  count = this%SubFaceCounts(faceNumber)
  
  end function pr_GetSubFaceCount

!------------------------------------------
  function pr_GetSubCellRowCount(this) result(fval)
  implicit none
  class(ModpathCellDataType) :: this
  integer :: fval
  
  fval = this%SubCellRowCount

  end function pr_GetSubCellRowCount

!------------------------------------------
  function pr_GetSubCellColumnCount(this) result(fval)
  implicit none
  class(ModpathCellDataType) :: this
  integer :: fval
  
  fval = this%SubCellColumnCount

  end function pr_GetSubCellColumnCount

!------------------------------------------
  function pr_GetSubCellCount(this) result(fval)
  implicit none
  class(ModpathCellDataType) :: this
  integer :: fval
  
  fval = this%SubCellColumnCount * this%SubCellRowCount

  end function pr_GetSubCellCount
  
!------------------------------------------
  function pr_GetFaceConnection(this,faceNumber,subFaceNumber) result(conn)
  implicit none
  class(ModpathCellDataType) :: this
  integer,intent(in) :: faceNumber,subFaceNumber
  integer :: conn
  
  conn = 0
  if(subFaceNumber .lt. 1) return
  
  select case (faceNumber)
    case (1)
      if(subFaceNumber .le. this%SubFaceCounts(1)) conn = this%SubFaceConn1(subFaceNumber)
    case (2)
      if(subFaceNumber .le. this%SubFaceCounts(2)) conn = this%SubFaceConn2(subFaceNumber)
    case (3)
      if(subFaceNumber .le. this%SubFaceCounts(3)) conn = this%SubFaceConn3(subFaceNumber)
    case (4)
      if(subFaceNumber .le. this%SubFaceCounts(4)) conn = this%SubFaceConn4(subFaceNumber)
    case (5)
      if(subFaceNumber .le. this%SubFaceCounts(5)) conn = this%SubFaceConn5(subFaceNumber)
    case (6)
      if(subFaceNumber .le. this%SubFaceCounts(6)) conn = this%SubFaceConn6(subFaceNumber)
    case default
      conn = 0
  end select
  
  end function pr_GetFaceConnection
  
!------------------------------------------
  function pr_GetFaceFlow(this,faceNumber,subFaceNumber) result(flow)
  implicit none
  class(ModpathCellDataType) :: this
  integer,intent(in) :: faceNumber,subFaceNumber
  doubleprecision :: flow 
  
  flow = 0.0d0
  if(subFaceNumber .lt. 1) return
  
  select case (faceNumber)
    case (1)
      if(subFaceNumber .le. this%SubFaceCounts(1)) flow = this%Q1(subFaceNumber)
    case (2)
      if(subFaceNumber .le. this%SubFaceCounts(2)) flow = this%Q2(subFaceNumber)
    case (3)
      if(subFaceNumber .le. this%SubFaceCounts(3)) flow = this%Q3(subFaceNumber)
    case (4)
      if(subFaceNumber .le. this%SubFaceCounts(4)) flow = this%Q4(subFaceNumber)
    case (5)
      if(subFaceNumber .le. this%SubFaceCounts(5)) flow = this%Q5(subFaceNumber)
    case (6)
      if(subFaceNumber .le. this%SubFaceCounts(6)) flow = this%Q6(subFaceNumber)
    case default
      flow = 0.0d0
  end select
  
  end function pr_GetFaceFlow
  
!------------------------------------------
  subroutine pr_SetSubCellFlows(this,subFlow1,subFlow2,subFlow3,subFlow4)
  implicit none
  class(ModpathCellDataType) :: this
  doubleprecision :: subFlow1,subFlow2,subFlow3,subFlow4
  
  this%SubCellFlows(1) = subFlow1
  this%SubCellFlows(2) = subFlow2
  this%SubCellFlows(3) = subFlow3
  this%SubCellFlows(4) = subFlow4
  
  end subroutine pr_SetSubCellFlows
  
!------------------------------------------
  function pr_GetSubCellFlow(this, n) result(flow)
  implicit none
  class(ModpathCellDataType) :: this
  doubleprecision :: flow
  integer :: n
  
  flow = this%SubCellFlows(n)
  
  end function pr_GetSubCellFlow

!------------------------------------------

  subroutine pr_ComputeSubCellFlows(this)
  implicit none 
  class(ModpathCellDataType) :: this
  doubleprecision,dimension(4) :: b
  doubleprecision :: rhs1,rhs2,rhs3,qfaces,qsrc,qsink,qsto
  
  rhs1 = 0d0
  rhs2 = 0d0
  rhs3 = 0d0
  
  if(this%GetSubCellCount() .eq. 1) return
  
  ! Initialize b
  b(1) = 0
  b(2) = 0
  b(3) = 0
  b(4) = 0
  
  ! Compute internal source/sink values and set the right hand side
  qsrc = this%SourceFlow / 4d0
  qsink = this%SinkFlow / 4d0
  qsto = this%StorageFlow / 4d0
  
  ! Sub-cell 1
  qfaces = 0d0
  qfaces = qfaces + this%GetSubCellBoundaryFlow(1, 1)
  qfaces = qfaces - this%GetSubCellBoundaryFlow(4, 1)
  qfaces = qfaces + this%GetSubCellBoundaryFlow(5, 1)
  qfaces = qfaces - this%GetSubCellBoundaryFlow(6, 1)
  b(1) = qfaces + qsrc + qsink + qsto
  
  ! Sub-cell 2
  qfaces = 0d0
  qfaces = qfaces - this%GetSubCellBoundaryFlow(2, 1)
  qfaces = qfaces - this%GetSubCellBoundaryFlow(4, 2)
  qfaces = qfaces + this%GetSubCellBoundaryFlow(5, 2)
  qfaces = qfaces - this%GetSubCellBoundaryFlow(6, 2)
  b(2) = qfaces + qsrc + qsink + qsto
  
  ! Sub-cell 3
  qfaces = 0d0
  qfaces = qfaces + this%GetSubCellBoundaryFlow(1, 2)
  qfaces = qfaces + this%GetSubCellBoundaryFlow(3, 1)
  qfaces = qfaces + this%GetSubCellBoundaryFlow(5, 3)
  qfaces = qfaces - this%GetSubCellBoundaryFlow(6, 3)
  b(3) = qfaces + qsrc + qsink + qsto
  
  ! Solve for subcell flows
  this%SubCellFlows(1) = 5d-1*(b(1) + 5d-1*(b(3) - b(2)))
  this%SubCellFlows(3) = -b(1) + this%SubCellFlows(1)
  this%SubCellFlows(2) = b(3) - this%SubCellFlows(3)
  this%SubCellFlows(4) = -b(2) - this%SubCellFlows(1)
  
  end subroutine pr_ComputeSubCellFlows

!------------------------------------------
  function pr_GetSubCellBoundaryFlow(this, faceNumber, subFaceNumber) result(flow)
  implicit none 
  class(ModpathCellDataType) :: this
  integer,intent(in) :: faceNumber,subFaceNumber
  doubleprecision :: flow
  
  select case(faceNumber)
       case (1)
            if(this%SubFaceCounts(1) .eq. 1) then
                 flow = this%Q1(1) / 2.0d0
            else
                 flow = this%Q1(subFaceNumber)
            end if
       case (2)
            if(this%SubFaceCounts(2) .eq. 1) then
                 flow = this%Q2(1) / 2.0d0
            else
                 flow = this%Q2(subFaceNumber)
            end if
       case (3)
            if(this%SubFaceCounts(3) .eq. 1) then
                 flow = this%Q3(1) / 2.0d0
            else
                 flow = this%Q3(subFaceNumber)
            end if
       case (4)
            if(this%SubFaceCounts(4) .eq. 1) then
                 flow = this%Q4(1) / 2.0d0
            else
                 flow = this%Q4(subFaceNumber)
            end if
       case (5)
            if(this%SubFaceCounts(5) .eq. 1) then
                 flow = this%Q5(1) / 4.0d0
            else
                 flow = this%Q5(subFaceNumber)
            end if
       case (6)
            if(this%SubFaceCounts(6) .eq. 1) then
                 flow = this%Q6(1) / 4.0d0
            else
                 flow = this%Q6(subFaceNumber)
            end if
  end select
  
  end function pr_GetSubCellBoundaryFlow

!------------------------------------------
  function pr_GetAveragedFaceFlow(this, faceNumber) result(flow)
  implicit none 
  class(ModpathCellDataType) :: this
  integer,intent(in) :: faceNumber
  doubleprecision :: flow
  
  flow = 0d0
  !select case (faceNumber)
  !  case (1)
  !    arraySize = size(this%Q1)
  !    do n = 1, arraySize
  !      flow = flow + this%Q1(n)
  !    end do
  !  case (2)
  !    arraySize = size(this%Q2)
  !    do n = 1, arraySize
  !      flow = flow + this%Q2(n)
  !    end do
  !  case (3)
  !    arraySize = size(this%Q3)
  !    do n = 1, arraySize
  !      flow = flow + this%Q3(n)
  !    end do
  !  case (4)
  !    arraySize = size(this%Q4)
  !    do n = 1, arraySize
  !      flow = flow + this%Q4(n)
  !    end do
  !  case (5)
  !    arraySize = size(this%Q5)
  !    do n = 1, arraySize
  !      flow = flow + this%Q5(n)
  !    end do
  !  case (6)
  !    arraySize = size(this%Q6)
  !    do n = 1, arraySize
  !      flow = flow + this%Q6(n)
  !    end do
  !  case default
  !    ! do nothing
  !end select

  ! RWPT: faster
  select case (faceNumber)
    case (1)
      flow = sum( this%Q1 )
    case (2)
      flow = sum( this%Q2 )
    case (3)
      flow = sum( this%Q3 )
    case (4)
      flow = sum( this%Q4 )
    case (5)
      flow = sum( this%Q5 )
    case (6)
      flow = sum( this%Q6 )
    case default
      ! do nothing
  end select
  

  end function pr_GetAveragedFaceFlow

!------------------------------------------
  subroutine pr_AssignAveragedFaceFlowArray(this, flows)
  implicit none 
  class(ModpathCellDataType) :: this
  doubleprecision,dimension(6) :: flows
  integer :: n
  
  do n = 1, 6
    flows(n) = this%GetAveragedFaceFlow(n)
  end do
  
  end subroutine pr_AssignAveragedFaceFlowArray


!------------------------------------------
  subroutine pr_FillMassSubCellDataBuffer(this, subCellData, subRow, subColumn, backwardTracking)
  implicit none 
  class(ModpathCellDataType) :: this
  type(ModpathSubCellDataType),intent(inout) :: subCellData
  integer,intent(in) :: subRow,subColumn
  logical,intent(in) :: backwardTracking
  doubleprecision,dimension(6) :: flows
  doubleprecision :: sign,xinc,yinc
  integer :: n,rowcolumn,count
  
  ! Reset the data in subCellData (Not strictly necessary, but useful for debugging purposes)
  call subCellData%Reset()
  
  call this%FillSubCellFaceFlowsBuffer(subRow, subColumn, flows)
  
  subCellData%DX = this%DX / dble(this%SubCellColumnCount)
  subCellData%DY = this%DY / dble(this%SubCellRowCount)
  ! For RWPT and drying cells generate blow up due to
  ! default value 1e-4 for dry cells.
  !subCellData%DZ = this%GetDZ()
  ! This form returns cell vertical size, 
  ! given by the grid, regardless of 
  ! saturation status. The latter is stored 
  ! as a cell variable for later use.
  subCellData%DZ = this%GetDZRW()

  ! RWPT
  ! dry and partiallyDry properties where defined when calling GetDZ()
  subCellData%dry          = this%dry
  subCellData%partiallyDry = this%partiallyDry
  subCellData%Head         = this%Head
  subCellData%Top          = this%Top
  subCellData%Bottom       = this%Bottom
  ! END RWPT
  
  sign = 1.0d0
  if(backwardTracking) sign = -sign
  subCellData%VX1 = sign * flows(1) / subCellData%DY /subCellData%DZ / this%Porosity / this%Retardation
  subCellData%VX2 = sign * flows(2) / subCellData%DY /subCellData%DZ / this%Porosity / this%Retardation
  subCellData%VY1 = sign * flows(3) / subCellData%DX /subCellData%DZ / this%Porosity / this%Retardation
  subCellData%VY2 = sign * flows(4) / subCellData%DX /subCellData%DZ / this%Porosity / this%Retardation
  subCellData%VZ1 = sign * flows(5) / subCellData%DX /subCellData%DY / this%Porosity / this%Retardation
  subCellData%VZ2 = sign * flows(6) / subCellData%DX /subCellData%DY / this%Porosity / this%Retardation
  
  subCellData%Row = subRow
  subCellData%Column = subColumn
  
  xinc = 1.0d0 / dble(this%SubCellColumnCount)
  subCellData%OffsetX(1) = (subColumn - 1) * xinc
  subCellData%OffsetX(2) = 1.0d0
  if(subColumn .lt. this%SubCellColumnCount) then
    subCellData%OffsetX(2) = subColumn * xinc
  end if
  
  yinc = 1.0d0 / dble(this%SubCellRowCount)
  subCellData%OffsetY(1) = (this%SubCellRowCount - subRow) * yinc
  subCellData%OffsetY(2) = 1.0d0
  if(subRow .gt. 1.0d0) then
    subCellData%OffsetY(2) = (this%SubCellRowCount - subRow + 1) * yinc
  end if
  
  subCellData%OffsetZ(1) = 0d0
  subCellData%OffsetZ(2) = 1.0d0
  
  ! Assign the connections for the 6 faces.
  ! All internal connections are set to -1.
  ! Boundary cells are set to the node number of the neighbor cell.
  ! Boundary faces that do not have adjacent neighbors are set to 0.
  
  do n = 1, 6
    ! Start by initializing all face connections to -1. 
    subCellData%Connection(n) = -1

    ! Start by initializing mass boundaries to zero
    subCellData%MassBoundary(n) = 0

  end do
  
  ! Assign the actual connection values to all of the faces that are not internal faces.
  if(subCellData%Row .eq. 1) then
    count = this%GetSubFaceCount(4)
    subCellData%Connection(4) = 0
    if(count .eq. 1) then
      subCellData%Connection(4) = this%SubFaceConn4(1)
      subCellData%MassBoundary(4) = this%MassBoundarySubFace4(1)
    else if(count .gt. 1) then
      subCellData%Connection(4) = this%SubFaceConn4(subCellData%Column)
      subCellData%MassBoundary(4) = this%MassBoundarySubFace4(subCellData%Column)
    end if
  end if
  
  if(subCellData%Row .eq. this%SubCellRowCount) then
    count = this%GetSubFaceCount(3)
    subCellData%Connection(3) = 0
    if(count .eq. 1) then
      subCellData%Connection(3) = this%SubFaceConn3(1)
      subCellData%MassBoundary(3) = this%MassBoundarySubFace3(1)
    else if(count .gt. 1) then
      subCellData%Connection(3) = this%SubFaceConn3(subCellData%Column)
      subCellData%MassBoundary(3) = this%MassBoundarySubFace3(subCellData%Column)
    end if
  end if
  
  if(subCellData%Column .eq. 1) then
    count = this%GetSubFaceCount(1)
    subCellData%Connection(1) = 0
    if(count .eq. 1) then
      subCellData%Connection(1) = this%SubFaceConn1(1)
      subCellData%MassBoundary(1) = this%MassBoundarySubFace1(1)
    else if(count .gt. 1) then
      subCellData%Connection(1) = this%SubFaceConn1(subCellData%Row)
      subCellData%MassBoundary(1) = this%MassBoundarySubFace1(subCellData%Row)
    end if
  end if
  
  if(subCellData%Column .eq. this%SubCellColumnCount) then
    count = this%GetSubFaceCount(2)
    subCellData%Connection(2) = 0
    if(count .eq. 1) then
      subCellData%Connection(2) = this%SubFaceConn2(1)
      subCellData%MassBoundary(2) = this%MassBoundarySubFace2(1)
    else if(count .gt. 1) then
      subCellData%Connection(2) = this%SubFaceConn2(subCellData%Row)
      subCellData%MassBoundary(2) = this%MassBoundarySubFace2(subCellData%Row)
    end if
  end if
  
  rowcolumn = this%SubCellRowCount * this%SubCellColumnCount
  if(this%GetSubFaceCount(5) .eq. 1) then
    subCellData%Connection(5) = this%SubFaceConn5(1)
    subCellData%MassBoundary(5) = this%MassBoundarySubFace5(1)
  else if(this%GetSubFaceCount(5) .eq. rowcolumn) then
    n = ((subCellData%Row - 1) * this%SubCellColumnCount) + subCellData%Column
    subCellData%Connection(5) = this%SubFaceConn5(n)
    subCellData%MassBoundary(5) = this%MassBoundarySubFace5(n)
  end if
  
  if(this%GetSubFaceCount(6) .eq. 1) then
    subCellData%Connection(6) = this%SubFaceConn6(1)
    subCellData%MassBoundary(6) = this%MassBoundarySubFace6(1)
  else if(this%GetSubFaceCount(6) .eq. rowcolumn) then
    n = ((subCellData%Row - 1) * this%SubCellColumnCount) + subCellData%Column
    subCellData%Connection(6) = this%SubFaceConn6(n)
    subCellData%MassBoundary(6) = this%MassBoundarySubFace6(n)
  end if


  ! RWPT
  ! Assign additional properties to sub cell buffer
  ! need it ?
  subCellData%Porosity    = this%Porosity 
  subCellData%Retardation = this%Retardation 


  ! Necessary for distributed dispersivities
  subCellData%alphaLH = this%alphaLH
  subCellData%alphaLV = this%alphaLV
  subCellData%alphaTH = this%alphaTH
  subCellData%alphaTV = this%alphaTV
  subCellData%dMEff   = this%dMEff


  end subroutine pr_FillMassSubCellDataBuffer


!------------------------------------------
  subroutine pr_FillSubCellDataBuffer(this, subCellData, subRow, subColumn, backwardTracking)
  implicit none 
  class(ModpathCellDataType) :: this
  type(ModpathSubCellDataType),intent(inout) :: subCellData
  integer,intent(in) :: subRow,subColumn
  logical,intent(in) :: backwardTracking
  doubleprecision,dimension(6) :: flows
  doubleprecision :: sign,xinc,yinc
  integer :: n,rowcolumn,count
  
  ! Reset the data in subCellData (Not strictly necessary, but useful for debugging purposes)
  call subCellData%Reset()
  
  call this%FillSubCellFaceFlowsBuffer(subRow, subColumn, flows)
  
  subCellData%DX = this%DX / dble(this%SubCellColumnCount)
  subCellData%DY = this%DY / dble(this%SubCellRowCount)
  subCellData%DZ = this%GetDZ()
  
  sign = 1.0d0
  if(backwardTracking) sign = -sign
  subCellData%VX1 = sign * flows(1) / subCellData%DY /subCellData%DZ / this%Porosity / this%Retardation
  subCellData%VX2 = sign * flows(2) / subCellData%DY /subCellData%DZ / this%Porosity / this%Retardation
  subCellData%VY1 = sign * flows(3) / subCellData%DX /subCellData%DZ / this%Porosity / this%Retardation
  subCellData%VY2 = sign * flows(4) / subCellData%DX /subCellData%DZ / this%Porosity / this%Retardation
  subCellData%VZ1 = sign * flows(5) / subCellData%DX /subCellData%DY / this%Porosity / this%Retardation
  subCellData%VZ2 = sign * flows(6) / subCellData%DX /subCellData%DY / this%Porosity / this%Retardation
  
  subCellData%Row = subRow
  subCellData%Column = subColumn
  
  xinc = 1.0d0 / dble(this%SubCellColumnCount)
  subCellData%OffsetX(1) = (subColumn - 1) * xinc
  subCellData%OffsetX(2) = 1.0d0
  if(subColumn .lt. this%SubCellColumnCount) then
    subCellData%OffsetX(2) = subColumn * xinc
  end if
  
  yinc = 1.0d0 / dble(this%SubCellRowCount)
  subCellData%OffsetY(1) = (this%SubCellRowCount - subRow) * yinc
  subCellData%OffsetY(2) = 1.0d0
  if(subRow .gt. 1.0d0) then
    subCellData%OffsetY(2) = (this%SubCellRowCount - subRow + 1) * yinc
  end if
  
  subCellData%OffsetZ(1) = 0d0
  subCellData%OffsetZ(2) = 1.0d0
  
  ! Assign the connections for the 6 faces.
  ! All internal connections are set to -1.
  ! Boundary cells are set to the node number of the neighbor cell.
  ! Boundary faces that do not have adjacent neighbors are set to 0.
  
  ! Start by initializing all face connections to -1. 
  do n = 1, 6
    subCellData%Connection(n) = -1
  end do
  
  ! Assign the actual connection values to all of the faces that are not internal faces.
  if(subCellData%Row .eq. 1) then
    count = this%GetSubFaceCount(4)
    subCellData%Connection(4) = 0
    if(count .eq. 1) then
      subCellData%Connection(4) = this%SubFaceConn4(1)
    else if(count .gt. 1) then
      subCellData%Connection(4) = this%SubFaceConn4(subCellData%Column)
    end if
  end if
  
  if(subCellData%Row .eq. this%SubCellRowCount) then
    count = this%GetSubFaceCount(3)
    subCellData%Connection(3) = 0
    if(count .eq. 1) then
      subCellData%Connection(3) = this%SubFaceConn3(1)
    else if(count .gt. 1) then
      subCellData%Connection(3) = this%SubFaceConn3(subCellData%Column)
    end if
  end if
  
  if(subCellData%Column .eq. 1) then
    count = this%GetSubFaceCount(1)
    subCellData%Connection(1) = 0
    if(count .eq. 1) then
      subCellData%Connection(1) = this%SubFaceConn1(1)
    else if(count .gt. 1) then
      subCellData%Connection(1) = this%SubFaceConn1(subCellData%Row)
    end if
  end if
  
  if(subCellData%Column .eq. this%SubCellColumnCount) then
    count = this%GetSubFaceCount(2)
    subCellData%Connection(2) = 0
    if(count .eq. 1) then
      subCellData%Connection(2) = this%SubFaceConn2(1)
    else if(count .gt. 1) then
      subCellData%Connection(2) = this%SubFaceConn2(subCellData%Row)
    end if
  end if
  
  rowcolumn = this%SubCellRowCount * this%SubCellColumnCount
  if(this%GetSubFaceCount(5) .eq. 1) then
    subCellData%Connection(5) = this%SubFaceConn5(1)
  else if(this%GetSubFaceCount(5) .eq. rowcolumn) then
    n = ((subCellData%Row - 1) * this%SubCellColumnCount) + subCellData%Column
    subCellData%Connection(5) = this%SubFaceConn5(n)
  end if
  
  if(this%GetSubFaceCount(6) .eq. 1) then
    subCellData%Connection(6) = this%SubFaceConn6(1)
  else if(this%GetSubFaceCount(6) .eq. rowcolumn) then
    n = ((subCellData%Row - 1) * this%SubCellColumnCount) + subCellData%Column
    subCellData%Connection(6) = this%SubFaceConn6(n)
  end if
  
  end subroutine pr_FillSubCellDataBuffer

!------------------------------------------
  function pr_GetSubCellData(this, subRow, subColumn, backwardTracking) result(subCellData)
  implicit none 
  class(ModpathCellDataType) :: this
  integer,intent(in) :: subRow,subColumn
  logical,intent(in) :: backwardTracking
  type(ModpathSubCellDataType) :: subCellData
  
  call this%FillSubCellDataBuffer(subCellData, subRow, subColumn, backwardTracking)

  end function pr_GetSubCellData
  
!------------------------------------------
  subroutine pr_FillSubCellFaceFlowsBuffer(this, subRow, subColumn, faceFlows)
  implicit none 
  class(ModpathCellDataType) :: this
  integer,intent(in) :: subRow,subColumn
  doubleprecision,dimension(6) :: faceFlows
  integer :: subCellNumber
  
  if(this%GetSubCellCount() .eq. 1) then
      faceFlows(1) = this%Q1(1)
      faceFlows(2) = this%Q2(1)
      faceFlows(3) = this%Q3(1)
      faceFlows(4) = this%Q4(1)
      faceFlows(5) = this%Q5(1)
      faceFlows(6) = this%Q6(1)
  else
      subCellNumber = (subRow - 1) * this%GetSubCellColumnCount() + subColumn
      select case (subCellNumber)
        case (1)
          faceFlows(1) = this%GetSubCellBoundaryFlow(1, 1)
          faceFlows(2) = this%SubCellFlows(1)
          faceFlows(3) = this%SubCellFlows(3)
          faceFlows(4) = this%GetSubCellBoundaryFlow(4, 1)
          faceFlows(5) = this%GetSubCellBoundaryFlow(5, 1)
          faceFlows(6) = this%GetSubCellBoundaryFlow(6, 1)
        case (2)
          faceFlows(1) = this%SubCellFlows(1)
          faceFlows(2) = this%GetSubCellBoundaryFlow(2, 1)
          faceFlows(3) = this%SubCellFlows(4)
          faceFlows(4) = this%GetSubCellBoundaryFlow(4, 2)
          faceFlows(5) = this%GetSubCellBoundaryFlow(5, 2)
          faceFlows(6) = this%GetSubCellBoundaryFlow(6, 2)
        case (3)
          faceFlows(1) = this%GetSubCellBoundaryFlow(1, 2)
          faceFlows(2) = this%SubCellFlows(2)
          faceFlows(3) = this%GetSubCellBoundaryFlow(3, 1)
          faceFlows(4) = this%SubCellFlows(3)
          faceFlows(5) = this%GetSubCellBoundaryFlow(5, 3)
          faceFlows(6) = this%GetSubCellBoundaryFlow(6, 3)
        case (4)
          faceFlows(1) = this%SubCellFlows(2)
          faceFlows(2) = this%GetSubCellBoundaryFlow(2, 2)
          faceFlows(3) = this%GetSubCellBoundaryFlow(3, 2)
          faceFlows(4) = this%SubCellFlows(4)
          faceFlows(5) = this%GetSubCellBoundaryFlow(5, 4)
          faceFlows(6) = this%GetSubCellBoundaryFlow(6, 4)
        case default
          ! for now, do nothing.
      end select
  end if
  
  end subroutine pr_FillSubCellFaceFlowsBuffer


  !-----------------------------------------------
  ! RWPT METHODS
  !-----------------------------------------------

  ! RWPT-USG
  function pr_GetNeighborSubCellIndexes( this, subRow, subColumn ) result( neighborSubCellIndexes )
      !-----------------------------------------------------------
      ! Reference cell locations in 3D neighbor cell buffer array
      ! 1 : dc1
      ! 2 : ic13
      ! 3 : ic14
      ! 4 : dc2
      ! 5 : ic23
      ! 6 : ic24
      ! 7 : dc3
      ! 8 : ic35
      ! 9 : ic36
      ! 10: dc4
      ! 11: ic45
      ! 12: ic46
      ! 13: dc5
      ! 14: ic51
      ! 15: ic52
      ! 16: dc6
      ! 17: ic61
      ! 18: ic62
      ! 
      ! In neighborSubCellIndexes: bufferIndex, subRow, subColumn
      ! 
      !---------------------------------------------------------
      implicit none
      class( ModpathCellDataType ) :: this
      integer, intent(in) :: subRow, subColumn
      integer, dimension(3,18) :: neighborSubCellIndexes
      !---------------------------------------------------------

      ! If current cell is not refined
      if( this%GetSubCellCount() .eq. 1 ) then 

          neighborSubCellIndexes(:,1)  = [ 1,1,1]
          neighborSubCellIndexes(:,2)  = [ 2,1,1]
          neighborSubCellIndexes(:,3)  = [ 3,1,1]
          neighborSubCellIndexes(:,4)  = [ 4,1,1]
          neighborSubCellIndexes(:,5)  = [ 5,1,1]
          neighborSubCellIndexes(:,6)  = [ 6,1,1]
          neighborSubCellIndexes(:,7)  = [ 7,1,1]
          neighborSubCellIndexes(:,8)  = [ 8,1,1]
          neighborSubCellIndexes(:,9)  = [ 9,1,1]
          neighborSubCellIndexes(:,10) = [10,1,1]
          neighborSubCellIndexes(:,11) = [11,1,1]
          neighborSubCellIndexes(:,12) = [12,1,1]
          neighborSubCellIndexes(:,13) = [13,1,1]
          neighborSubCellIndexes(:,14) = [14,1,1]
          neighborSubCellIndexes(:,15) = [15,1,1]
          neighborSubCellIndexes(:,16) = [16,1,1]
          neighborSubCellIndexes(:,17) = [17,1,1]
          neighborSubCellIndexes(:,18) = [18,1,1]

          ! Done
          return

      end if

      ! If current cell is refined then 
      ! subRow and subColumn determine indexes from neighbor cells/subcells 
      ! to be used in corner interpolation process 

      if ( subRow .eq. 1 ) then
          ! If subcell from first subRow 
          if ( subColumn .eq. 1 ) then
              ! If subcell from first subColumn
              neighborSubCellIndexes(:,1)  = [ 1,1,2]
              neighborSubCellIndexes(:,2)  = [ 1,2,2]
              neighborSubCellIndexes(:,3)  = [ 3,2,2]
              neighborSubCellIndexes(:,4)  = [ 0,1,2]
              neighborSubCellIndexes(:,5)  = [ 0,2,2]
              neighborSubCellIndexes(:,6)  = [10,2,2]
              neighborSubCellIndexes(:,7)  = [ 0,2,1]
              neighborSubCellIndexes(:,8)  = [13,2,1]
              neighborSubCellIndexes(:,9)  = [16,2,1]
              neighborSubCellIndexes(:,10) = [10,2,1]
              neighborSubCellIndexes(:,11) = [11,2,1]
              neighborSubCellIndexes(:,12) = [12,2,1]
              neighborSubCellIndexes(:,13) = [13,1,1]
              neighborSubCellIndexes(:,14) = [14,1,2]
              neighborSubCellIndexes(:,15) = [13,1,2]
              neighborSubCellIndexes(:,16) = [16,1,1]
              neighborSubCellIndexes(:,17) = [17,1,2]
              neighborSubCellIndexes(:,18) = [16,1,2]
          else
              ! If subcell from second subColumn
              neighborSubCellIndexes(:,1)  = [ 0,1,1]
              neighborSubCellIndexes(:,2)  = [ 0,2,1]
              neighborSubCellIndexes(:,3)  = [10,2,1]
              neighborSubCellIndexes(:,4)  = [ 4,1,1]
              neighborSubCellIndexes(:,5)  = [ 4,2,1]
              neighborSubCellIndexes(:,6)  = [ 6,2,1]
              neighborSubCellIndexes(:,7)  = [ 0,2,2]
              neighborSubCellIndexes(:,8)  = [13,2,2]
              neighborSubCellIndexes(:,9)  = [16,2,2]
              neighborSubCellIndexes(:,10) = [10,2,2]
              neighborSubCellIndexes(:,11) = [11,2,2]
              neighborSubCellIndexes(:,12) = [12,2,2]
              neighborSubCellIndexes(:,13) = [13,1,2]
              neighborSubCellIndexes(:,14) = [13,1,1]
              neighborSubCellIndexes(:,15) = [15,1,1]
              neighborSubCellIndexes(:,16) = [16,1,2]
              neighborSubCellIndexes(:,17) = [16,1,1]
              neighborSubCellIndexes(:,18) = [18,1,1]
          end if        
      else
          ! If subcell from second subRow 
          if ( subColumn .eq. 1 ) then 
              ! If subcell from first subColumn
              neighborSubCellIndexes(:,1)  = [ 1,2,2]
              neighborSubCellIndexes(:,2)  = [ 2,1,2]
              neighborSubCellIndexes(:,3)  = [ 1,1,2]
              neighborSubCellIndexes(:,4)  = [ 0,2,2]
              neighborSubCellIndexes(:,5)  = [ 7,1,2]
              neighborSubCellIndexes(:,6)  = [ 0,1,2]
              neighborSubCellIndexes(:,7)  = [ 7,1,1]
              neighborSubCellIndexes(:,8)  = [ 8,1,1]
              neighborSubCellIndexes(:,9)  = [ 9,1,1]
              neighborSubCellIndexes(:,10) = [ 0,1,1]
              neighborSubCellIndexes(:,11) = [13,1,1]
              neighborSubCellIndexes(:,12) = [16,1,1]
              neighborSubCellIndexes(:,13) = [13,2,1]
              neighborSubCellIndexes(:,14) = [14,2,2]
              neighborSubCellIndexes(:,15) = [13,2,2]
              neighborSubCellIndexes(:,16) = [16,2,1]
              neighborSubCellIndexes(:,17) = [17,2,2]
              neighborSubCellIndexes(:,18) = [16,2,2]
          else
              ! If subcell from second subColumn
              neighborSubCellIndexes(:,1)  = [ 0,2,1]
              neighborSubCellIndexes(:,2)  = [ 7,1,1]
              neighborSubCellIndexes(:,3)  = [ 0,1,1]
              neighborSubCellIndexes(:,4)  = [ 4,2,1]
              neighborSubCellIndexes(:,5)  = [ 5,1,1]
              neighborSubCellIndexes(:,6)  = [ 4,1,1]
              neighborSubCellIndexes(:,7)  = [ 7,1,2]
              neighborSubCellIndexes(:,8)  = [ 8,1,2]
              neighborSubCellIndexes(:,9)  = [ 9,1,2]
              neighborSubCellIndexes(:,10) = [ 0,1,2]
              neighborSubCellIndexes(:,11) = [13,1,2]
              neighborSubCellIndexes(:,12) = [16,1,2]
              neighborSubCellIndexes(:,13) = [13,2,2]
              neighborSubCellIndexes(:,14) = [13,2,1]
              neighborSubCellIndexes(:,15) = [15,2,1]
              neighborSubCellIndexes(:,16) = [16,2,2]
              neighborSubCellIndexes(:,17) = [16,2,1]
              neighborSubCellIndexes(:,18) = [18,2,1]
          end if        
      end if 


      ! Done
      return


  end function pr_GetNeighborSubCellIndexes 



  ! RWPT-USG
  subroutine pr_FillSubCellFaceAreas(this, areas, skipSubCells)
  !---------------------------------------------------------
  !
  !---------------------------------------------------------
  implicit none
  class( ModpathCellDataType ) :: this
  doubleprecision, dimension(3), intent(inout) :: areas
  doubleprecision :: dx, dy, dz
  logical, intent(in) :: skipSubCells
  !---------------------------------------------------------

    if ( .not. skipSubCells ) then
      dx = this%DX / dble(this%SubCellColumnCount)
      dy = this%DY / dble(this%SubCellRowCount)
    else
      dx = this%DX
      dy = this%DY
    end if
    ! areas are used for flow purposes so 
    ! this function is not a problem in RWPT.
    ! Verify that it does not defines cell properties...

    ! If dry cell, this value is really small 1e-4
    !dz = this%GetDZ() 
    dz = this%GetDZRW() 

    areas(1) = dy*dz ! These areas are intended for balances in corners 
    areas(2) = dx*dz
    areas(3) = dx*dy
      
    return


  end subroutine pr_FillSubCellFaceAreas 


  ! RWPT-USG
  subroutine pr_FillSubCellFaceFlows(this, subRow, subColumn, faceFlows, skipSubCells)
  !---------------------------------------------------------
  !
  !---------------------------------------------------------
  implicit none 
  class(ModpathCellDataType) :: this
  integer,intent(in) :: subRow,subColumn
  doubleprecision,dimension(6) :: faceFlows
  logical, intent(in) :: skipSubCells
  !---------------------------------------------------------
 
    
    if ( .not. skipSubCells ) then
      ! Get face flows as usual   
      call this%FillSubCellFaceFlowsBuffer( subRow, subColumn, faceFlows )
    else
      ! Skip subcell indexation and bring cell face flows
      faceFlows(1) = this%GetAveragedFaceFlow(1)
      faceFlows(2) = this%GetAveragedFaceFlow(2)
      faceFlows(3) = this%GetAveragedFaceFlow(3)
      faceFlows(4) = this%GetAveragedFaceFlow(4)
      faceFlows(5) = this%GetAveragedFaceFlow(5)
      faceFlows(6) = this%GetAveragedFaceFlow(6)
    end if

  
  end subroutine pr_FillSubCellFaceFlows


  ! RWPT-USG
  function pr_GetVolume( this, skipSubCells ) result( volume )
  !---------------------------------------------------------
  !
  !---------------------------------------------------------
  implicit none
  class( ModpathCellDataType ) :: this
  ! input
  logical, intent(in) :: skipSubCells
  ! output
  doubleprecision :: volume
  ! local
  doubleprecision :: dx, dy, dz
  !---------------------------------------------------------

    if ( .not. skipSubCells ) then
      dx = this%DX / dble(this%SubCellColumnCount)
      dy = this%DY / dble(this%SubCellRowCount)
    else
      dx = this%DX
      dy = this%DY
    end if
    !dz = this%GetDZ()
    dz = this%GetDZRW()

    volume = dx*dy*dz
      
    return

  end function pr_GetVolume 


  function pr_GetWaterVolume( this ) result( waterVolume )
  !---------------------------------------------------------
  !
  !---------------------------------------------------------
  implicit none
  class( ModpathCellDataType ) :: this
  ! output
  doubleprecision :: waterVolume
  ! local
  doubleprecision :: dx, dy, dz
  !---------------------------------------------------------

    waterVolume = 0d0

    dx = this%DX
    dy = this%DY
    dz = this%GetDZRW() ! brings the saturation
    
    ! If dry no water
    if ( this%dry ) return

    waterVolume = this%Porosity*dx*dy*dz
      
    return

  end function pr_GetWaterVolume 


  subroutine pr_VerifyDryCell(this)
  !---------------------------------------------------------
  !
  !---------------------------------------------------------
  implicit none
  class(ModpathCellDataType) :: this
  doubleprecision :: dz
  doubleprecision, dimension(6) :: flows
  !---------------------------------------------------------
  
  dz = this%Top - this%Bottom
  this%dry = .false.
  this%partiallyDry = .false.

  ! If the layer is convertible, set dz = Head - Bottom if Head < Top
  if(this%LayerType .eq. 1) then
    if(this%Head .lt. this%Top) then 
      ! if dz < 0, means that the cell is completely dry
      ! if dz > 0, but Head < Top, cell is partially dry
      ! this second case can be used for displacing particles by RWPT in 
      ! that are within the region that is partially saturated.
      dz = this%Head - this%Bottom
      this%dry = .false.
      this%partiallyDry = .true.
      ! If dz < 0, set dz to an arbitrary, small positive value (MODPATH default)
      if(dz .le. 0.0d0) then 
        dz = 1.0d-4
        this%dry = .true.
        this%partiallyDry = .false.

        ! The alternative to handle MF6 models with Newton-Raphson 
        ! If a cell has head<bot, but there are flows in it, 
        ! mark only as partiallyDry in order to allow the 
        ! the displacement of particles following the 
        ! cell flows.
        call this%AssignAveragedFaceFlowArray( flows )
        if ( any(flows.ne.0d0) ) then 
          this%dry = .false.
          this%partiallyDry = .true.
        end if 
      end if 
    end if 
  end if
 
  end subroutine pr_VerifyDryCell


  function pr_GetDZRW(this) result(dz)
  !-----------------------------------------------
  ! Return dz as domain size when dry
  !-----------------------------------------------
  implicit none
  class(ModpathCellDataType) :: this
  doubleprecision :: dz, dzc
  doubleprecision, dimension(6) :: flows
  !-----------------------------------------------

  dz  = this%Top - this%Bottom
  dzc = this%Top - this%Bottom
  this%dry = .false.
  this%partiallyDry = .false.

  ! If the layer is convertible, set dz = Head - Bottom if Head < Top
  if(this%LayerType .eq. 1) then
    if(this%Head .lt. this%Top) then 
      ! if dz < 0, means that the cell is completely dry
      ! if dz > 0, but Head < Top, cell is partially dry
      ! this second case can be used for displacing particles by RWPT
      ! that are within the region that is partially saturated.
      dz = this%Head - this%Bottom
      this%dry = .false.
      this%partiallyDry = .true.
      ! If dz < 0, set dz to an arbitrary, small positive value (MODPATH default)
      if ( &
        ( dz .le. 0d0 ) .or. &
        ( (this%Head - this%Bottom)/dzc .lt. 0.01d0) ) then 
        dz = this%Top - this%Bottom ! to avoid blow up dzrw/dz
        this%dry = .true.
        this%partiallyDry = .false.

        ! The alternative to handle MF6 models with Newton-Raphson 
        ! If a cell has head<bot, but there are flows in it, 
        ! mark only as partiallyDry in order to allow the 
        ! the displacement of particles following the 
        ! cell flows.
        call this%AssignAveragedFaceFlowArray( flows )
        if ( any(flows.ne.0d0) ) then 
          this%dry = .false.
          this%partiallyDry = .true.
        end if 
      end if 
    end if
  end if
 

  end function pr_GetDZRW




end module ModpathCellDataModule
