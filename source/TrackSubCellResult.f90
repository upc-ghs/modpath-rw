module TrackSubCellResultModule
  use ParticleLocationModule,only : ParticleLocationType
  use ModpathSubCellDataModule,only : ModpathSubCellDataType
  implicit none
  
! Set default access status to private
  private

  type,public :: TrackSubCellResultType 
    integer :: CellNumber,Row,Column
    type(ParticleLocationType) :: InitialLocation
    type(ParticleLocationType) :: FinalLocation
    integer :: ExitFace = 0
    integer :: ExitFaceConnection = 0
!    logical :: InternalExitFace = .false.
    doubleprecision :: MaximumTime = 0d0
    integer :: Status = 0
  contains
    procedure :: Reset=>pr_Reset
    procedure :: Status_Undefined=>TrackSubCellResultType_Undefined
    procedure :: Status_ExitAtCellFace=>TrackSubCellResultType_ExitAtCellFace
    procedure :: Status_ExitAtInternalFace=>TrackSubCellResultType_ExitAtInternalFace
    procedure :: Status_ReachedMaximumTime=>TrackSubCellResultType_ReachedMaximumTime
    procedure :: Status_NoExitPossible=>TrackSubCellResultType_NoExitPossible
    procedure :: InitializeLocation=>pr_InitializeLocation
  end type
  
contains

subroutine pr_Reset(this)
  implicit none
  class(TrackSubCellResultType) :: this

  this%CellNumber = 0
  this%ExitFace = 0
  this%ExitFaceConnection = 0
  this%MaximumTime = 0.0d0
  this%Status = 0
  call this%InitialLocation%Reset()
  call this%FinalLocation%Reset()

end subroutine

!---------------------------------------------------
  function TrackSubCellResultType_Undefined(this) result(status)
  implicit none
  class(TrackSubCellResultType) :: this
  integer :: status
  
  status = 0
  
  end function TrackSubCellResultType_Undefined

!---------------------------------------------------
  function TrackSubCellResultType_ExitAtCellFace(this) result(status)
  implicit none
  class(TrackSubCellResultType) :: this
  integer :: status
  
  status = 1
  
  end function TrackSubCellResultType_ExitAtCellFace

!---------------------------------------------------
  function TrackSubCellResultType_ExitAtInternalFace(this) result(status)
  implicit none
  class(TrackSubCellResultType) :: this
  integer :: status
  
  status = 4
  
  end function TrackSubCellResultType_ExitAtInternalFace

!---------------------------------------------------
  function TrackSubCellResultType_ReachedMaximumTime(this) result(status)
  implicit none
  class(TrackSubCellResultType) :: this
  integer :: status
  
  status = 2
  
  end function TrackSubCellResultType_ReachedMaximumTime

!---------------------------------------------------
  function TrackSubCellResultType_NoExitPossible(this) result(status)
  implicit none
  class(TrackSubCellResultType) :: this
  integer :: status
  
  status = 3
  
  end function TrackSubCellResultType_NoExitPossible


  subroutine pr_InitializeLocation( this, initialLocation )
  !---------------------------------------------------
  !
  !---------------------------------------------------
    class(TrackSubCellResultType) :: this
    type(ParticleLocationType),intent(in) :: initialLocation
    integer :: cellNumber
    doubleprecision :: initialX,initialY,initialZ,initialTime
  !---------------------------------------------------
    
    call this%Reset()
    
    cellNumber  = initialLocation%CellNumber
    initialX    = initialLocation%LocalX
    initialY    = initialLocation%LocalY
    initialZ    = initialLocation%LocalZ
    initialTime = initialLocation%TrackingTime
    
    this%CellNumber                   = cellNumber
    this%InitialLocation%CellNumber   = cellNumber
    this%InitialLocation%LocalX       = initialX
    this%InitialLocation%LocalY       = initialY
    this%InitialLocation%LocalZ       = initialZ
    this%InitialLocation%TrackingTime = initialTime
    this%FinalLocation%LocalX         = initialX
    this%FinalLocation%LocalY         = initialY
    this%FinalLocation%LocalZ         = initialZ
    this%FinalLocation%TrackingTime   = initialTime
    this%Status = this%Status_Undefined()

  end subroutine pr_InitializeLocation


end module TrackSubCellResultModule
