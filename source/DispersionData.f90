module DispersionDataModule

  implicit none
  
  ! Set default access status to private
  private

  type, public :: DispersionDataType
    integer            :: id
    character(len=300) :: stringid
    ! 1: linear, 2: nonlinear
    integer            :: modelKind = 1

    ! Dispersivities for linear dispersion model
    doubleprecision,dimension(:),allocatable :: AlphaL
    doubleprecision,dimension(:),allocatable :: AlphaTH
    doubleprecision,dimension(:),allocatable :: AlphaTV
    ! Effective molecular diffusion, corrected by tortuosity,
    ! in case is spatially distributed 
    doubleprecision,dimension(:),allocatable :: DMEff
    ! Ideally there is something to identify whether
    ! effective molecular diffusion is distributed or not

    ! Parameters for nonlinear dispersion model 
    doubleprecision :: dmaqueous
    doubleprecision,dimension(:),allocatable :: BetaL
    doubleprecision,dimension(:),allocatable :: BetaTH
    doubleprecision,dimension(:),allocatable :: BetaTV
    doubleprecision,dimension(:),allocatable :: Delta
    doubleprecision,dimension(:),allocatable :: DGrain

    logical :: initialized =.false.

  contains
    procedure :: Reset => pr_Reset
    procedure :: InitializeByModelKind => pr_InitializeByModelKind
  end type


contains


  subroutine pr_Reset( this )
    !-------------------------------------------------------------------------
    ! Specifications
    !-------------------------------------------------------------------------
    implicit none
    class(DispersionDataType) :: this
    !-------------------------------------------------------------------------

    this%initialized = .false.
    this%id          = 0
    this%stringid    = '' 
    this%dmaqueous   = 0d0
    if(allocated(this%AlphaL)) deallocate(this%AlphaL)
    if(allocated(this%AlphaTH)) deallocate(this%AlphaTH)
    if(allocated(this%AlphaTV)) deallocate(this%AlphaTV)
    if(allocated(this%DMEff)) deallocate(this%DMEff)
    if(allocated(this%BetaL)) deallocate(this%BetaL)
    if(allocated(this%BetaTH)) deallocate(this%BetaTH)
    if(allocated(this%BetaTV)) deallocate(this%BetaTV)
    if(allocated(this%Delta)) deallocate(this%Delta)
    if(allocated(this%DGrain)) deallocate(this%DGrain)

  end subroutine pr_Reset


  subroutine pr_InitializeByModelKind( this, cellCount ) 
    !-------------------------------------------------------------------------
    ! Specifications
    !-------------------------------------------------------------------------
    implicit none
    class(DispersionDataType) :: this
    integer :: cellCount
    !-------------------------------------------------------------------------

    ! Linear
    select case( this%modelKind )
      case (1)
        if(allocated(this%DMEff)) deallocate(this%DMEff)
        if(allocated(this%AlphaL)) deallocate(this%AlphaL)
        if(allocated(this%AlphaTH)) deallocate(this%AlphaTH)
        if(allocated(this%AlphaTV)) deallocate(this%AlphaTV)
        allocate(this%DMEff(cellCount))
        allocate(this%AlphaL(cellCount))
        allocate(this%AlphaTH(cellCount))
        allocate(this%AlphaTV(cellCount))
      case (2)
        if(allocated(this%DMEff)) deallocate(this%DMEff)
        if(allocated(this%BetaL)) deallocate(this%BetaL)
        if(allocated(this%BetaTH)) deallocate(this%BetaTH)
        if(allocated(this%BetaTV)) deallocate(this%BetaTV)
        if(allocated(this%Delta)) deallocate(this%Delta)
        if(allocated(this%DGrain)) deallocate(this%DGrain)
        allocate(this%DMEff(cellCount))
        allocate(this%BetaL(cellCount))
        allocate(this%BetaTH(cellCount))
        allocate(this%BetaTV(cellCount))
        allocate(this%Delta(cellCount))
        allocate(this%DGrain(cellCount))
    end select

  end subroutine pr_InitializeByModelKind


end module DispersionDataModule
