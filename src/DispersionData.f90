module DispersionDataModule
  implicit none
  
  ! Set default access status to private
  private
  !
  type, public :: DispersionDataType
    !
    integer            :: id
    character(len=20)  :: stringid
    !
    ! -- 0: linear isotropic 
    ! -- 1: linear axisymmetric 
    integer :: modelKind = 1
    !  
    ! -- 0: uniform 
    ! -- 1: distributed, u3d reader
    integer :: parametersFormat 
    logical :: uniformParameters 
    !  
    ! -- Dispersivities for linear dispersion model
    doubleprecision,dimension(:),allocatable :: AlphaL  ! AlphaL is pointing here
    doubleprecision,dimension(:),allocatable :: AlphaLV
    doubleprecision,dimension(:),allocatable :: AlphaTH ! AlphaT is pointing here 
    doubleprecision,dimension(:),allocatable :: AlphaTV
    doubleprecision,dimension(:),allocatable :: DMEff   ! Effective, corrected by tortuosity
    !
    ! -- Parameters for nonlinear dispersion model (not implemented)
    doubleprecision :: dmaqueous = 0d0
    doubleprecision,dimension(:),allocatable :: BetaL
    doubleprecision,dimension(:),allocatable :: BetaTH
    doubleprecision,dimension(:),allocatable :: BetaTV
    doubleprecision,dimension(:),allocatable :: Delta
    doubleprecision,dimension(:),allocatable :: DGrain
    !
    ! -- initialized flag
    logical :: initialized =.false.
    !
  contains
    procedure :: Reset => pr_Reset
    procedure :: InitializeParameters => pr_InitializeParameters
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
    if(allocated(this%AlphaLV)) deallocate(this%AlphaLV)
    if(allocated(this%AlphaTH)) deallocate(this%AlphaTH)
    if(allocated(this%AlphaTV)) deallocate(this%AlphaTV)
    if(allocated(this%DMEff)) deallocate(this%DMEff)
    if(allocated(this%BetaL)) deallocate(this%BetaL)
    if(allocated(this%BetaTH)) deallocate(this%BetaTH)
    if(allocated(this%BetaTV)) deallocate(this%BetaTV)
    if(allocated(this%Delta)) deallocate(this%Delta)
    if(allocated(this%DGrain)) deallocate(this%DGrain)

  end subroutine pr_Reset


  subroutine pr_InitializeParameters( this, cellCount ) 
    !-------------------------------------------------------------------------
    ! Specifications
    !-------------------------------------------------------------------------
    implicit none
    class(DispersionDataType) :: this
    integer :: cellCount
    !-------------------------------------------------------------------------

    select case(this%parametersFormat)
    ! spatially uniform
    case(0)
     ! Set the flag
     this%uniformParameters = .true.
     !
     select case( this%modelKind )
     case (0) ! -- Linear isotropic 
       if(allocated(this%DMEff)) deallocate(this%DMEff)
       if(allocated(this%AlphaL)) deallocate(this%AlphaL)
       if(allocated(this%AlphaTH)) deallocate(this%AlphaTH)
       allocate(this%DMEff(1))
       allocate(this%AlphaL(1))
       allocate(this%AlphaTH(1))
     case (1) ! -- Linear axisymmetric
       if(allocated(this%DMEff)) deallocate(this%DMEff)
       if(allocated(this%AlphaL)) deallocate(this%AlphaL)
       if(allocated(this%AlphaLV)) deallocate(this%AlphaLV)
       if(allocated(this%AlphaTH)) deallocate(this%AlphaTH)
       if(allocated(this%AlphaTV)) deallocate(this%AlphaTV)
       allocate(this%DMEff(1))
       allocate(this%AlphaL(1))
       allocate(this%AlphaLV(1))
       allocate(this%AlphaTH(1))
       allocate(this%AlphaTV(1))
     end select

    ! distributed
    case(1)
     ! Set the flag
     this%uniformParameters = .false.

     select case( this%modelKind )
       case (0) ! -- Linear isotropic 
         if(allocated(this%DMEff)) deallocate(this%DMEff)
         if(allocated(this%AlphaL)) deallocate(this%AlphaL)
         if(allocated(this%AlphaTH)) deallocate(this%AlphaTH)
         allocate(this%DMEff(cellCount))
         allocate(this%AlphaL(cellCount))
         allocate(this%AlphaTH(cellCount))
       case (1) ! -- Linear axisymmetric
         if(allocated(this%DMEff)) deallocate(this%DMEff)
         if(allocated(this%AlphaL)) deallocate(this%AlphaL)
         if(allocated(this%AlphaLV)) deallocate(this%AlphaLV)
         if(allocated(this%AlphaTH)) deallocate(this%AlphaTH)
         if(allocated(this%AlphaTV)) deallocate(this%AlphaTV)
         allocate(this%DMEff(cellCount))
         allocate(this%AlphaL(cellCount))
         allocate(this%AlphaLV(cellCount))
         allocate(this%AlphaTH(cellCount))
         allocate(this%AlphaTV(cellCount))
     end select
    end select

  end subroutine pr_InitializeParameters


end module DispersionDataModule
