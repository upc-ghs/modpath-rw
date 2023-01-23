module DIspersionDataModule

  implicit none
  
  ! Set default access status to private
  private

  type, public :: DispersionDataType
    integer :: id
    character(len=300) :: stringid
    integer :: dispersionModel = 0

    ! Dispersivities for linear dispersion model
    doubleprecision,dimension(:),allocatable :: AlphaL
    doubleprecision,dimension(:),allocatable :: AlphaTH
    doubleprecision,dimension(:),allocatable :: AlphaTV
    ! Effective molecular diffusion, corrected by tortuosity,
    ! in case is spatially distributed 
    doubleprecision,dimension(:),allocatable :: DMEff
    ! Ideally there is something to identify whether
    ! effective molecular diffusion is distributed or not


    ! Transport properties
    doubleprecision :: dAqueous
    doubleprecision :: aqueousDiffusion
    doubleprecision :: poreDiffusion  ! or effective diffusion
    doubleprecision :: effectiveDiffusion  ! Aqueous diffusion with tortuosity correction

    logical :: initialized =.false.

  contains
    procedure :: Initialize => pr_Initialize
    procedure :: Reset => pr_Reset
  end type


contains


    subroutine pr_Initialize( this, cellCount ) 
      !-------------------------------------------------------------------------
      ! Specifications
      !-------------------------------------------------------------------------
      implicit none
      class(DispersionDataType) :: this
      integer :: cellCount
      !-------------------------------------------------------------------------

      if(allocated(this%AlphaL)) deallocate(this%AlphaL)
      if(allocated(this%AlphaTH)) deallocate(this%AlphaTH)
      if(allocated(this%AlphaTV)) deallocate(this%AlphaTV)
      allocate(this%AlphaL(cellCount))
      allocate(this%AlphaTH(cellCount))
      allocate(this%AlphaTV(cellCount))

      this%initialized = .true.


    end subroutine


    subroutine pr_Reset( this )
      !-------------------------------------------------------------------------
      ! Specifications
      !-------------------------------------------------------------------------
      implicit none
      class(DispersionDataType) :: this
      !-------------------------------------------------------------------------

      this%initialized = .false.
      this%id       = 0
      this%stringid = '' 
      if(allocated(this%AlphaL)) deallocate(this%AlphaL)
      if(allocated(this%AlphaTH)) deallocate(this%AlphaTH)
      if(allocated(this%AlphaTV)) deallocate(this%AlphaTV)
      if(allocated(this%DMEff)) deallocate(this%DMEff)

    end subroutine pr_Reset


end module DispersionDataModule
