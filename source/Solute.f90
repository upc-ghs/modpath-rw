module SoluteModule

  implicit none
  
  ! Set default access status to private
  private

  type, public :: SoluteType
    integer :: id
    character(len=300) :: stringid
    ! Transport properties
    doubleprecision :: dAqueous
    doubleprecision :: aqueousDiffusion
    doubleprecision :: poreDiffusion  ! or effective diffusion
    doubleprecision :: effectiveDiffusion  ! Aqueous diffusion with tortuosity correction
    logical :: initialized =.false.
    integer :: dispersionModel = 0
    integer :: nParticleGroups = 0
    integer, dimension(:), allocatable :: pGroups
    ! Dispersivities ?
    ! Not necessarily consistent with definition 
    ! of dispersivity as medium property, but
    ! allows to implement species specific dispersion
    doubleprecision,dimension(:),allocatable :: AlphaLong
    doubleprecision,dimension(:),allocatable :: AlphaTran
    doubleprecision :: betaLong, betaTrans

  contains
    procedure :: Initialize => pr_Initialize
    procedure :: Reset => pr_Reset
  end type


contains

!!!!!!!! Some subroutines

    subroutine pr_Initialize( this, cellCount ) 
      !*************************************************************************
      !
      !*************************************************************************
      ! Specifications
      !-------------------------------------------------------------------------
      implicit none
      class(SoluteType) :: this
      integer :: cellCount
      !-------------------------------------------------------------------------

      if(allocated(this%AlphaLong)) deallocate(this%AlphaLong)
      if(allocated(this%AlphaTran)) deallocate(this%AlphaTran)
      allocate(this%AlphaTran(cellCount))
      allocate(this%AlphaLong(cellCount))

      this%initialized = .true.


    end subroutine


    subroutine pr_Reset( this )
      !*************************************************************************
      !
      !*************************************************************************
      ! Specifications
      !-------------------------------------------------------------------------
      implicit none
      class(SoluteType) :: this
      !-------------------------------------------------------------------------

      this%dAqueous = 0d0
      this%id       = 0
      this%stringid = '' 
      if(allocated(this%AlphaLong)) deallocate(this%AlphaLong)
      if(allocated(this%AlphaTran)) deallocate(this%AlphaTran)
      this%initialized = .false.

    end subroutine pr_Reset


end module SoluteModule
