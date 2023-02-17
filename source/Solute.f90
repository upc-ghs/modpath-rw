module SoluteModule
  use DispersionDataModule,only : DispersionDataType

  implicit none
  
  ! Set default access status to private
  private

  type, public :: SoluteType
    ! Id's
    integer                            :: id      ! internal
    integer                            :: userid  ! given by the user ( used ? )
    character(len=20)                  :: stringid
    logical                            :: initialized =.false.
    ! Disperion model foreign key
    class(DispersionDataType), pointer :: dispersion
    integer                            :: dispersionId 
    character(len=20)                  :: dispersionStringId
    ! PGroups foreign key
    integer                            :: nParticleGroups = 0
    integer, dimension(:), allocatable :: pGroups
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
      class(SoluteType) :: this
      integer :: cellCount
      !-------------------------------------------------------------------------

      this%initialized = .true.

    end subroutine


    subroutine pr_Reset( this )
      !-------------------------------------------------------------------------
      ! Specifications
      !-------------------------------------------------------------------------
      implicit none
      class(SoluteType) :: this
      !-------------------------------------------------------------------------

      this%id       = 0
      this%stringid = '' 
      this%initialized = .false.

    end subroutine pr_Reset


end module SoluteModule
