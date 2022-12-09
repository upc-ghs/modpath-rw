module ObservationModule
  implicit none
  
  ! Set default access status to private
  private

  type, public :: ObservationType
    integer, allocatable, dimension(:) :: cells
    integer                            :: nCells
    integer, allocatable, dimension(:) :: nRecordsCell
    integer                            :: id
    integer                            :: style
    integer                            :: cellOption
    integer                            :: timeOption
    character(len=300)                 :: outputFileName
    integer                            :: outputUnit
  end type


contains

  ! Some subroutines


end module ObservationModule
