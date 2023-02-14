module ObservationModule
  implicit none
  
  ! Set default access status to private
  private

  type, public :: ObservationType
    integer, allocatable, dimension(:) :: cells
    integer                            :: nCells
    integer, allocatable, dimension(:) :: nRecordsCell
    integer                            :: nAuxRecords = 0
    integer                            :: id
    character(len=20)                  :: stringid
    integer                            :: style
    character(len=20)                  :: stylestringid
    integer                            :: cellOption
    integer                            :: timeOption
    integer                            :: outputUnit
    character(len=200)                 :: outputFileName
    integer                            :: recOutputUnit
    character(len=200)                 :: recOutputFileName
    integer                            :: auxOutputUnit
    character(len=200)                 :: auxOutputFileName
  end type


contains

  ! Some subroutines


end module ObservationModule
