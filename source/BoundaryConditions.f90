module BoundaryConditionsModule
  implicit none
  
  ! Set default access status to private
  private

  type, public :: PrescribedType
    !---------------------------------------------
    ! Manages prescribed concentrations 
    !--------------------------------------------
    integer                            :: id
    character(len=200)                 :: stringid
    integer, allocatable, dimension(:) :: cells
    integer                            :: nCells
    integer, allocatable, dimension(:) :: nRecordsCell
    integer                            :: nAuxRecords = 0
    integer                            :: style
    integer                            :: cellOption
    integer                            :: timeOption
    character(len=300)                 :: outputFileName
    integer                            :: outputUnit
    integer                            :: auxOutputUnit
    character(len=300)                 :: auxOutputFileName

    ! Prototyping
    doubleprecision :: particlesMass
    doubleprecision :: concentration
    integer         :: iNPX, iNPY, iNPZ
    integer         :: nParticles
    doubleprecision :: totalPrescribedMass
    doubleprecision :: totalMassCounter = 0d0
   
    ! It requires a soluteId or solute foreign key

  end type


contains

  ! Some subroutines


end module BoundaryConditionsModule
