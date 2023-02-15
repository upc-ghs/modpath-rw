module ParticleModule
  implicit none
  
! Set default access status to private
  private

type,public :: ParticleType
  integer :: CellNumber, Layer, Face, Group, ID, Status, SequenceNumber
  doubleprecision :: LocalX, LocalY, LocalZ, GlobalZ, TrackingTime
  integer :: InitialCellNumber, InitialFace, Drape, InitialLayer
  doubleprecision :: InitialLocalX, InitialLocalY, InitialLocalZ, InitialGlobalZ, InitialTrackingTime
  doubleprecision,dimension(:),allocatable :: ExitVelocity
  ! GPKDE-RECONSTRUCTION
  doubleprecision :: GlobalX, GlobalY
  ! RWPT
  doubleprecision :: Mass   = 1d0
  integer         :: Solute = 1
  ! may be deprecated
  doubleprecision :: DAqueous

end type








contains




end module ParticleModule
