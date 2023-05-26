module ParticleModule
  implicit none

  ! Set default access status to private
  private

  type,public :: ParticleType
    ! RW: The Face property is a potential candidate for deprecation
    integer :: CellNumber, Layer, Face, Group, ID, Status, SequenceNumber
    doubleprecision :: LocalX, LocalY, LocalZ, GlobalZ, TrackingTime
    integer :: InitialCellNumber, InitialFace, Drape, InitialLayer
    doubleprecision :: InitialLocalX, InitialLocalY, InitialLocalZ, InitialGlobalZ, InitialTrackingTime
    ! GPKDE-RECONSTRUCTION
    doubleprecision :: GlobalX, GlobalY
    ! RWPT
    doubleprecision :: Mass = 1d0
  end type


contains


end module ParticleModule
