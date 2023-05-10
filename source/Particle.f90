module ParticleModule
  use PrecisionModule, only : fp
  implicit none

  ! Set default access status to private
  private

  type,public :: ParticleType
    ! RW: The Face property is a potential candidate for deprecation
    integer :: CellNumber, Layer, Face, Group, ID, Status, SequenceNumber
    real(fp) :: LocalX, LocalY, LocalZ, GlobalZ, TrackingTime
    integer :: InitialCellNumber, InitialFace, Drape, InitialLayer
    real(fp) :: InitialLocalX, InitialLocalY, InitialLocalZ, InitialGlobalZ, InitialTrackingTime
    ! GPKDE-RECONSTRUCTION
    real(fp) :: GlobalX, GlobalY
    ! RWPT
    real(fp) :: Mass = 1.0_fp
  end type


contains


end module ParticleModule
