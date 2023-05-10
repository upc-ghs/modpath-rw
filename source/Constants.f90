module ConstantsModule
  use PrecisionModule, only : fp 
  implicit none
  !---------------------------------------------------
  real(fp), parameter, public :: fZERO       = 0.0_fp
  real(fp), parameter, public :: fONE        = 1.0_fp
  real(fp), parameter, public :: fTWO        = 2.0_fp
  real(fp), parameter, public :: fTHREE      = 3.0_fp
  real(fp), parameter, public :: fFOUR       = 4.0_fp
  real(fp), parameter, public :: fFIVE       = 5.0_fp
  real(fp), parameter, public :: fSIX        = 6.0_fp
  real(fp), parameter, public :: fEIGHT      = 8.0_fp
  real(fp), parameter, public :: pi          = fFOUR*atan(fONE)
  real(fp), parameter, public :: sqrtPi      = sqrt(fFOUR*atan(fONE))
  real(fp), parameter, public :: sqrtTwo     = sqrt(fTWO)
  real(fp), parameter, public :: sqrtEightPi = sqrt(fEIGHT*pi)
end module ConstantsModule
