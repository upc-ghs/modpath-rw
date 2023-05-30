module PrecisionModule
  use iso_fortran_env, only : real32, real64
  implicit none
  !---------------------------------------------------
  ! fp: float precision
#ifdef REAL32
  integer,parameter,public :: fp = real32 ! 4 bytes
#else
  integer,parameter,public :: fp = real64 ! 8 bytes
#endif

end module PrecisionModule
