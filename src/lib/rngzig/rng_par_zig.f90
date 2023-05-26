! Marsaglia & Tsang generator for random normals & random exponentials.
! Translated from C by Alan Miller (amiller@bigpond.net.au)

! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
! random variables', J. Statist. Software, v5(8).

! This is an electronic journal which can be downloaded from:
! http://www.jstatsoft.org/v05/i08

! N.B. It is assumed that all integer(int32)s are 32-bit.
! N.B. The value of M2 has been halved to compensate for the lack of
!      unsigned integer(int32)s in Fortran.

! Latest version - 1 January 2001

! Parallel version - October 2006

! This version has been customised for parallel processing use,
! specifically with OpenMP.  Each thread uses its own pseudo-random 
! sequence. (Gib Bogle)

! Changed to subroutines, cleaned tab characters, 
! unified to lower case, new interface (Vladimir Fuka) -  May 2014

! Changed the underlying integer generator from 32-bit SHR3 to
! 64-bit xoroshiro128+ (David Blackman and Sebastiano Vigna. 
!  Scrambled linear pseudorandom number generators, 2018) translated from C
! from http://xoshiro.di.unimi.it/ . The advantage is faster generating, 
! longer independent series and jump functions available for seeding 
! parallel sequences. Generation of uniform floating point numbers
! now using better and faster method from 
! https://experilous.com/1/blog/post/perfect-fast-random-floating-point-numbers
! Also unified all pre-fixes to rng_. 
! The code does not assume any default integer and can be compiled with
! larger default integer kinds.
! (Vladimir Fuka) - June 2018

!--------------------------------------------------------------------------

module rng_int_to_real
  use iso_fortran_env, only: int32, int64, real32, real64

contains

  ! Generate real numbers from interval [0, 1) from random integers
  ! Translated from C rand_float_co() and rand_float_co() from
  ! https://experilous.com/1/blog/post/perfect-fast-random-floating-point-numbers#half-open-range
  
  function rng_int32_to_real32(i) result(res)
    real(real32) :: res
    integer(int32), value :: i

    i  = ior(int(Z'3F800000',int32), shiftr(i, 9))
    res = transfer(i, 1.0_real32) - 1;
  end function
 
  function rng_int64_to_real32(i) result(res)
    real(real32) :: res
    integer(int64), value :: i
    integer(int32) :: tmp

    tmp = transfer(shiftr(i,32), tmp)
    tmp  = ior(int(Z'3F800000',int32), shiftr(tmp, 9))
    res = transfer(tmp, 1.0_real32) - 1;
  end function
 
  function rng_int64_to_real64(i) result(res)
    real(real64) :: res
    integer(int64), value :: i
    !integer(int32) :: tmp
    
    i  = ior(int(Z'3FF0000000000000',int64), shiftr(i, 12))
    res = transfer(i, 1.0_real64) - 1;
  end function
 
end module

module rng_par_zig

  use iso_fortran_env, only: int32, int64, real32, real64
  use rng_int_to_real

  implicit none

  private

  public  :: rng_init, rng_uni, rng_norm, rng_exp, rng_jump
  
  integer,  parameter  :: sp = real32, dp = real64
  
  real(dp), parameter  :: m1=2147483648.0_dp,   m2=2147483648.0_dp,       &
                          half=0.5_dp
  real(dp)             :: dn0=3.442619855899_dp, tn0=3.442619855899_dp,   &
                          vn=0.00991256303526217_dp,                      &
                          q,                    de0=7.697117470131487_dp, &
                          te0=7.697117470131487_dp,                       &
                          ve=0.003949659822581572_dp

  integer, save :: rng_n = 0, rng_step
  integer(int32), allocatable, save ::  rng_kn(:,:), rng_ke(:,:)
  integer(int64), allocatable, save ::  rng_jsr(:)
  real(dp), allocatable, save :: rng_wn(:,:), rng_fn(:,:), rng_we(:,:), rng_fe(:,:)

  interface rng_uni
    module procedure rng_uni
    module procedure rng_uni_32
    module procedure rng_uni_ser
    module procedure rng_uni_32_ser
  end interface

  interface rng_norm
    module procedure rng_rnor
    module procedure rng_rnor_32
    module procedure rng_rnor_ser
    module procedure rng_rnor_32_ser
  end interface

  interface rng_exp
    module procedure rng_rexp
    module procedure rng_rexp_32
    module procedure rng_rexp_ser
    module procedure rng_rexp_32_ser
  end interface
  
  interface rng_init
    module procedure rng_zigset
  end interface
  
  interface rng_xoroshiro128plus
    module procedure rng_xoroshiro128plus
    module procedure rng_xoroshiro128plus_kpar
    module procedure rng_xoroshiro128plus_kpar_32
  end interface
  
  interface rng_jump
    module procedure rng_jump
    module procedure rng_jump_n
  end interface

  ! Define if signed integer addition must not overflow.
  ! Unsigned integer addition can only be done in C.
!#ifdef STRICT_INTEGER_OVERFLOW    
!  interface 
!  function sum_and_overflow(a, b) result(res) bind(C, name="sum_and_overflow")
!    use, intrinsic :: iso_c_binding
!    integer(c_int32_t) :: res
!    integer(c_int32_t), value :: a, b
!  end function
!  end interface
!#endif    


contains


  subroutine rng_zigset( npar, rng_jsrseed, grainsize)

    integer, intent(in)  :: npar
    integer(int64), intent(in)  :: rng_jsrseed(2)
    integer, intent(in), optional  :: grainsize

    integer(int64) :: s(2)
    integer :: i, kpar
    real(dp) dn, tn, de, te

    rng_n = npar

    if (present(grainsize)) then
      rng_step = grainsize
    else
      rng_step = 32
    end if

    ! First we need to allocate all the non-volatile arrays with the size npar
    allocate(rng_jsr(0:npar*rng_step))
    allocate(rng_kn(0:127,0:npar-1))
    allocate(rng_ke(0:255,0:npar-1))
    allocate(rng_wn(0:127,0:npar-1))
    allocate(rng_fn(0:127,0:npar-1))
    allocate(rng_we(0:255,0:npar-1))
    allocate(rng_fe(0:255,0:npar-1))
    
    s = rng_jsrseed

    ! Now treat each instance separately
    do kpar = 0,npar-1
      !  Set the seed
      rng_jsr(kpar*rng_step:kpar*rng_step+1) = s
      
      ! The next instance's seed is generated by the jump function for xoroshiro128+
      call rng_jump(s)

      !  Tables for RNOR
      dn = dn0
      tn = tn0
      q = vn*exp(half*dn*dn)
      rng_kn(0,kpar) = int((dn/q)*m1, int32)
      rng_kn(1,kpar) = 0
      rng_wn(0,kpar) = q/m1
      rng_wn(127,kpar) = dn/m1
      rng_fn(0,kpar) = 1.0_dp
      rng_fn(127,kpar) = exp( -half*dn*dn )
      do  i = 126, 1, -1
          dn = sqrt( -2.0_dp * log( vn/dn + exp( -half*dn*dn ) ) )       ! dn
          rng_kn(i+1,kpar) = int((dn/tn)*m1, int32)
          tn = dn                                                        ! tn
          rng_fn(i,kpar) = exp(-half*dn*dn)
          rng_wn(i,kpar) = dn/m1
      end do

      !  Tables for Rexp
      de = de0
      te = te0
      q = ve*exp( de )
      rng_ke(0,kpar) = int((de/q)*m2, int32)
      rng_ke(1,kpar) = 0
      rng_we(0,kpar) = q/m2
      rng_we(255,kpar) = de/m2
      rng_fe(0,kpar) = 1.0_dp
      rng_fe(255,kpar) = exp( -de )
      do  i = 254, 1, -1
          de = -log( ve/de + exp( -de ) )                                ! de
          rng_ke(i+1,kpar) = int(m2 * (de/te), int32)
          te = de                                                        ! te
          rng_fe(i,kpar) = exp( -de )
          rng_we(i,kpar) = de/m2
      end do
    enddo
  end subroutine rng_zigset


  
  ! Generate random 64-bit integers
  subroutine rng_xoroshiro128plus_kpar(ival, kpar)
    integer(int64), intent(out)   :: ival
    integer, intent(in) :: kpar
    
    call rng_xoroshiro128plus(ival, rng_jsr(kpar*rng_step:kpar*rng_step+1))
  end subroutine

  ! Generate random 32-bit integers by taking the *higher* 32-bits of the result
  subroutine rng_xoroshiro128plus_kpar_32(ival, kpar)
    integer(int32), intent(out)   :: ival
    integer, intent(in) :: kpar
    integer(int64) :: tmp
        
    call rng_xoroshiro128plus(tmp, rng_jsr(kpar*rng_step:kpar*rng_step+1))
    
    ival = transfer(shiftr(tmp,32), ival)
  end subroutine

  ! Generate random 64-bit integers
  subroutine rng_xoroshiro128plus(ival, s) 
    integer(int64), intent(out)   :: ival
    integer(int64), intent(inout) :: s(2)
    integer(int64) :: s1, s2
    
    s1 = s(1)
    s2 = s(2)
    
!#ifdef STRICT_INTEGER_OVERFLOW    
!    ival = sum_and_overflow(s1, s2)
!#else
    ival = s1 + s2
!#endif
    s2 = ieor(s2, s1)
    s1 = ieor( ieor(rotl(s1, 24), s2), shiftl(s2, 16))
    s2 = rotl(s2, 37)
 
    s(1) = s1
    s(2) = s2

  contains
    function rotl(x, k)
      integer(int64) :: rotl
      integer(int64) :: x
      integer :: k
    
      rotl = ior( shiftl(x, k), shiftr(x, 64-k))
    end function
  end subroutine

  ! Jump by 2^64 xoroshiro128plus calls
  subroutine rng_jump(s)
    integer(int64), intent(inout) :: s(2)
    ! The first constant is df900294d8f554a5 using not() to avoid overflow
    integer(int64), parameter :: jmp(2) = [ not(int(Z'206ffd6b270aab5a', int64)), &
                                            int(Z'170865df4b3201fc', int64)]
    integer(int64) :: s1, s2, dummy
    
    integer :: i, b
    
    s1 = 0
    s2 = 0
    
    do i = 1, size(jmp)
      do b = 0, 63
        if (iand(jmp(i), shiftl(1_int64,b))/=0) then
          s1 = ieor(s1, s(1))
          s2 = ieor(s2, s(2))
        end if  
        call rng_xoroshiro128plus(dummy, s)
      end do
    end do
    s(1) = s1
    s(2) = s2
  end subroutine
  
  ! Do n jumps defined above
  subroutine rng_jump_n(s, n)
    integer(int64), intent(inout) :: s(2)
    integer, intent(in) :: n
    integer :: i
    
    do i = 1, n
      call rng_jump(s)
    end do
  end subroutine

  !  Generate uniformly distributed random numbers, sequence kpar
  subroutine rng_uni(fn_val, kpar)
    integer :: kpar
    real(dp), intent(out) ::  fn_val
    integer(int64) :: x

    if (kpar >= rng_n) then
        write(*,*) 'RNG Error: thread number',kpar, 'exceeds initialized max: ',rng_n-1
        write(*,*) 'or RNG not initialized.'
        stop
    endif
    
    call rng_xoroshiro128plus(x, kpar)
    
    fn_val = rng_int64_to_real64(x)
  end subroutine rng_uni



  !  Generate uniformly distributed random numbers, sequence kpar
  subroutine rng_uni_32(fn_val, kpar)
    integer :: kpar
    real(sp), intent(out) ::  fn_val
    integer(int64) :: x

    if (kpar >= rng_n) then
        write(*,*) 'RNG Error: thread number',kpar, 'exceeds initialized max: ',rng_n-1
        write(*,*) 'or RNG not initialized.'
        stop
    endif
    
    call rng_xoroshiro128plus(x, kpar)
    
    fn_val = rng_int64_to_real32(x)
  end subroutine rng_uni_32



  !  Generate random normals, sequence kpar
  subroutine rng_rnor(fn_val, kpar)
    real(dp), intent(out) ::  fn_val
    integer :: kpar

    real(dp), parameter  ::  r = 3.442620_dp
    real(dp)             ::  x, y, z
    integer(int32) :: iz, hz

    if (kpar >= rng_n) then
        write(*,*) 'RNG Error: thread number',kpar, 'exceeds initialized max: ',rng_n-1
        write(*,*) 'or RNG not initialized.'
        stop
    endif
    
    call rng_xoroshiro128plus(hz, kpar)
    iz = iand( hz, 127 )
    
    if( abs( hz ) < rng_kn(iz,kpar) ) then
        fn_val = hz * rng_wn(iz,kpar)
    else
        do
          if( iz == 0 ) then
              do
                call rng_uni(z, kpar)
                x = -0.2904764_dp * log( z )
                call rng_uni(z, kpar)
                y = -log( z )
                if( y+y >= x*x ) exit
              end do
              fn_val = r+x
              if( hz <= 0 ) fn_val = -fn_val
              return
          end if
          
          x = hz * rng_wn(iz,kpar)
          
          call rng_uni(z, kpar)
          if( rng_fn(iz,kpar) + z*(rng_fn(iz-1,kpar)-rng_fn(iz,kpar)) < exp(-half*x*x) ) then
              fn_val = x
              return
          end if
          
          call rng_xoroshiro128plus(hz, kpar)
          iz = iand( hz, 127 )
          if( abs( hz ) < rng_kn(iz,kpar) ) then
              fn_val = hz * rng_wn(iz,kpar)
              return
          end if
        end do
    end if
  end subroutine rng_rnor



  !  Generate random exponentials, sequence kpar
  subroutine rng_rexp(fn_val, kpar)
    real(dp), intent(out)  ::  fn_val
    integer :: kpar

    real(dp)  ::  x, y
    integer(int32) :: iz, jz

    if (kpar >= rng_n) then
        write(*,*) 'RNG Error: thread number',kpar, 'exceeds initialized max: ',rng_n-1
        write(*,*) 'or RNG not initialized.'
        stop
    endif

    call rng_xoroshiro128plus(jz, kpar)
    iz = iand( jz, 255 )
    if( abs( jz ) < rng_ke(iz,kpar) ) then
        fn_val = abs(jz) * rng_we(iz,kpar)
        return
    end if

    do
        if( iz == 0 ) then
          call rng_uni(y, kpar)
          fn_val = 7.69711_dp - log( y )
          return
        end if
        x = abs( jz ) * rng_we(iz,kpar)
        
        call rng_uni(y, kpar)
        
        if( rng_fe(iz,kpar) + y * (rng_fe(iz-1,kpar) - rng_fe(iz,kpar)) < exp( -x ) ) then
          fn_val = x
          return
        end if
        
        call rng_xoroshiro128plus(jz, kpar)
        iz = iand( jz, 255 )
        if( abs( jz ) < rng_ke(iz,kpar) ) then
          fn_val = abs( jz ) * rng_we(iz,kpar)
          return
        end if
    end do
  end subroutine rng_rexp
  
  subroutine rng_rnor_32(x, kpar)
    real(sp), intent(out) :: x
    integer :: kpar
    real(dp) :: x64
    call rng_rnor(x64, kpar)
    x = real(x64, sp)
  end subroutine
  
  subroutine rng_rexp_32(x, kpar)
    real(sp), intent(out) :: x
    integer :: kpar
    real(dp) :: x64
    call rng_rexp(x64, kpar)
    x = real(x64, sp)
  end subroutine
  
  subroutine rng_uni_ser(x)
    real(dp), intent(out) :: x

    call rng_uni(x, 0)
  end subroutine
  
  subroutine rng_rnor_ser(x)
    real(dp), intent(out) :: x

    call rng_rnor(x, 0)
  end subroutine
  
  subroutine rng_rexp_ser(x)
    real(dp), intent(out) :: x

    call rng_rexp(x, 0)
  end subroutine
  
  subroutine rng_uni_32_ser(x)
    real(sp), intent(out) :: x
    real(dp) :: x64
    call rng_uni(x64, 0)
    x = real(x64, sp)
  end subroutine
  
  subroutine rng_rnor_32_ser(x)
    real(sp), intent(out) :: x
    real(dp) :: x64
    call rng_rnor(x64, 0)
    x = real(x64, sp)
  end subroutine
  
  subroutine rng_rexp_32_ser(x)
    real(sp), intent(out) :: x
    real(dp) :: x64
    call rng_rexp(x64, 0)
    x = real(x64, sp)
  end subroutine
  
end module rng_par_zig
