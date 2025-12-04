module CompilerVersion
  use iso_fortran_env, only: compiler_version
  implicit none
  private
  ! -- compiler version
  character(len=10) :: ccompiler
  character(len=20) :: cversion
  character(len=20) :: cdate
  integer :: icompiler = 0
  integer :: iversion = 0
  integer :: imajor = 0
  integer :: iminor = 0
  integer :: imicro = 0
  public :: get_compiler
  public :: get_compiler_txt

contains
  
  subroutine get_compiler_txt(txt)
  character(len=90), intent(inout) :: txt
        
    ! -- set variables
#ifdef __GFORTRAN__ 
    icompiler = 1
    cversion = __VERSION__
    cdate = __DATE__ // ' ' // __TIME__
#endif
#ifdef __INTEL_COMPILER
    icompiler = 2
    iversion = __INTEL_COMPILER
    cdate = __DATE__ // ' ' // __TIME__
    imicro = __INTEL_COMPILER_UPDATE
#endif
    !
    if (icompiler < 1) then
      ccompiler = 'UNKNOWN'
      cversion = '??.??'
      cdate = '??? ?? ???? ??:??:??'
    else if (icompiler == 1) then
      ccompiler = 'GFORTRAN'
    else if (icompiler == 2) then
      ccompiler = 'INTEL'
    end if
    !
    ! -- write compiler version 
    write (txt, '(a,3(1x,a))') &
      'Program compiled', trim(adjustl(cdate)), &
      'with', trim(adjustl(compiler_version()))

  end subroutine get_compiler_txt

  subroutine get_compiler(compiler)
  character(len=10), intent(inout) :: compiler
    
    ! -- set variables
#ifdef __GFORTRAN__ 
    icompiler = 1
    cversion = __VERSION__
    cdate = __DATE__ // ' ' // __TIME__
#endif
#ifdef __INTEL_COMPILER
    icompiler = 2
    iversion = __INTEL_COMPILER
    cdate = __DATE__ // ' ' // __TIME__
    imicro = __INTEL_COMPILER_UPDATE
#endif

    if (icompiler < 1) then
      compiler = 'UNKNOWN'
    else if (icompiler == 1) then
      compiler = 'GFORTRAN'
    else if (icompiler == 2) then
      compiler = 'INTEL'
    end if
    
  end subroutine get_compiler

end module CompilerVersion
