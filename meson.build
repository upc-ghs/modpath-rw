## Project ##
project(
  'MODPATH-RW',
  'fortran',
  version: '1.0.0',
  meson_version: '>= 0.59.0',
  default_options : [
    'b_vscrt=static_from_buildtype', # Link runtime libraries statically on Windows
    'optimization=2',
    'fortran_std=f2008',
    'warning_level=0',
])

## Parameters ##

# Optimization
if get_option('optimization') == '3'
  error('Only optimization levels <= 2 are supported')
endif
if get_option('optimization') == '2'
  profile = 'release'
else
  profile = 'develop'
endif
message('The used profile is:', profile)

# OpenMP
omp = dependency('openmp', language:'fortran')
dependencies = [ ]
if omp.found()
  with_omp = true
  dependencies = [ omp ]
else
  with_omp = false
endif 
message('Compile with OpenMP:', with_omp)

# Compiler
fc    = meson.get_compiler('fortran')
fc_id = fc.get_id()
compile_args = []
link_args    = []

# Command line options for gfortran
if fc_id == 'gcc'
  # General options
  compile_args += [
                   '-fbacktrace',
                   '-fall-intrinsics',
                   '-pedantic',
                   '-cpp',
                   '-Wline-truncation',
                   '-Wtabs',
                   '-Wunused-label',
                   '-Wunused-variable',
                   '-Wcharacter-truncation',
                   ]
  if with_omp
    compile_args += ['-fopenmp']
  endif

  # Options specific to profile
  if profile == 'release'
    compile_args += ['-ffpe-summary=overflow', '-ffpe-trap=overflow,zero,invalid']
  elif profile == 'develop'
    compile_args += ['-fcheck=all', '-ffpe-trap=overflow,zero,invalid,denormal']
  endif

  # Define OS with gfortran for OS specific code
  # These are identical to pre-defined macros available with ifort
  system = build_machine.system()
  if system == 'linux'
    compile_args += '-D__linux__'
  elif system == 'darwin'
    compile_args += '-D__APPLE__'
  elif system == 'windows'
    compile_args += '-D_WIN32'
  endif
endif

# Command line options for ifort
if fc_id == 'intel-cl'
  # windows
  compile_args += ['/fpe:0',              # Activate all floating point exceptions
                   '/heap-arrays:0',
                   '/traceback',
                   '/fpp',                # Activate preprocessing
                   '/Qdiag-disable:7416', # f2008 warning
                   '/Qdiag-disable:7025', # f2008 warning
                   '/Qdiag-disable:5268', # Line too long
                  ]
  link_args += ['/ignore:4217',  # access through ddlimport might be inefficient
                '/ignore:4286'   # same as 4217, but more general
               ]
  if with_omp
    compile_args += ['/Qopenmp'] # OpenMP
  endif

elif fc_id == 'intel'
  # linux and macOS
  compile_args += ['-fpe0',              # Activate all floating point exceptions
                   '-no-heap-arrays',
                   '-traceback',
                   '-fpp',                # Activate preprocessing
                   '-diag-disable:7416',  # f2008 warning
                   '-diag-disable:7025',  # f2008 warning
                   '-diag-disable:5268',  # Line too long
                  ]
  link_args += '-static-intel'

  if with_omp
    compile_args += ['-qopenmp'] # OpenMP
  endif

endif

add_project_arguments(fc.get_supported_arguments(compile_args), language: 'fortran')
add_project_link_arguments(fc.get_supported_arguments(link_args), language: 'fortran')

# Load the src
subdir('src')

# Tests to evaluate installation success
test('Test version', theexe, args : ['-v',])
test('Test help   ', theexe, args : ['-h',])
