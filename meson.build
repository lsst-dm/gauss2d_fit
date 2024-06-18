project(
  'libgauss2dfit',
  'cpp',
  version : files('VERSION'),
  license : 'GPL-3.0-or-later',
  default_options : ['cpp_std=c++17',],
)

eups = find_program('eups', required: false)
use_eups = eups.found()
if use_eups
	use_eups = run_command('eups', 'list', 'gauss2d_fit', check: false).returncode() == 0
endif

gsl_dep = dependency('gsl', required: false)

if gsl_dep.found()
    add_project_arguments('-DLSST_GAUSS2D_FIT_HAS_GSL', language : 'cpp')
endif

if use_eups
    gauss2d_eupsdir = run_command('sh', '-c', 'echo $GAUSS2D_DIR', check: true).stdout().strip() / 'build-release'
    cpp = meson.get_compiler('cpp')
    gauss2d_dep = declare_dependency(
        dependencies: [cpp.find_library('gauss2d', dirs: [gauss2d_eupsdir / 'lib64/lsst']),],
        include_directories: gauss2d_eupsdir / 'include',
    )

    parameters_eupsdir = run_command('sh', '-c', 'echo $MODELFIT_PARAMETERS_DIR', check: true).stdout().strip() / 'build-release'
    parameters_dep = declare_dependency(
        include_directories: parameters_eupsdir / 'include',
    )
else
    gauss2d_dep = dependency('lsst_gauss2d')
    parameters_dep = dependency('lsst_modelfit_parameters')
endif

inc = include_directories('include')

depends = declare_dependency(
    include_directories: inc,
    dependencies: [gauss2d_dep, gsl_dep, parameters_dep],
)

subdir('include')
subdir('src')
subdir('tests')

if not use_eups
    pkg_mod = import('pkgconfig')
    if pkg_mod.found()
        pkg_mod.generate(libraries : gauss2dfitlib,
                         version : '0.1',
                         name : 'libgauss2d_fit',
                         filebase : 'lsst_gauss2d_fit',
                         description : 'Create, manipulate and evaluate 2D Gaussian mixtures and images thereof.')
    endif
endif

# Build doxygen docs (done by eupspkg install if use_eups)
if not use_eups
    pymod = import('python')
    python = pymod.find_installation(
        'python3',
        required : false,
        modules: ['lsst.sconsUtils'],
    )
    has_py3 = python.found()
    scons = find_program('scons', required: false)
    has_scons = scons.found()
    if has_scons and has_py3
        message('Running SCons to build docs; see ' + meson.build_root() + '/meson-logs/meson-log.txt for output')
        run_command('sh', '-c', 'cd $MESON_SOURCE_ROOT && echo $PWD && scons', check: false, capture: true)
    endif
endif
