project(
  'libgauss2dfit',
  'cpp',
  version : run_command('cat', meson.source_root() + '/VERSION', check: true).stdout().strip(),
  license : 'GPL-3.0-or-later',
  default_options : ['cpp_std=c++17',],
)

eups = find_program('eups', required: false)
use_eups = eups.found()
if use_eups
	use_eups = run_command('eups', 'list', 'gauss2dfit', check: false).returncode() == 0
endif

gsl_dep = dependency('gsl', required: false)

if use_eups
    gauss2d_eupsdir = run_command('eups', 'list', '-d', 'gauss2d', check: true).stdout().strip()
    cpp = meson.get_compiler('cpp')
    gauss2d_dep = declare_dependency(
        dependencies: [cpp.find_library('gauss2d', dirs: [gauss2d_eupsdir + '/build-release/lib64']),],
        include_directories: gauss2d_eupsdir + '/build-release/include',
    )

    parameters_eupsdir = run_command('eups', 'list', '-d', 'parameters', check: true).stdout().strip()
    parameters_dep = declare_dependency(
        include_directories: parameters_eupsdir + '/build-release/include',
    )
else
    gauss2d_dep = dependency('gauss2d')
    parameters_dep = dependency('parameters')
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
                         name : 'libgauss2dfit',
                         filebase : 'gauss2dfit',
                         description : 'Create, manipulate and evaluate 2D Gaussian mixtures and images thereof.')
    endif
endif

subdir('doc')
