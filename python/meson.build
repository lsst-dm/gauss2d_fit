project(
    'gauss2dfit',
    'cpp',
    default_options : ['cpp_std=c++17',],
    version : run_command('cat', meson.source_root() + '/VERSION', check: true).stdout().strip(),
)

add_project_arguments('-DVERSION=@0@'.format(meson.project_version()), language : 'cpp')

eups = find_program('eups', required: false)
use_eups = eups.found()

if use_eups
	gauss2dfit_eups = run_command('eups', 'list', '-d', 'gauss2dfit', check: false)
	use_eups = gauss2dfit_eups.returncode() == 0
endif

gsl_dep = dependency('gsl', required: false)

if gsl_dep.found()
    add_project_arguments('-DGAUSS2D_FIT_HAS_GSL', language : 'cpp')
endif

if use_eups
    gauss2d_eupsdir = run_command('sh', '-c', 'echo $GAUSS2D_DIR', check: true).stdout().strip()
    cpp = meson.get_compiler('cpp')
    gauss2d_dep = declare_dependency(
        dependencies: [cpp.find_library('gauss2d', dirs: [gauss2d_eupsdir + '/build-release/lib64']),],
        include_directories: gauss2d_eupsdir + '/build-release/include',
    )
    gauss2d_python_dep = declare_dependency(
        include_directories: gauss2d_eupsdir + '/python/lib/build-release/include',
    )
    gauss2dfit_eupsdir = gauss2dfit_eups.stdout().strip()
    parameters_eupsdir = run_command('sh', '-c', 'echo $MODELFIT_PARAMETERS_DIR', check: true).stdout().strip()
    # TODO: Can this dependency be extracted from '../meson.build' somehow
    gauss2dfit_dep = declare_dependency(
        dependencies: [cpp.find_library('gauss2dfit', dirs: [gauss2dfit_eupsdir + '/build-release/lib64']),],
        include_directories: [
            gauss2dfit_eupsdir + '/build-release/include',
            parameters_eupsdir + '/build-release/include',
        ]
    )
else
    gauss2d_dep = dependency('gauss2d')
    gauss2d_python_dep = dependency('gauss2d_python')
    gauss2dfit_dep = dependency('gauss2dfit')
endif

pybind11_dep = dependency('pybind11')
pymod = import('python')
python = pymod.find_installation(
    'python3',
    required : true,
    modules: ['gauss2d', 'pytest'],
)
subdir('gauss2d/fit')
subdir('tests')

if not use_eups
	pip_cmd = [
	    'python', '-m', 'pip', 'install', meson.source_root(),
	     '-t', python.get_install_dir(),
	     '--no-dependencies',
	     '--no-build-isolation',
	     '--upgrade',
     ]
	message('Will run on install:')
	message(pip_cmd)
	meson.add_install_script(pip_cmd)
endif
