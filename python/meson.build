project(
    'libgauss2dfit_python',
    'cpp',
    default_options : ['cpp_std=c++17',],
    license : 'GPL-3.0-or-later',
    version : run_command('cat', meson.source_root() + '/VERSION', check: true).stdout().strip(),
)

add_project_arguments('-DVERSION=@0@'.format(meson.project_version()), language : 'cpp')

eups = find_program('eups', required: false)
use_eups = eups.found()

if use_eups
	gauss2dfit_eupsdir = run_command('sh', '-c', 'echo $GAUSS2D_FIT_DIR', check: false, capture: true)
	use_eups = gauss2dfit_eupsdir.returncode() == 0
endif

gsl_dep = dependency('gsl', required: false)

if gsl_dep.found()
    add_project_arguments('-DLSST_GAUSS2D_FIT_HAS_GSL', language : 'cpp')
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
    gauss2dfit_eupsdir = gauss2dfit_eups.stdout().strip() / 'build-release/'
    parameters_eupsdir = run_command(
        'sh', '-c', 'echo $MODELFIT_PARAMETERS_DIR', check: true
    ).stdout().strip() / 'build-release'
    # TODO: Can this dependency be extracted from '../meson.build' somehow
    gauss2dfit_dep = declare_dependency(
        dependencies: [cpp.find_library('gauss2dfit', dirs: [gauss2dfit_eupsdir / 'lib64']),],
        include_directories: [
            gauss2dfit_eupsdir / 'include',
            parameters_eupsdir / 'include',
        ]
    )
else
    gauss2d_dep = dependency('lsst_gauss2d')
    gauss2d_python_dep = dependency('lsst_gauss2d_python')
    gauss2dfit_dep = dependency('lsst_gauss2d_fit')
endif

pybind11_dep = dependency('pybind11')
pymod = import('python')
python = pymod.find_installation(
    'python3',
    required : true,
    modules: ['lsst.gauss2d', 'pytest'],
)
subdir('src')

if not use_eups
    # A hacky way to pip install (to say the least) the Python module
    # this is "necessary" to allow its use as a dependency in other packages
    # (e.g. in requirements.txt, pyproject.toml)
    # because meson doesn't generate a dist-info
    # the downside is there are now two ways to manage/install the package
    # pip doesn't seem to find the pybind extension module automagically either
    # (perhaps some config is missing in pyproject.toml)
    subdir_module = meson.source_root() / 'src/lsst/gauss2d/fit'
	pip_cmd = [
	    'python', '-m', 'pip', 'install', meson.source_root(),
	     '-t', python.get_install_dir(),
	     '--no-dependencies',
	     '--no-build-isolation',
	     '--force-reinstall',
     ]
	message('Will run on install:')
	message(pip_cmd)
	add_link_cmd = ['ln', '-fs', module.full_path(), subdir_module]
	meson.add_install_script(add_link_cmd)
	meson.add_install_script(pip_cmd)
	remove_link_cmd = ['rm', subdir_module / module.full_path().split('/')[-1]]
	meson.add_install_script(remove_link_cmd)
endif

# Run tests - not sure why workdir is required, but it is
test('pytest', python, args : ['-m', 'pytest', 'tests'], workdir : meson.project_source_root())