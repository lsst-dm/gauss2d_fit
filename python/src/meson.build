subdir_gauss2d_fit = 'lsst/gauss2d/fit'

sources = [
    subdir_gauss2d_fit / 'centroidparameters.cc',
    subdir_gauss2d_fit / 'channel.cc',
    subdir_gauss2d_fit / 'chromatic.cc',
    subdir_gauss2d_fit / 'component.cc',
    subdir_gauss2d_fit / 'componentmixture.cc',
    subdir_gauss2d_fit / 'data.cc',
    subdir_gauss2d_fit / 'ellipticalcomponent.cc',
    subdir_gauss2d_fit / 'fractionalintegralmodel.cc',
    subdir_gauss2d_fit / 'gaussiancomponent.cc',
    subdir_gauss2d_fit / 'gaussianmodelintegral.cc',
    subdir_gauss2d_fit / 'gaussianparametricellipse.cc',
    subdir_gauss2d_fit / 'gaussianprior.cc',
    subdir_gauss2d_fit / 'gsl.cc',
    subdir_gauss2d_fit / 'gslsersicmixinterpolator.cc',
    subdir_gauss2d_fit / 'integralmodel.cc',
    subdir_gauss2d_fit / 'interpolation.cc',
    subdir_gauss2d_fit / 'linearintegralmodel.cc',
    subdir_gauss2d_fit / 'linearsersicmixinterpolator.cc',
    subdir_gauss2d_fit / 'model.cc',
    subdir_gauss2d_fit / 'observation.cc',
    subdir_gauss2d_fit / 'param_filter.cc',
    subdir_gauss2d_fit / 'parameters.cc',
    subdir_gauss2d_fit / 'parametric.cc',
    subdir_gauss2d_fit / 'parametricgaussian1d.cc',
    subdir_gauss2d_fit / 'parametricellipse.cc',
    subdir_gauss2d_fit / 'parametricmodel.cc',
    subdir_gauss2d_fit / 'prior.cc',
    subdir_gauss2d_fit / 'psfmodel.cc',
    subdir_gauss2d_fit / 'sersicmix.cc',
    subdir_gauss2d_fit / 'sersicmixcomponent.cc',
    subdir_gauss2d_fit / 'sersicparametricellipse.cc',
    subdir_gauss2d_fit / 'shapeprior.cc',
    subdir_gauss2d_fit / 'source.cc',
    subdir_gauss2d_fit / 'pybind11.cc',
]

kwargs = {}

eups = find_program('eups', required: false)
use_eups = eups.found()
if use_eups
	use_eups = run_command('eups', 'list', '-d', 'gauss2d_fit', check: false).returncode() == 0
endif

if use_eups
    # meson.current_build_dir() would be ideal but it seems to delete the file?
    kwargs = {'install_dir': meson.build_root() / 'lib' / subdir_gauss2d_fit}
else
    # Folder relative to site-packages to install to
    kwargs = {'subdir': subdir_gauss2d_fit}
endif

module = python.extension_module(
    '_gauss2d_fit', # Name of the module
    sources,  # List of sources
    install: true,  # Install it
    dependencies: [ # List of dependencies
        python.dependency(), # Add Python.h as dependency
        gauss2d_dep,
        gauss2d_python_dep,
        gauss2dfit_dep,
        pybind11_dep,
    ],
    kwargs: kwargs
)

# Pure Python sources
python_sources = [
    subdir_gauss2d_fit / '__init__.py',
    subdir_gauss2d_fit / 'utils.py',
]

if use_eups
    install_data(python_sources, kwargs: kwargs)
else
    # Install pure Python
    python.install_sources(
        python_sources,
        pure: false, # Will be installed next to binaries
        subdir: 'lsst',  # Folder relative to site-packages to install to
        preserve_path: true,
    )
endif
