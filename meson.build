project(
  'libgauss2dfit',
  'cpp',
  version : '0.1.0',
  license : 'GPL-3.0-or-later',
  default_options : ['cpp_std=c++17',],
)

gauss2d_dep = dependency('gauss2d')
gsl_dep = dependency('GSL')
parameters_dep = dependency('parameters')

inc = include_directories('include')

depends = declare_dependency(
  include_directories: inc,
  dependencies: [gauss2d_dep, gsl_dep, parameters_dep],
)

subdir('include')
subdir('src')
subdir('tests')

pkg_mod = import('pkgconfig')
pkg_mod.generate(libraries : gauss2dfitlib,
                 version : '0.1',
                 name : 'libgauss2dfit',
                 filebase : 'gauss2dfit',
                 description : 'Create, manipulate and evaluate 2D Gaussian mixtures and images thereof.')
