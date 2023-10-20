Gauss2DFit
#######

.. todo image:: https://img.shields.io/pypi/v/gauss2dfit.svg
   .. todo   :target: https://pypi.python.org/pypi/gauss2dfit

.. todo image:: https://img.shields.io/pypi/pyversions/gauss2dfit.svg
   .. todo   :target: https://pypi.python.org/pypi/gauss2dfit

*gauss2dfit* provides python bindings for the C++ libgauss2dfit library.
Python 3.8+, `pybind11 <https://github.com/pybind/pybind11>`_, and numpy are 
requirements, as numpy arrays are used to store images for tests.

With libgauss2dfit installed, these bindings can be built using meson:

meson --prefix=~/.local build && meson compile -C build && meson install -C build

... where the prefix argument is your desired installation path.

You may also need to configure pkg-config, e.g. prepend 
PKG_CONFIG_PATH=~/.local/lib64/pkgconfig to the build command for a local
and/or non-root installation.

Alternatively, it can be built for local development using poetry like so:

poetry build && python3 -m pip install .

... using the include pyproject.toml.

.. todo *gauss2dfit* is available in `PyPI <https://pypi.python.org/pypi/gauss2dfit>`_
   .. and thus can be easily installed via::

.. pip install gauss2dfit
