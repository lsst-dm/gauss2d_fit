# Change Log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](https://semver.org/).

### [0.1.4] 2024-07-30

* Fixed: Put correct version in pkg-config file for standalone C++ builds
* Added: -ffp-contract=off to default C++ compiler options
* Fixed: Allow tiny (1e-25) deviations from equality in some Python tests
* Fixed: Python formatting
* Added: Black and ruff configurations to pyproject.toml
* Fixed: Removed incorrect subdir for standalone Python build
* See [DM-45473](https://rubinobs.atlassian.net/browse/DM-45473) for details.

### [0.1.3] 2024-07-29

* Fixed: test_channel is no longer sensitive to test order.
* Fixed: test_model no longer calls compute_hessian with overly verbose print=True.
* See [DM-45392](https://rubinobs.atlassian.net/browse/DM-45392) for details.

### [0.1.2] 2024-07-23

* Added: Shell script for pytest to handle MacOS SIP stripping DYLD_LIBRARY_PATH
* Fixed: Added missing LSST_LIBRARY_PATH entry in ups/gauss2d_fit.table
* Fixed: Added pytest.ini_options to python/pyproject.toml
* Fixed: Replaced another pow10 that was missed in previous update
* Removed: lsst.gauss2d from required Python modules in python/meson.build
* This is a hotfix for [DM-45346](https://rubinobs.atlassian.net/browse/DM-45346).

### [0.1.1] 2024-07-22

* Added: doctest upgrade to 2.4.11
* Fixed: Replace pow10(x) with pow(10.0, x) in JanskyToABMagTransform
* Fixed: Added small tolerance to Python loglike checks expecting zero (to be investigated
  in [DM-45308](https://rubinobs.atlassian.net/browse/DM-45308))
* Fixed: lib64 refs in python/meson.build and meson.build replaced with lib.
* Fixed: --libdir "lib" added to eups meson config command.
* See [DM-45346](https://rubinobs.atlassian.net/browse/DM-45346) for details.

### [0.1.0] 2024-07-09

* Changed: Initial release, forked to https://github.com/lsst/gauss2d_fit.
* See [DM-43906](https://rubinobs.atlassian.net/browse/DM-43906) for details.

[0.1.4]: https://github.com/lsst-dm/gauss2d_fit/compare/0.1.3...0.1.4

[0.1.3]: https://github.com/lsst-dm/gauss2d_fit/compare/0.1.2...0.1.3

[0.1.2]: https://github.com/lsst-dm/gauss2d_fit/compare/0.1.1...0.1.2

[0.1.1]: https://github.com/lsst-dm/gauss2d_fit/compare/0.1.0...0.1.1

[0.1.0]: https://github.com/lsst-dm/gauss2d_fit/compare/a42ec007c...0.1.0
