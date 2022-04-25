# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2022-04-25
### Added
- `PcSaftFunctional` now always uses `Joback` as ideal gas model if parameters are available. [#25](https://github.com/feos-org/feos-pcsaft/pull/25)
- Added entropy scaling correlation parameters of viscosity and thermal conductivity for the homo-segmented group contribution method in `PcSaftRecord::from_segments`. [#30](https://github.com/feos-org/feos-pcsaft/pull/30)
- Added parameters published by Loetgering-Lin et al. for entropy scaling of viscosity via homo GC method. [#30](https://github.com/feos-org/feos-pcsaft/pull/30)
- Added `pyproject.toml`. [#31]

### Changed
- Added optional arguments to the constructor of `PcSaftFunctional` in Python to make it more analogous to `PcSaft`. [#34](https://github.com/feos-org/feos-pcsaft/pull/34)
- Building Pc-SAFT parameters from segments does not check anymore, whether multiple polar or associating groups are present. [#33](https://github.com/feos-org/feos-pcsaft/pull/33)
- Moved the creation of the python module to the build_wheel auxilliary crate, so that only the relevant structs and macros are available for the dependents. [#29](https://github.com/feos-org/feos-pcsaft/pull/29)

### Packaging
- updated `feos-core` and `feos-dft` to `0.2`. [#30](https://github.com/feos-org/feos-pcsaft/pull/30)

## [0.1.0] - 2022-01-12
### Added
- Initial release
