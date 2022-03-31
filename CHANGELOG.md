# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- `PcSaftFunctional` now always uses `Joback` as ideal gas model if parameters are available. [#25](https://github.com/feos-org/feos-pcsaft/pull/25)

### Changed
- Added optional arguments to the constructor of `PcSaftFunctional` in Python to make it more analogous to `PcSaft`. [#34](https://github.com/feos-org/feos-pcsaft/pull/34)
- Building Pc-SAFT parameters from segments does not check anymore, whether multiple polar or associating groups are present. [#33](https://github.com/feos-org/feos-pcsaft/pull/33)

## [0.1.0] - 2022-01-12
### Added
- Initial release
