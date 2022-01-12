# FeOs - PC-SAFT

[![crate](https://img.shields.io/crates/v/feos-pcsaft.svg)](https://crates.io/crates/feos-pcsaft)
[![documentation](https://docs.rs/feos-pcsaft/badge.svg)](https://docs.rs/feos-pcsaft)
[![documentation](https://img.shields.io/badge/docs-github--pages-blue)](https://feos-org.github.io/feos/)

Implementation of the PC(P)-SAFT equation of state[^gross2001][^gross2005][^gross2006] and corresponding Helmholtz energy functional[^sauer2016] within the FeOs project. This project contains a Rust implementation as well as bindings to Python.

## Usage in Python

If you want to use `feos-pcsaft` in Python, take a look at the [`feos`-repository](https://github.com/feos-org/feos). `feos` contains multiple equation of state implementations in a single, easy-to-use Python package.

## FeOs

> FeOs is a framework for equations of state and classical density function theory

You can learn more about the principles behind `FeOs` [here](https://feos-org.github.io/feos/).


## Parameters
[parameters](parameters/README.md)

[^gross2001]: [J. Gross and G. Sadowski (2001). *Industrial & Engineering Chemistry Research*, 40(4), 1244-1260.](https://doi.org/10.1021/ie0003887)
[^gross2005]: [J. Gross (2005). *AIChE Journal*, 51(9), 2556-2568.](https://doi.org/10.1002/aic.10502)
[^gross2006]: [J. Gross and J. Vrabec (2006). *AIChE Journal*, 52(3), 1194-1204.](https://doi.org/10.1002/aic.10683)
[^sauer2016]: [E. Sauer and J. Gross (2017). *Industrial & Engineering Chemistry Research*, 56(14), 4119-4135](https://doi.org/10.1021/acs.iecr.6b04551)

## Installation

Add this to your `Cargo.toml`

```toml
[dependencies]
feos-pcsaft = "0.1"
```

## Test building python wheel

From within a Python virtual environment with `maturin` installed, type

```
maturin build --release --out dist --no-sdist -m build_wheel/Cargo.toml
```
