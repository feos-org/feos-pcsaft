use feos_core::*;
use feos_dft::adsorption::*;
use feos_dft::fundamental_measure_theory::FMTVersion;
use feos_dft::interface::*;
use feos_dft::python::*;
use feos_dft::solvation::*;
use feos_dft::*;
use feos_pcsaft::python::*;
use feos_pcsaft::{PcSaftFunctional, PcSaftOptions};
use numpy::*;
use pyo3::exceptions::{PyIndexError, PyValueError};
use pyo3::prelude::*;
use quantity::python::*;
use quantity::si::*;
use std::collections::HashMap;
use std::rc::Rc;

/// PC-SAFT Helmholtz energy functional.
///
/// Parameters
/// ----------
/// parameters: PcSaftParameters
///     The set of PC-SAFT parameters.
/// fmt_version: FMTVersion, optional
///     The specific variant of the FMT term. Defaults to FMTVersion.WhiteBear
/// max_eta : float, optional
///     Maximum packing fraction. Defaults to 0.5.
/// max_iter_cross_assoc : unsigned integer, optional
///     Maximum number of iterations for cross association. Defaults to 50.
/// tol_cross_assoc : float
///     Tolerance for convergence of cross association. Defaults to 1e-10.
/// dq_variant : {'dq35', 'dq44'}, optional
///     Combination rule used in the dipole/quadrupole term. Defaults to 'dq35'
///
/// Returns
/// -------
/// PcSaftFunctional
#[pyclass(name = "PcSaftFunctional", unsendable)]
#[pyo3(
    text_signature = "(parameters, fmt_version, max_eta, max_iter_cross_assoc, tol_cross_assoc, dq_variant)"
)]
#[derive(Clone)]
pub struct PyPcSaftFunctional(pub Rc<DFT<PcSaftFunctional>>);

#[pymethods]
impl PyPcSaftFunctional {
    #[new]
    #[args(
        fmt_version = "FMTVersion::WhiteBear",
        max_eta = "0.5",
        max_iter_cross_assoc = "50",
        tol_cross_assoc = "1e-10",
        dq_variant = "\"dq35\""
    )]
    fn new(
        parameters: PyPcSaftParameters,
        fmt_version: FMTVersion,
        max_eta: f64,
        max_iter_cross_assoc: usize,
        tol_cross_assoc: f64,
        dq_variant: &str,
    ) -> Self {
        let options = PcSaftOptions {
            max_eta,
            max_iter_cross_assoc,
            tol_cross_assoc,
            dq_variant: dq_variant.into(),
        };
        Self(Rc::new(PcSaftFunctional::with_options(
            parameters.0,
            fmt_version,
            options,
        )))
    }
}

impl_equation_of_state!(PyPcSaftFunctional);

impl_state!(DFT<PcSaftFunctional>, PyPcSaftFunctional);
impl_state_molarweight!(DFT<PcSaftFunctional>, PyPcSaftFunctional);
impl_phase_equilibrium!(DFT<PcSaftFunctional>, PyPcSaftFunctional);

impl_planar_interface!(PcSaftFunctional);
impl_surface_tension_diagram!(PcSaftFunctional);

impl_pore!(PcSaftFunctional, PyPcSaftFunctional);
impl_adsorption!(PcSaftFunctional, PyPcSaftFunctional);

impl_pair_correlation!(PcSaftFunctional);
impl_solvation_profile!(PcSaftFunctional);

#[pymodule]
pub fn dft(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPcSaftFunctional>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseDiagram>()?;
    m.add_class::<PyPhaseEquilibrium>()?;
    m.add_class::<PyPlanarInterface>()?;
    m.add_class::<Geometry>()?;
    m.add_class::<PyPore1D>()?;
    m.add_class::<PyPore3D>()?;
    m.add_class::<PyPairCorrelation>()?;
    m.add_class::<PyExternalPotential>()?;
    m.add_class::<PyAdsorption1D>()?;
    m.add_class::<PyAdsorption3D>()?;
    m.add_class::<PySurfaceTensionDiagram>()?;
    m.add_class::<PyDFTSolver>()?;
    m.add_class::<PySolvationProfile>()?;
    m.add_class::<FMTVersion>()?;
    Ok(())
}
