use super::parameters::PyPcSaftParameters;
use crate::eos::{PcSaft, PcSaftOptions};
use feos_core::python::{PyContributions, PyVerbosity};
use feos_core::utils::{
    DataSet, EquilibriumLiquidDensity, Estimator, LiquidDensity, VaporPressure,
};
use feos_core::*;
use numpy::convert::ToPyArray;
use numpy::{PyArray1, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use quantity::python::*;
use quantity::si::*;
use std::collections::HashMap;
use std::rc::Rc;

/// Initialize PC-SAFT equation of state.
///
/// Parameters
/// ----------
/// parameters : PcSaftParameters
///     The parameters of the PC-Saft equation of state to use.
/// max_eta : float, optional
///     Maximum packing fraction. Defaults to 0.5.
/// max_iter_cross_assoc : unsigned integer, optional
///     Maximum number of iterations for cross association. Defaults to 50.
/// tol_cross_assoc : float
///     Tolerance for convergence of cross association. Defaults to 1e-10.
///
/// Returns
/// -------
/// PcSaft
///     The PC-SAFT equation of state that can be used to compute thermodynamic
///     states.
#[pyclass(name = "PcSaft", unsendable)]
#[pyo3(
    text_signature = "(parameters, max_eta, max_iter_cross_assoc, tol_cross_assoc, joback_parameters)"
)]
#[derive(Clone)]
pub struct PyPcSaft(pub Rc<PcSaft>);

#[pymethods]
impl PyPcSaft {
    #[new]
    #[args(
        max_eta = "0.5",
        max_iter_cross_assoc = "50",
        tol_cross_assoc = "1e-10"
    )]
    fn new(
        parameters: PyPcSaftParameters,
        max_eta: f64,
        max_iter_cross_assoc: usize,
        tol_cross_assoc: f64,
    ) -> Self {
        let options = PcSaftOptions {
            max_eta,
            max_iter_cross_assoc,
            tol_cross_assoc,
        };
        Self(Rc::new(PcSaft::with_options(parameters.0.clone(), options)))
    }
}

impl_equation_of_state!(PyPcSaft);
impl_virial_coefficients!(PyPcSaft);

impl_state!(PcSaft, PyPcSaft);
impl_state_molarweight!(PcSaft, PyPcSaft);
impl_state_entropy_scaling!(PcSaft, PyPcSaft);
impl_vle_state!(PcSaft, PyPcSaft);
impl_estimator!(PcSaft, PyPcSaft);

#[pymodule]
pub fn eos(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPcSaft>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseDiagramPure>()?;
    m.add_class::<PyPhaseDiagramBinary>()?;
    m.add_class::<PyPhaseDiagramHetero>()?;
    m.add_class::<PyPhaseEquilibrium>()?;

    let utils = PyModule::new(py, "utils")?;
    utils.add_class::<PyDataSet>()?;
    utils.add_class::<PyEstimator>()?;
    m.add_submodule(utils)?;
    Ok(())
}
