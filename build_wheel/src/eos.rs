use feos_core::*;
use feos_pcsaft::python::PyPcSaftParameters;
use feos_pcsaft::{PcSaft, PcSaftOptions};
use numpy::convert::ToPyArray;
use numpy::{PyArray1, PyArray2};
use pyo3::exceptions::{PyIndexError, PyValueError};
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
/// dq_variant : {'dq35', 'dq44'}, optional
///     Combination rule used in the dipole/quadrupole term. Defaults to 'dq35'
///
/// Returns
/// -------
/// PcSaft
///     The PC-SAFT equation of state that can be used to compute thermodynamic
///     states.
#[pyclass(name = "PcSaft", unsendable)]
#[pyo3(text_signature = "(parameters, max_eta, max_iter_cross_assoc, tol_cross_assoc, dq_variant)")]
#[derive(Clone)]
pub struct PyPcSaft(pub Rc<PcSaft>);

#[pymethods]
impl PyPcSaft {
    #[new]
    #[args(
        max_eta = "0.5",
        max_iter_cross_assoc = "50",
        tol_cross_assoc = "1e-10",
        dq_variant = "\"dq35\""
    )]
    fn new(
        parameters: PyPcSaftParameters,
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
        Self(Rc::new(PcSaft::with_options(parameters.0, options)))
    }
}

impl_equation_of_state!(PyPcSaft);
impl_virial_coefficients!(PyPcSaft);

impl_state!(PcSaft, PyPcSaft);
impl_state_molarweight!(PcSaft, PyPcSaft);
impl_state_entropy_scaling!(PcSaft, PyPcSaft);
impl_phase_equilibrium!(PcSaft, PyPcSaft);

#[pymodule]
pub fn eos(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPcSaft>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseDiagram>()?;
    m.add_class::<PyPhaseEquilibrium>()?;
    Ok(())
}
