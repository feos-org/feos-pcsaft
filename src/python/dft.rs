use super::parameters::*;
use crate::dft::PcSaftFunctional;
use feos_core::python::{PyContributions, PyVerbosity};
use feos_core::utils::{
    DataSet, EquilibriumLiquidDensity, Estimator, LiquidDensity, VaporPressure,
};
use feos_core::*;
use feos_dft::adsorption::*;
use feos_dft::interface::*;
use feos_dft::python::*;
use feos_dft::solvation::*;
use feos_dft::*;
use numpy::*;
use pyo3::exceptions::PyValueError;
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
///
/// Returns
/// -------
/// PcSaftFunctional
#[pyclass(name = "PcSaftFunctional", unsendable)]
#[pyo3(text_signature = "(parameters)")]
#[derive(Clone)]
pub struct PyPcSaftFunctional(pub Rc<DFT<PcSaftFunctional>>);

#[pymethods]
impl PyPcSaftFunctional {
    #[new]
    fn new(parameters: PyPcSaftParameters) -> Self {
        Self(Rc::new(PcSaftFunctional::new(parameters.0.clone())))
    }

    /// PCP SAFT Helmholtz energy functional without simplifications
    /// for pure components.
    ///
    /// Parameters
    /// ----------
    /// parameters: PcSaftParameters
    ///     The set of SAFT parameters.
    /// fmt_version: FMTVersion
    ///     Specify the FMT term.
    ///
    /// Returns
    /// -------
    /// PcSaftFunctional
    #[staticmethod]
    #[pyo3(text_signature = "(parameters, fmt_version)")]
    fn new_full(parameters: PyPcSaftParameters, fmt_version: PyFMTVersion) -> Self {
        Self(Rc::new(PcSaftFunctional::new_full(
            parameters.0.clone(),
            fmt_version.0.clone(),
        )))
    }
}

impl_equation_of_state!(PyPcSaftFunctional);

impl_state!(DFT<PcSaftFunctional>, PyPcSaftFunctional);
impl_state_molarweight!(DFT<PcSaftFunctional>, PyPcSaftFunctional);
impl_vle_state!(DFT<PcSaftFunctional>, PyPcSaftFunctional);

impl_estimator!(DFT<PcSaftFunctional>, PyPcSaftFunctional);

impl_planar_interface!(PcSaftFunctional);
impl_surface_tension_diagram!(PcSaftFunctional);

impl_pore!(PcSaftFunctional, PyPcSaftFunctional);
impl_adsorption!(PcSaftFunctional, PyPcSaftFunctional);

impl_pair_correlation!(PcSaftFunctional);
impl_solvation_profile!(PcSaftFunctional);

#[pymodule]
pub fn dft(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyPcSaftFunctional>()?;
    m.add_class::<PyState>()?;
    m.add_class::<PyPhaseDiagramPure>()?;
    m.add_class::<PyPhaseDiagramBinary>()?;
    m.add_class::<PyPhaseDiagramHetero>()?;
    m.add_class::<PyPhaseEquilibrium>()?;
    m.add_class::<PyPlanarInterface>()?;
    m.add_class::<PyGeometry>()?;
    m.add_class::<PyPore1D>()?;
    m.add_class::<PyPore3D>()?;
    m.add_class::<PyPairCorrelation>()?;
    m.add_class::<PyExternalPotential>()?;
    m.add_class::<PyAdsorption1D>()?;
    m.add_class::<PyAdsorption3D>()?;
    m.add_class::<PySurfaceTensionDiagram>()?;
    m.add_class::<PyDFTSolver>()?;
    m.add_class::<PySolvationProfile>()?;
    m.add_class::<PyFMTVersion>()?;

    let utils = PyModule::new(py, "utils")?;
    utils.add_class::<PyDataSet>()?;
    utils.add_class::<PyEstimator>()?;
    m.add_submodule(utils)?;
    Ok(())
}
