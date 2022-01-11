use feos_pcsaft::python::feos_pcsaft;
use pyo3::prelude::*;

#[pymodule]
pub fn build_wheel(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    feos_pcsaft(py, m)
}
