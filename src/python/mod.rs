use feos_core::python::joback::PyJobackRecord;
use feos_core::python::parameter::*;
use feos_core::python::*;
use pyo3::prelude::*;
use pyo3::wrap_pymodule;
use quantity::python::PyInit_quantity;

mod eos;
pub use eos::*;
mod dft;
use dft::*;
mod parameters;
pub use parameters::*;

#[pymodule]
pub fn feos_pcsaft(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyIdentifier>()?;
    m.add_class::<PyVerbosity>()?;
    m.add_class::<PyContributions>()?;
    m.add_class::<PyChemicalRecord>()?;
    m.add_class::<PyJobackRecord>()?;

    m.add_class::<PyPcSaftRecord>()?;
    m.add_class::<PyPureRecord>()?;
    m.add_class::<PySegmentRecord>()?;
    m.add_class::<PyBinaryRecord>()?;
    m.add_class::<PyBinarySegmentRecord>()?;
    m.add_class::<PyPcSaftParameters>()?;

    m.add_wrapped(wrap_pymodule!(eos))?;
    m.add_wrapped(wrap_pymodule!(dft))?;
    m.add_wrapped(wrap_pymodule!(quantity))?;

    py.run(
        "\
import sys
sys.modules['feos_pcsaft.eos'] = eos
sys.modules['feos_pcsaft.eos.utils'] = eos.utils
sys.modules['feos_pcsaft.dft'] = dft
sys.modules['feos_pcsaft.dft.utils'] = dft.utils
quantity.SINumber.__module__ = 'feos_pcsaft.si'
quantity.SIArray1.__module__ = 'feos_pcsaft.si'
quantity.SIArray2.__module__ = 'feos_pcsaft.si'
quantity.SIArray3.__module__ = 'feos_pcsaft.si'
quantity.SIArray4.__module__ = 'feos_pcsaft.si'
sys.modules['feos_pcsaft.si'] = quantity
    ",
        None,
        Some(m.dict()),
    )?;
    Ok(())
}
