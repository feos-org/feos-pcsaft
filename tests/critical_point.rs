use approx::assert_relative_eq;
use feos_core::parameter::{IdentifierOption, Parameter};
use feos_core::State;
use feos_pcsaft::{PcSaft, PcSaftParameters};
use ndarray::arr1;
use quantity::si::*;
use std::error::Error;
use std::rc::Rc;

#[test]
fn test_critical_point_pure() -> Result<(), Box<dyn Error>> {
    let params = PcSaftParameters::from_json(
        &["propane"],
        "tests/test_parameters.json",
        None,
        IdentifierOption::Name,
    )?;
    let saft = Rc::new(PcSaft::new(params));
    let t = 300.0 * KELVIN;
    let cp = State::critical_point(&saft, None, Some(t), Default::default())?;
    assert_relative_eq!(cp.temperature, 375.12441 * KELVIN, max_relative = 1e-8);
    assert_relative_eq!(
        cp.density,
        4733.00377 * MOL / METER.powi(3),
        max_relative = 1e-6
    );
    Ok(())
}

#[test]
fn test_critical_point_mix() -> Result<(), Box<dyn Error>> {
    let params = PcSaftParameters::from_json(
        &["propane", "butane"],
        "tests/test_parameters.json",
        None,
        IdentifierOption::Name,
    )?;
    let saft = Rc::new(PcSaft::new(params));
    let t = 300.0 * KELVIN;
    let moles = arr1(&[1.5, 1.5]) * MOL;
    let cp = State::critical_point(&saft, Some(&moles), Some(t), Default::default())?;
    assert_relative_eq!(cp.temperature, 407.93481 * KELVIN, max_relative = 1e-8);
    assert_relative_eq!(
        cp.density,
        4265.50745 * MOL / METER.powi(3),
        max_relative = 1e-6
    );
    Ok(())
}
