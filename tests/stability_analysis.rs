use feos_core::parameter::{IdentifierOption, Parameter};
use feos_core::{DensityInitialization, PhaseEquilibrium, State};
use feos_pcsaft::{PcSaft, PcSaftParameters};
use ndarray::arr1;
use quantity::si::*;
use std::error::Error;
use std::rc::Rc;

#[test]
fn test_stability_analysis() -> Result<(), Box<dyn Error>> {
    let params = PcSaftParameters::from_json(
        vec!["water_np", "hexane"],
        "tests/test_parameters.json",
        None,
        IdentifierOption::Name,
    )?;
    let mix = Rc::new(PcSaft::new(Rc::new(params)));
    let unstable = State::new_npt(
        &mix,
        300.0 * KELVIN,
        1.0 * BAR,
        &(arr1(&[0.5, 0.5]) * MOL),
        DensityInitialization::Liquid,
    )?;
    let check = unstable.stability_analysis(Default::default())?;
    assert!(check.len() > 0);

    let params = PcSaftParameters::from_json(
        vec!["propane", "butane"],
        "tests/test_parameters.json",
        None,
        IdentifierOption::Name,
    )?;
    let mix = Rc::new(PcSaft::new(Rc::new(params)));
    let vle = PhaseEquilibrium::bubble_point_tx(
        &mix,
        300.0 * KELVIN,
        Some(6.0 * BAR),
        &arr1(&[0.5, 0.5]),
        None,
        Default::default(),
    )?;
    let vapor_check = vle.vapor().stability_analysis(Default::default())?;
    let liquid_check = vle.liquid().stability_analysis(Default::default())?;
    assert!(vapor_check.len() == 0);
    assert!(liquid_check.len() == 0);
    Ok(())
}
