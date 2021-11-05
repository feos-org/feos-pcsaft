use crate::parameters::PcSaftParameters;
use feos_core::joback::Joback;
use feos_core::{
    Contributions, EntropyScaling, EosError, EosResult, EosUnit, EquationOfState, HelmholtzEnergy,
    IdealGasContribution, MolarWeight, State,
};
use ndarray::Array1;
use quantity::si::*;
use std::f64::consts::{FRAC_PI_6, PI};
use std::rc::Rc;

pub(crate) mod association;
pub(crate) mod dispersion;
pub(crate) mod hard_chain;
pub(crate) mod hard_sphere;
pub(crate) mod polar;
mod qspr;
use association::{Association, CrossAssociation};
use dispersion::Dispersion;
use hard_chain::HardChain;
use hard_sphere::HardSphere;
use polar::{DQVariants, Dipole, DipoleQuadrupole, Quadrupole};
use qspr::QSPR;

#[allow(clippy::upper_case_acronyms)]
enum IdealGasContributions {
    QSPR(QSPR),
    Joback(Joback),
}

#[derive(Copy, Clone)]
pub struct PcSaftOptions {
    pub max_eta: f64,
    pub max_iter_cross_assoc: usize,
    pub tol_cross_assoc: f64,
}

impl Default for PcSaftOptions {
    fn default() -> Self {
        Self {
            max_eta: 0.5,
            max_iter_cross_assoc: 50,
            tol_cross_assoc: 1e-10,
        }
    }
}

pub struct PcSaft {
    parameters: Rc<PcSaftParameters>,
    options: PcSaftOptions,
    contributions: Vec<Box<dyn HelmholtzEnergy>>,
    ideal_gas: IdealGasContributions,
}

impl PcSaft {
    pub fn new(parameters: PcSaftParameters) -> Self {
        Self::with_options(parameters, PcSaftOptions::default())
    }

    pub fn with_options(parameters: PcSaftParameters, options: PcSaftOptions) -> Self {
        let parameters = Rc::new(parameters);
        let mut contributions: Vec<Box<dyn HelmholtzEnergy>> = Vec::with_capacity(7);
        contributions.push(Box::new(HardSphere {
            parameters: parameters.clone(),
        }));
        contributions.push(Box::new(HardChain {
            parameters: parameters.clone(),
        }));
        contributions.push(Box::new(Dispersion {
            parameters: parameters.clone(),
        }));
        if parameters.ndipole > 0 {
            contributions.push(Box::new(Dipole {
                parameters: parameters.clone(),
            }));
        };
        if parameters.nquadpole > 0 {
            contributions.push(Box::new(Quadrupole {
                parameters: parameters.clone(),
            }));
        };
        if parameters.ndipole > 0 && parameters.nquadpole > 0 {
            contributions.push(Box::new(DipoleQuadrupole {
                parameters: parameters.clone(),
                variant: DQVariants::DQ35,
            }));
        };
        match parameters.nassoc {
            0 => (),
            1 => contributions.push(Box::new(Association {
                parameters: parameters.clone(),
            })),
            _ => contributions.push(Box::new(CrossAssociation {
                parameters: parameters.clone(),
                max_iter: options.max_iter_cross_assoc,
                tol: options.tol_cross_assoc,
            })),
        };

        let joback_records = parameters.joback_records.clone();

        Self {
            parameters: parameters.clone(),
            options,
            contributions,
            ideal_gas: joback_records.map_or(
                IdealGasContributions::QSPR(QSPR { parameters }),
                |joback_records| IdealGasContributions::Joback(Joback::new(joback_records)),
            ),
        }
    }
}

impl EquationOfState for PcSaft {
    fn components(&self) -> usize {
        self.parameters.pure_records.len()
    }

    fn subset(&self, component_list: &[usize]) -> Self {
        Self::with_options(self.parameters.subset(component_list), self.options)
        // match &self.ideal_gas {
        //     IdealGasContributions::QSPR(_) => {
        //         Self::with_options(self.parameters.subset(component_list), self.options, None)
        //     }
        //     IdealGasContributions::Joback(joback) => Self::with_options(
        //         self.parameters.subset(component_list),
        //         self.options,
        //         Some(joback.parameters.subset(component_list)),
        //     ),
        // }
    }

    fn compute_max_density(&self, moles: &Array1<f64>) -> f64 {
        self.options.max_eta * moles.sum()
            / (FRAC_PI_6 * &self.parameters.m * self.parameters.sigma.mapv(|v| v.powi(3)) * moles)
                .sum()
    }

    fn residual(&self) -> &[Box<dyn HelmholtzEnergy>] {
        &self.contributions
    }

    fn ideal_gas(&self) -> &dyn IdealGasContribution {
        match &self.ideal_gas {
            IdealGasContributions::QSPR(qspr) => qspr,
            IdealGasContributions::Joback(joback) => joback,
        }
    }
}

impl MolarWeight<SIUnit> for PcSaft {
    fn molar_weight(&self) -> SIArray1 {
        self.parameters.molarweight.clone() * GRAM / MOL
    }
}

fn omega11(t: f64) -> f64 {
    1.06036 * t.powf(-0.15610)
        + 0.19300 * (-0.47635 * t).exp()
        + 1.03587 * (-1.52996 * t).exp()
        + 1.76474 * (-3.89411 * t).exp()
}

fn omega22(t: f64) -> f64 {
    1.16145 * t.powf(-0.14874) + 0.52487 * (-0.77320 * t).exp() + 2.16178 * (-2.43787 * t).exp()
        - 6.435e-4 * t.powf(0.14874) * (18.0323 * t.powf(-0.76830) - 7.27371).sin()
}

impl EntropyScaling<SIUnit, Self> for PcSaft {
    fn viscosity_reference(&self, state: &State<SIUnit, Self>) -> EosResult<SINumber> {
        let p = &self.parameters;
        let mw = &p.molarweight;
        let ce: Array1<SINumber> = (0..self.components())
            .map(|i| {
                let tr = (state.temperature / p.epsilon_k[i] / KELVIN)
                    .into_value()
                    .unwrap();
                5.0 / 16.0
                    * (mw[i] * GRAM / MOL * KB / NAV * state.temperature / PI)
                        .sqrt()
                        .unwrap()
                    / omega22(tr)
                    / (p.sigma[i] * ANGSTROM).powi(2)
            })
            .collect();
        let x = &state.molefracs;
        let mut ce_mix = 0.0 * MILLI * PASCAL * SECOND;
        for i in 0..self.components() {
            let denom: f64 = (0..self.components())
                .map(|j| {
                    x[j] * (1.0
                        + (ce[i] / ce[j]).into_value().unwrap().sqrt()
                            * (mw[j] / mw[i]).powf(1.0 / 4.0))
                    .powi(2)
                        / (8.0 * (1.0 + mw[i] / mw[j])).sqrt()
                })
                .sum();
            ce_mix += ce[i] * x[i] / denom
        }
        Ok(ce_mix)
    }

    fn viscosity_correlation(&self, s_res: f64, x: &Array1<f64>) -> EosResult<f64> {
        let coefficients = self
            .parameters
            .viscosity
            .as_ref()
            .expect("Missing viscosity coefficients.");
        let m = (x * &self.parameters.m).sum();
        let s = s_res / m;
        let pref = (x * &self.parameters.m).mapv(|v| v / m);
        let a: f64 = (&coefficients.row(0) * x).sum();
        let b: f64 = (&coefficients.row(1) * &pref).sum();
        let c: f64 = (&coefficients.row(2) * &pref).sum();
        let d: f64 = (&coefficients.row(3) * &pref).sum();
        Ok(a + b * s + c * s.powi(2) + d * s.powi(3))
    }

    fn diffusion_reference(&self, state: &State<SIUnit, Self>) -> EosResult<SINumber> {
        if self.components() != 1 {
            return Err(EosError::IncompatibleComponents(self.components(), 1));
        }
        let p = &self.parameters;
        let res: Array1<SINumber> = (0..self.components())
            .map(|i| {
                let tr = (state.temperature / p.epsilon_k[i] / KELVIN)
                    .into_value()
                    .unwrap();
                3.0 / 8.0 / (p.sigma[i] * ANGSTROM).powi(2) / omega11(tr) / (state.density * NAV)
                    * (state.temperature * RGAS / PI / (p.molarweight[i] * GRAM / MOL))
                        .sqrt()
                        .unwrap()
            })
            .collect();
        Ok(res[0])
    }

    fn diffusion_correlation(&self, s_res: f64, x: &Array1<f64>) -> EosResult<f64> {
        if self.components() != 1 {
            return Err(EosError::IncompatibleComponents(self.components(), 1));
        }
        let coefficients = self
            .parameters
            .diffusion
            .as_ref()
            .expect("Missing diffusion coefficients.");
        let m = (x * &self.parameters.m).sum();
        let s = s_res / m;
        let pref = (x * &self.parameters.m).mapv(|v| v / m);
        let a: f64 = (&coefficients.row(0) * x).sum();
        let b: f64 = (&coefficients.row(1) * &pref).sum();
        let c: f64 = (&coefficients.row(2) * &pref).sum();
        let d: f64 = (&coefficients.row(3) * &pref).sum();
        let e: f64 = (&coefficients.row(4) * &pref).sum();
        Ok(a + b * s - c * (1.0 - s.exp()) * s.powi(2) - d * s.powi(4) - e * s.powi(8))
    }

    // fn thermal_conductivity_reference(
    //     &self,
    //     state: &State<SIUnit, E>,
    // ) -> EosResult<SINumber> {
    //     if self.components() != 1 {
    //         return Err(EosError::IncompatibleComponents(self.components(), 1));
    //     }
    //     let p = &self.parameters;
    //     let res: Array1<SINumber> = (0..self.components())
    //         .map(|i| {
    //             let tr = (state.temperature / p.epsilon_k[i] / KELVIN)
    //                 .into_value()
    //                 .unwrap();
    //             let cp = State::critical_point_pure(&state.eos, Some(state.temperature)).unwrap();
    //             let s_res_cp_reduced = cp
    //                 .entropy(Contributions::Residual)
    //                 .to_reduced(SIUnit::reference_entropy())
    //                 .unwrap();
    //             let s_res_reduced = cp
    //                 .entropy(Contributions::Residual)
    //                 .to_reduced(SIUnit::reference_entropy())
    //                 .unwrap();
    //             let ref_ce = 0.083235
    //                 * ((state.temperature / KELVIN).into_value().unwrap()
    //                     / (p.molarweight[0]
    //                     / p.m[0]))
    //                     .sqrt()
    //                 / p.sigma[0]
    //                 / p.sigma[0]
    //                 / omega22(tr)
    //                 * p.m[0];
    //             let alpha_visc = (-s_res_reduced / s_res_cp_reduced).exp();
    //             let ref_ts = (-0.0167141 * tr / p.m[0] + 0.0470581 * (tr / p.m[0]).powi(2))
    //                 * (p.m[0] * p.m[0] * p.sigma[i].powi(3) * p.epsilon_k[0])
    //                 / 100000.0;
    //             (ref_ce + ref_ts * alpha_visc) * WATT / METER / KELVIN
    //         })
    //         .collect();
    //     Ok(res[0])
    // }

    // Equation 11 of DOI: 10.1021/acs.iecr.9b03998
    fn thermal_conductivity_reference(&self, state: &State<SIUnit, Self>) -> EosResult<SINumber> {
        if self.components() != 1 {
            return Err(EosError::IncompatibleComponents(self.components(), 1));
        }
        let p = &self.parameters;
        let res: Array1<SINumber> = (0..self.components())
            .map(|i| {
                let tr = (state.temperature / p.epsilon_k[i] / KELVIN)
                    .into_value()
                    .unwrap();
                let ce = 83.235
                    * f64::powf(10.0, -1.5)
                    * ((state.temperature / KELVIN).into_value().unwrap() / p.molarweight[0]
                        * p.m[0])
                        .sqrt()
                    / (p.sigma[0] * p.sigma[0])
                    / omega22(tr);
                ce * WATT / METER / KELVIN
                    + state.density
                        * self.diffusion_reference(state).unwrap()
                        * self
                            .diffusion_correlation(
                                state
                                    .molar_entropy(Contributions::Residual)
                                    .to_reduced(SIUnit::reference_molar_entropy())
                                    .unwrap(),
                                &state.molefracs,
                            )
                            .unwrap()
                        * (state.c_v(Contributions::Total) - 1.5 * RGAS)
            })
            .collect();
        Ok(res[0])
    }

    fn thermal_conductivity_correlation(&self, s_res: f64, x: &Array1<f64>) -> EosResult<f64> {
        if self.components() != 1 {
            return Err(EosError::IncompatibleComponents(self.components(), 1));
        }
        let coefficients = self
            .parameters
            .thermal_conductivity
            .as_ref()
            .expect("Missing thermal conductivity coefficients");
        let m = (x * &self.parameters.m).sum();
        let s = s_res / m;
        let pref = (x * &self.parameters.m).mapv(|v| v / m);
        let a: f64 = (&coefficients.row(0) * x).sum();
        let b: f64 = (&coefficients.row(1) * &pref).sum();
        let c: f64 = (&coefficients.row(2) * &pref).sum();
        let d: f64 = (&coefficients.row(3) * &pref).sum();
        Ok(a + b * s + c * (1.0 - s.exp()) + d * s.powi(2))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parameters::utils::{
        butane_parameters, propane_butane_parameters, propane_parameters,
    };
    use approx::assert_relative_eq;
    use feos_core::{Contributions, DensityInitialization, PhaseEquilibrium, State};
    use ndarray::arr1;
    use quantity::si::{BAR, KELVIN, METER, PASCAL, RGAS, SECOND};

    #[test]
    fn ideal_gas_pressure() {
        let e = Rc::new(PcSaft::new(propane_parameters()));
        let t = 200.0 * KELVIN;
        let v = 1e-3 * METER.powi(3);
        let n = arr1(&[1.0]) * MOL;
        let s = State::new_nvt(&e, t, v, &n).unwrap();
        let p_ig = s.total_moles * RGAS * t / v;
        assert_relative_eq!(s.pressure(Contributions::IdealGas), p_ig, epsilon = 1e-10);
        assert_relative_eq!(
            s.pressure(Contributions::IdealGas) + s.pressure(Contributions::Residual),
            s.pressure(Contributions::Total),
            epsilon = 1e-10
        );
    }

    #[test]
    fn ideal_gas_heat_capacity_joback() {
        let e = Rc::new(PcSaft::new(propane_parameters()));
        let t = 200.0 * KELVIN;
        let v = 1e-3 * METER.powi(3);
        let n = arr1(&[1.0]) * MOL;
        let s = State::new_nvt(&e, t, v, &n).unwrap();
        let p_ig = s.total_moles * RGAS * t / v;
        assert_relative_eq!(s.pressure(Contributions::IdealGas), p_ig, epsilon = 1e-10);
        assert_relative_eq!(
            s.pressure(Contributions::IdealGas) + s.pressure(Contributions::Residual),
            s.pressure(Contributions::Total),
            epsilon = 1e-10
        );
    }

    #[test]
    fn new_tpn() {
        let e = Rc::new(PcSaft::new(propane_parameters()));
        let t = 300.0 * KELVIN;
        let p = BAR;
        let m = arr1(&[1.0]) * MOL;
        let s = State::new_npt(&e, t, p, &m, DensityInitialization::None);
        let p_calc = if let Ok(state) = s {
            state.pressure(Contributions::Total)
        } else {
            0.0 * PASCAL
        };
        assert_relative_eq!(p, p_calc, epsilon = 1e-6);
    }

    #[test]
    fn vle_pure_t() {
        let e = Rc::new(PcSaft::new(propane_parameters()));
        let t = 300.0 * KELVIN;
        let vle = PhaseEquilibrium::pure_t(&e, t, None, Default::default());
        if let Ok(v) = vle {
            assert_relative_eq!(
                v.vapor().pressure(Contributions::Total),
                v.liquid().pressure(Contributions::Total),
                epsilon = 1e-6
            )
        }
    }

    #[test]
    fn critical_point() {
        let e = Rc::new(PcSaft::new(propane_parameters()));
        let t = 300.0 * KELVIN;
        let cp = State::critical_point(&e, None, Some(t), Default::default());
        if let Ok(v) = cp {
            assert_relative_eq!(v.temperature, 375.1244078318015 * KELVIN, epsilon = 1e-8)
        }
    }

    #[test]
    fn speed_of_sound() {
        let e = Rc::new(PcSaft::new(propane_parameters()));
        let t = 300.0 * KELVIN;
        let p = BAR;
        let m = arr1(&[1.0]) * MOL;
        let s = State::new_npt(&e, t, p, &m, DensityInitialization::None).unwrap();
        assert_relative_eq!(
            s.speed_of_sound(),
            245.00185709137546 * METER / SECOND,
            epsilon = 1e-4
        )
    }

    #[test]
    fn mix_single() {
        let e1 = Rc::new(PcSaft::new(propane_parameters()));
        let e2 = Rc::new(PcSaft::new(butane_parameters()));
        let e12 = Rc::new(PcSaft::new(propane_butane_parameters()));
        let t = 300.0 * KELVIN;
        let v = 0.02456883872966545 * METER.powi(3);
        let m1 = arr1(&[2.0]) * MOL;
        let m1m = arr1(&[2.0, 0.0]) * MOL;
        let m2m = arr1(&[0.0, 2.0]) * MOL;
        let s1 = State::new_nvt(&e1, t, v, &m1).unwrap();
        let s2 = State::new_nvt(&e2, t, v, &m1).unwrap();
        let s1m = State::new_nvt(&e12, t, v, &m1m).unwrap();
        let s2m = State::new_nvt(&e12, t, v, &m2m).unwrap();
        assert_relative_eq!(
            s1.pressure(Contributions::Total),
            s1m.pressure(Contributions::Total),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            s2.pressure(Contributions::Total),
            s2m.pressure(Contributions::Total),
            epsilon = 1e-12
        )
    }

    #[test]
    fn viscosity() -> EosResult<()> {
        let e = Rc::new(PcSaft::new(propane_parameters()));
        let t = 300.0 * KELVIN;
        let p = BAR;
        let n = arr1(&[1.0]) * MOL;
        let s = State::new_npt(&e, t, p, &n, DensityInitialization::None).unwrap();
        assert_relative_eq!(
            s.viscosity()?,
            0.00797 * MILLI * PASCAL * SECOND,
            epsilon = 1e-5
        );
        assert_relative_eq!(
            s.ln_viscosity_reduced()?,
            (s.viscosity()? / e.viscosity_reference(&s)?)
                .into_value()
                .unwrap()
                .ln(),
            epsilon = 1e-15
        );
        Ok(())
    }

    #[test]
    fn diffusion() -> EosResult<()> {
        let e = Rc::new(PcSaft::new(propane_parameters()));
        let t = 300.0 * KELVIN;
        let p = BAR;
        let n = arr1(&[1.0]) * MOL;
        let s = State::new_npt(&e, t, p, &n, DensityInitialization::None).unwrap();
        assert_relative_eq!(
            s.diffusion()?,
            0.01505 * (CENTI * METER).powi(2) / SECOND,
            epsilon = 1e-5
        );
        assert_relative_eq!(
            s.ln_diffusion_reduced()?,
            (s.diffusion()? / e.diffusion_reference(&s)?)
                .into_value()
                .unwrap()
                .ln(),
            epsilon = 1e-15
        );
        Ok(())
    }
}
