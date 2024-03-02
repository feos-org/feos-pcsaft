#![warn(clippy::all)]
#![allow(clippy::too_many_arguments)]

#[cfg(feature = "campd")]
mod campd;
#[cfg(feature = "dft")]
mod dft;
mod eos;
mod parameters;

#[cfg(feature = "campd")]
pub use campd::{PcSaftFixedPropertyModel, PcSaftPropertyModel};
#[cfg(feature = "dft")]
pub use dft::PcSaftFunctional;
pub use eos::{PcSaft, PcSaftOptions};
pub use parameters::{PcSaftParameters, PcSaftRecord};

#[cfg(feature = "python")]
pub mod python;
