#![warn(clippy::all)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::suspicious_operation_groupings)]

mod dft;
mod eos;
mod parameters;

pub use dft::PcSaftFunctional;
pub use eos::{PcSaft, PcSaftOptions};
pub use parameters::{PcSaftParameters, PcSaftRecord};

#[cfg(feature = "python")]
pub mod python;
