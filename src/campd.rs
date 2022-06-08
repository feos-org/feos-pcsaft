use crate::eos::PcSaft;
use crate::parameters::{PcSaftParameters, PcSaftRecord};
use feos_campd::{PropertyModel, SegmentAndBondCount};
use feos_core::joback::JobackRecord;
use feos_core::parameter::{Parameter, ParameterError, PureRecord, SegmentRecord};
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::rc::Rc;

#[derive(Clone, Serialize, Deserialize)]
pub struct PcSaftPropertyModel(Vec<SegmentRecord<PcSaftRecord, JobackRecord>>);

impl PcSaftPropertyModel {
    pub fn new<P: AsRef<Path>>(file: P) -> Result<Self, ParameterError> {
        Ok(Self(SegmentRecord::from_json(file)?))
    }
}

impl PropertyModel<SegmentAndBondCount> for PcSaftPropertyModel {
    type Eos = PcSaft;

    fn build_eos(&self, chemical_record: SegmentAndBondCount) -> Result<PcSaft, ParameterError> {
        Ok(PcSaft::new(Rc::new(PcSaftParameters::from_segments(
            vec![chemical_record],
            self.0.clone(),
            None,
        )?)))
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct PcSaftFixedPropertyModel(PureRecord<PcSaftRecord, JobackRecord>);

impl PcSaftFixedPropertyModel {
    pub fn new(pure_record: PureRecord<PcSaftRecord, JobackRecord>) -> Self {
        Self(pure_record)
    }
}

impl PropertyModel<()> for PcSaftFixedPropertyModel {
    type Eos = PcSaft;

    fn build_eos(&self, _: ()) -> Result<PcSaft, ParameterError> {
        Ok(PcSaft::new(Rc::new(PcSaftParameters::new_pure(
            self.0.clone(),
        ))))
    }
}
