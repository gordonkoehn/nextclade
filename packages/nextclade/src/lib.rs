pub mod align;
pub mod alphabet;
pub mod analyze;
pub mod constants;
pub mod coord;
pub mod features;
pub mod gene;
pub mod graph;
pub mod io;
pub mod qc;
pub mod run;
pub mod sort;
pub mod translate;
pub mod tree;
pub mod types;
pub mod utils;

#[cfg(test)]
mod tests {
  use crate::utils::global_init::global_init;
  use ctor::ctor;

  #[ctor]
  fn init() {
    global_init();
  }
}

use std::any::Any;
use std::fmt::Debug;

use crate::align::gap_open::get_gap_open_close_scores_flat;
use crate::align::insertions_strip::get_aa_insertions;
use crate::align::params::AlignPairwiseParams;
use crate::alphabet::nuc::to_nuc_seq;
use crate::analyze::aa_changes_find::aa_changes_find;
use crate::analyze::aa_changes_find_for_cds::{AaChangesParams, FindAaChangesOutput};
use crate::analyze::nuc_alignment::NucAlignment;
use crate::analyze::nuc_changes::{find_nuc_changes, FindNucChangesOutput};
use crate::coord::coord_map_global::CoordMapGlobal;
use crate::coord::range::NucRefGlobalRange;
use crate::gene::gene_map;
use crate::translate::frame_shifts_flatten::frame_shifts_flatten;
use crate::translate::translate_genes::{translate_genes, Translation};
use crate::translate::translate_genes_ref;

use std::time::Instant;

use serde::{Deserialize, Serialize};

use pyo3::prelude::*;

#[derive(Serialize, Deserialize, Debug)]
pub struct FullAligmentWithId {
  pub read_id: String,
  pub qry_seq_unaligned: Vec<alphabet::nuc::Nuc>,
  pub qry_seq: Vec<alphabet::nuc::Nuc>,
  pub translation: Translation,
  pub aa_insertions: Vec<align::insertions_strip::AaIns>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct AaAlignment {
  pub qry_seq: Vec<alphabet::nuc::Nuc>,
  pub translation: Translation,
  pub aa_insertions: Vec<align::insertions_strip::AaIns>,
}

fn perform_translation(ref_seq: &str, qry_seq: &str, gene_ref: &str) -> Result<AaAlignment, String> {
  let ref_seq = to_nuc_seq(ref_seq).map_err(|e| format!("Error converting to nucleotide sequence: {}", e))?;
  let qry_seq = to_nuc_seq(qry_seq).map_err(|e| format!("Error converting to nucleotide sequence: {}", e))?;

  let gene_map = gene_map::GeneMap::from_path(gene_ref).map_err(|e| format!("Error loading gene map: {}", e))?;
  let params_alignment = AlignPairwiseParams::default();

  let ref_translation = translate_genes_ref::translate_genes_ref(&ref_seq, &gene_map, &params_alignment)
    .map_err(|e| format!("Error translating genes: {}", e))?;

  let coord_map_global = CoordMapGlobal::new(&ref_seq);

  let FindNucChangesOutput {
    substitutions,
    deletions,
    alignment_range,
  } = find_nuc_changes(&qry_seq, &ref_seq);

  let aln = NucAlignment::new(&ref_seq, &qry_seq, &alignment_range);

  let gap_open_close_aa = get_gap_open_close_scores_flat(&ref_seq, &params_alignment);

  let translation = translate_genes(
    &qry_seq,
    &ref_seq,
    &ref_translation,
    &gene_map,
    &coord_map_global,
    &alignment_range,
    &gap_open_close_aa,
    &params_alignment,
  )
  .map_err(|e| format!("Error translating genes: {}", e))?;

  let aa_insertions = get_aa_insertions(&translation);

  Ok(AaAlignment {
    qry_seq,
    translation,
    aa_insertions,
  })
}

#[pyfunction]
fn translate_aa_align(ref_seq: &str, qry_seq: &str, gene_ref: &str) -> PyResult<String> {
  let aa_alignment =
    perform_translation(ref_seq, qry_seq, gene_ref).map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))?;

  Ok(format!(
    "TRANSLATION:\n{:?}\n\nAA INSERTIONS:\n{:?}",
    aa_alignment.translation, aa_alignment.aa_insertions
  ))
}

fn perform_full_translation_with_id(
  read_id: String,
  ref_seq: &str,
  qry_seq_unaligned: &str,
  qry_seq: &str,
  gene_ref: &str,
) -> Result<FullAligmentWithId, String> {
  let aa_alignment = perform_translation(ref_seq, qry_seq, gene_ref)?;

  let qry_seq_unaligned = to_nuc_seq(qry_seq_unaligned).map_err(|e| {
    format!(
      "Error converting unaligned query sequence to nucleotide sequence: {}",
      e
    )
  })?;

  Ok(FullAligmentWithId {
    read_id,
    qry_seq_unaligned,
    qry_seq: aa_alignment.qry_seq,
    translation: aa_alignment.translation,
    aa_insertions: aa_alignment.aa_insertions,
  })
}

#[pyfunction]
fn translate_full_align_with_id(
  read_id: String,
  ref_seq: &str,
  qry_seq_unaligned: &str,
  qry_seq: &str,
  gene_ref: &str,
) -> PyResult<String> {
  let full_alignment = perform_full_translation_with_id(read_id, ref_seq, qry_seq_unaligned, qry_seq, gene_ref)
    .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e))?;

  let ndjson = serde_json::to_string(&full_alignment)
    .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;

  Ok(ndjson)
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn nextclade(_py: Python, m: &PyModule) -> PyResult<()> {
  m.add_function(wrap_pyfunction!(translate_aa_align, m)?)?;
  m.add_function(wrap_pyfunction!(translate_full_align_with_id, m)?)?;
  Ok(())
}
