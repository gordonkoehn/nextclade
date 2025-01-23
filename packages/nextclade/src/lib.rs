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

use crate::alphabet::nuc::to_nuc_seq;
use crate::gene::gene_map;

use pyo3::prelude::*;

/// Formats the sum of two numbers as string.
#[pyfunction]
fn translate_aa_align(ref_seq: &str, gene_ref: &str) -> PyResult<String> {
  let ref_seq = to_nuc_seq(ref_seq);
  // get gene_map - try to get via the run_args
  let gene_map = match gene_map::GeneMap::from_path(gene_ref) {
    Ok(map) => map,
    Err(e) => {
      return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
        "Error loading gene map: {}",
        e
      )))
    }
  };
  /*
   let ref_translation =
     translate_genes_ref(&ref_seq, &gene_map, &params.alignment).wrap_err("When translating reference sequence")?;
  */
  let a = 1;
  let b = 2;
  Ok((a + b).to_string())
  // need query sequence
  // reference sequence

  // reference translation
  // gene map
  // global coordinates map
  // alignment range
  // gap open close aa
  // alignment params

  /*   let translation = translate_genes(
    &alignment.qry_seq,
    &alignment.ref_seq,
    ref_translation,
    gene_map,
    &coord_map_global,
    &alignment_range,
    gap_open_close_aa,
    &params.alignment,
  )?; */
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn nextclade(_py: Python, m: &PyModule) -> PyResult<()> {
  m.add_function(wrap_pyfunction!(translate_aa_align, m)?)?;
  Ok(())
}
