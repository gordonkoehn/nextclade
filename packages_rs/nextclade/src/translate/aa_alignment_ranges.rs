use crate::io::gene_map::GeneMap;
use crate::translate::coord_map::local_to_codon_range_exclusive;
use crate::translate::translate_genes::Translation;
use crate::utils::range::{intersect, NucRefGlobalRange, NucRefLocalRange, PositionLike};
use eyre::Report;

pub fn calculate_aa_alignment_ranges_in_place(
  global_alignment_range: &NucRefGlobalRange,
  translation: &mut Translation,
  gene_map: &GeneMap,
) -> Result<(), Report> {
  // For each translated CDS
  translation.iter_cdses_mut().try_for_each(|(_, cds_tr)| {
    // Get CDS annotation
    let cds = gene_map.get_cds(&cds_tr.name)?;

    let mut aa_alignment_ranges = vec![];
    let mut prev_segment_end = cds.segments[0].range.begin.as_usize();

    // For each segment
    for segment in &cds.segments {
      // Trim segment to include only what's inside alignment
      let included_range_global = intersect(global_alignment_range, &segment.range);
      if !included_range_global.is_empty() {
        // Convert to coordinates local to CDS (not local to segment!)
        let included_range_local = NucRefLocalRange::from_range(included_range_global - prev_segment_end as isize);
        aa_alignment_ranges.push(local_to_codon_range_exclusive(&included_range_local));
      }
      // CDS consists of concatenated segments; remember by how far we went along the CDS so far
      prev_segment_end += segment.len();
    }

    // Record computed AA alignment ranges on CDS translation
    // TODO: avoid mutable code. Move calculation of AA alignment ranges to where translation is.
    //   This requires global_alignment_range to be available there, but it is only computed much later currently.
    cds_tr.alignment_ranges = aa_alignment_ranges;

    Ok::<(), Report>(())
  })?;

  Ok(())
}
