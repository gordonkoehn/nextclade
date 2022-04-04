use crate::io::aa::Aa;
use crate::io::letter::Letter;
use crate::io::nuc::Nuc;
use crate::translate::translate_genes::Translation;
use crate::utils::error::keep_ok;
use crate::utils::range::Range;
use eyre::Report;
use map_in_place::MapVecInPlace;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LetterRange<L: Letter<L>> {
  pub begin: usize,
  pub end: usize,
  pub letter: L,
}

impl<L: Letter<L>> LetterRange<L> {
  pub fn contains_pos(&self, x: usize) -> bool {
    x >= self.begin && x < self.end
  }
}

impl<L: Letter<L>> LetterRange<L> {
  #[inline]
  pub fn len(&self) -> usize {
    self.end - self.begin
  }

  #[inline]
  pub fn is_empty(&self) -> bool {
    self.len() == 0
  }

  #[inline]
  pub fn to_range(&self) -> Range {
    Range {
      begin: self.begin,
      end: self.end,
    }
  }
}

pub type NucRange = LetterRange<Nuc>;
pub type AaRange = LetterRange<Aa>;

// Finds contiguous ranges (segments) in the sequence, such that for every character inside every range,
// the predicate function returns true and every range contains only the same letter.
//
// The predicate is a function that takes a character and returns boolean.
//
// For example if predicate returns `true` for characters A and C, this function will find ranges `AAAA` and `CCCCC`,
// but not `ZZZ` or `ACCCAC`.
pub fn find_letter_ranges_by<L: Letter<L>>(seq: &[L], pred: impl Fn(L) -> bool) -> Vec<LetterRange<L>> {
  let len = seq.len();

  let mut result = Vec::<LetterRange<L>>::new();
  let mut i = 0_usize;
  let mut begin = 0_usize;
  let mut found_maybe = Option::<L>::default();
  while i < len {
    let letter = seq[i];

    // Find beginning of a range
    if pred(letter) {
      begin = i;
      found_maybe = Some(letter);
    }

    match found_maybe {
      // If there's a current range we are working on (for which we found a `begin`), extend it
      Some(found) => {
        // Rewind forward until we find a mismatch
        while i < len && seq[i] == found {
          i += 1;
        }

        // We found the end of the current range, so now it's complete
        let end = i;

        // Remember the range
        result.push(LetterRange::<L> { begin, end, letter });

        found_maybe = None;
      }
      None => {
        if i < len {
          i += 1;
        }
      }
    }
  }
  result
}

/// Finds contiguous ranges (segments) consisting of a given nucleotide in the sequence.
pub fn find_letter_ranges<L: Letter<L>>(qry_aln: &[L], letter: L) -> Vec<LetterRange<L>> {
  find_letter_ranges_by(qry_aln, |candidate| candidate == letter)
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct GeneAaRange {
  pub gene_name: String,
  pub letter: Aa,
  pub ranges: Vec<AaRange>,
  pub length: usize,
}

impl GeneAaRange {
  pub fn contains_pos(&self, pos: usize) -> bool {
    self.ranges.iter().any(|range| range.contains_pos(pos))
  }
}

/// Finds contiguous ranges (segments) consisting of a given amino acid in the sequence.
pub fn find_aa_letter_ranges(translations: &[Result<Translation, Report>], letter: Aa) -> Vec<GeneAaRange> {
  keep_ok(translations)
    .map(|Translation { gene_name, seq, .. }| {
      let ranges = find_letter_ranges_by(seq, |candidate| candidate == letter);
      let length = ranges.iter().map(LetterRange::len).sum();
      GeneAaRange {
        gene_name: gene_name.clone(),
        letter,
        ranges,
        length,
      }
    })
    .collect()
}
