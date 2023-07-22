use crate::align::band_2d::full_matrix;
use crate::align::band_2d::Stripe;
use crate::align::params::AlignPairwiseParams;
use crate::align::seed_match::seed_match;
use crate::align::seed_match2::{get_seed_matches2, CodonSpacedIndex, SeedMatch2};
use crate::alphabet::letter::Letter;
use crate::alphabet::nuc::Nuc;
use crate::make_error;
use crate::utils::collections::first;
use eyre::Report;
use log::trace;
use num_traits::abs;
use num_traits::clamp;
use num_traits::clamp_min;
use std::cmp::max;
use std::cmp::min;
use std::collections::BTreeMap;

/// generate a vector of query sequence positions that are followed by at least `seed_length`
/// valid characters. Positions in this vector are thus "good" positions to start a query k-mer.
fn get_map_to_good_positions<L: Letter<L>>(qry_seq: &[L], seed_length: usize) -> Vec<usize> {
  let qry_len = qry_seq.len();

  let mut map_to_good_positions = Vec::<usize>::with_capacity(qry_len);
  let mut distance_to_last_bad_pos: i32 = 0;

  for (i, letter) in qry_seq.iter().enumerate() {
    // TODO: Exclude ambiguous letters
    if letter.is_unknown() {
      distance_to_last_bad_pos = 0;
    } else if distance_to_last_bad_pos >= seed_length as i32 {
      map_to_good_positions.push(i - seed_length);
    }
    distance_to_last_bad_pos += 1;
  }

  map_to_good_positions
}

#[derive(Debug, Clone, Copy)]
pub struct SeedMatch {
  pub qry_pos: usize,
  pub ref_pos: usize,
  pub score: usize,
}

/// Determine seed matches between query and reference sequence. will only attempt to
/// match k-mers without ambiguous characters. Search is performed via a left-to-right
/// search starting and the previous valid seed and extending at most to the maximally
/// allowed insertion/deletion (shift) distance.
pub fn get_seed_matches<L: Letter<L>>(
  qry_seq: &[L],
  ref_seq: &[L],
  params: &AlignPairwiseParams,
) -> (Vec<SeedMatch>, i32) {
  let mut seed_matches = Vec::<SeedMatch>::new();

  // get list of valid k-mer start positions
  let map_to_good_positions = get_map_to_good_positions(qry_seq, params.seed_length);
  let n_good_positions = map_to_good_positions.len();
  if n_good_positions < params.seed_length {
    return (seed_matches, 0);
  }

  // use 1/seed_spacing for long sequences, min_seeds otherwise
  let n_seeds = if ref_seq.len() > (params.min_seeds * params.seed_spacing) as usize {
    (ref_seq.len() as f32 / params.seed_spacing as f32) as i32
  } else {
    params.min_seeds
  };

  // distance of first seed from the end of the sequence (third of seed spacing)
  let margin = (ref_seq.len() as f32 / (n_seeds * 3) as f32).round() as i32;

  // Generate kmers equally spaced on the query
  let effective_margin = (margin as f32).min(n_good_positions as f32 / 4.0);
  let seed_cover = n_good_positions as f32 - 2.0 * effective_margin;
  let kmer_spacing = (seed_cover - 1.0) / ((n_seeds - 1) as f32);

  // loop over seeds and find matches, store in seed_matches
  let mut start_pos = 0; // start position of ref search
  let mut end_pos = ref_seq.len(); // end position of ref search
  let qry_pos = 0;

  for ni in 0..n_seeds {
    // pick index of of seed in map
    let good_position_index = (effective_margin + (kmer_spacing * ni as f32)).round() as usize;
    // get new query kmer-position
    let qry_pos_old = qry_pos;
    let qry_pos = map_to_good_positions[good_position_index];
    // increment upper bound for search in reference
    end_pos += qry_pos - qry_pos_old;

    //extract seed and find match
    let seed = &qry_seq[qry_pos..(qry_pos + params.seed_length)];
    let tmp_match = seed_match(seed, ref_seq, start_pos, end_pos, params.mismatches_allowed);

    // Only use seeds with at most allowed_mismatches
    if tmp_match.score >= params.seed_length - params.mismatches_allowed {
      // if this isn't the first match, check that reference position of current match after previous
      // if previous seed matched AFTER current seed, remove previous seed
      if let Some(prev_match) = seed_matches.last() {
        if tmp_match.ref_pos > prev_match.ref_pos {
          start_pos = prev_match.ref_pos;
        } else {
          // warn!("Crossed over seed removed. {:?}", prev_match);
          seed_matches.pop();
        }
      }

      let seed_match = SeedMatch {
        qry_pos,
        ref_pos: tmp_match.ref_pos,
        score: tmp_match.score,
      };

      // check that current seed matches AFTER previous seed (-2 if already triggered above)
      // and add current seed to list of matches
      if seed_matches
        .last()
        .map_or(true, |prev_match| prev_match.ref_pos < tmp_match.ref_pos)
      {
        //ensure seed positions increase strictly monotonically
        seed_matches.push(seed_match);
        // } else {
        //   warn!(
        //     "Seed not used because of identical ref_pos with previous seed: {:?}",
        //     seed_match
        //   );
      }
      // increment the "reference search end-pos" as the current reference + maximally allowed indel
      end_pos = tmp_match.ref_pos + params.max_indel;
    }
  }
  (seed_matches, n_seeds)
}

/// Determine rough positioning of qry to reference sequence by approximate seed matching
/// Returns vector of stripes, that is a band within which the alignment is expected to lie
pub fn seed_alignment(
  qry_seq: &[Nuc],
  ref_seq: &[Nuc],
  seed_index: &CodonSpacedIndex,
  params: &AlignPairwiseParams,
) -> Result<Vec<Stripe>, Report> {
  let qry_len = qry_seq.len();
  let ref_len = ref_seq.len();

  if ref_len + qry_len < (10 * params.seed_length) {
    // for very short sequences, use full square
    let stripes = full_matrix(ref_len, qry_len);
    trace!("Band construction: Short qry&ref sequence (< 5*seed_length), thus using full matrix");
    Ok(stripes)
  } else {
    // otherwise, determine seed matches roughly regularly spaced along the query sequence
    let seed_matches = get_seed_matches2(qry_seq, ref_seq, seed_index, params)?;
    create_stripes(
      &seed_matches,
      qry_len as isize,
      ref_len as isize,
      params.terminal_bandwidth as isize,
      params.excess_bandwidth as isize,
      params.max_indel as isize,
      params.allowed_mismatches as isize,
    )
  }
}

fn abs_shift(seed1: &SeedMatch2, seed2: &SeedMatch2) -> isize {
  abs(seed2.offset - seed1.offset)
}

/// Takes in seed matches and returns a vector of stripes
/// Stripes define the query sequence range for each reference position
pub fn create_stripes(
  chain: &[SeedMatch2],
  qry_len: isize,
  ref_len: isize,
  terminal_bandwidth: isize,
  excess_bandwidth: isize,
  max_indel: isize,
  allowed_mismatches: isize,
) -> Result<Vec<Stripe>, Report> {
  // This function steps through the chained seeds and determines and appropriate band
  // defined via stripes in query coordinates. These bands will later be chopped to reachable ranges

  // the broad idea is the following:
  // pre: deal with a special case the beginning and allow for terminal bandwidth
  // within: for each pair of chained seed (current and next), add a Trapezoid centered at the junction
  //         and a body Trapezoid for next. Extend the gap Trapezoid into the current and next
  // post: deal with the terminal trapezoid and allow of terminal bandwidth

  let mut broad_seeds = Vec::<TrapezoidDirectParams>::with_capacity(2 * chain.len() + 2);

  let mut current_seed = &chain[0];
  let mut current_seed_end = (current_seed.ref_pos + current_seed.length) as isize;
  // make initial trapezoid starting at 0 and extending into match by terminal_bandwidth
  let mut look_back_length = terminal_bandwidth;
  let mut look_forward_length = terminal_bandwidth;
  let mut current_broad_seed = TrapezoidDirectParams {
    ref_start: 0,
    ref_end: min(current_seed.ref_pos as isize + look_forward_length, ref_len),
    left_offset: current_seed.offset - terminal_bandwidth,
    right_offset: current_seed.offset + terminal_bandwidth,
  };

  // add body for first seed match starting were the previous ends
  if current_seed_end > current_broad_seed.ref_end {
    broad_seeds.push(current_broad_seed);
    current_broad_seed = TrapezoidDirectParams {
      ref_start: current_broad_seed.ref_end,
      ref_end: min(current_seed_end, ref_len + 1),
      left_offset: current_seed.offset - allowed_mismatches,
      right_offset: current_seed.offset + allowed_mismatches,
    };
  }

  // loop over remaining seeds in chain
  for chain_index in 1..chain.len() {
    let next_seed = &chain[chain_index];
    let mean_offset = (next_seed.offset + current_seed.offset) / 2; // offset of gap seed
    let shift = abs_shift(current_seed, next_seed) / 2; // distance from mean offset
    look_back_length = shift + excess_bandwidth;
    look_forward_length = shift + excess_bandwidth;
    // rewind the broad seeds until the ref_start of the last one preceeds the one too add
    while current_broad_seed.ref_start > max(0, current_seed_end - look_back_length) {
      current_broad_seed = if let Some(tmp_seed) = broad_seeds.pop() {
        tmp_seed
      } else {
        // we rewound all the way to the beginning, add a new terminal
        TrapezoidDirectParams {
          ref_start: 0,
          ref_end: min(current_seed.ref_pos as isize + look_forward_length, ref_len + 1),
          left_offset: current_seed.offset - look_back_length,
          right_offset: current_seed.offset + look_back_length,
        }
      };
      look_back_length = max(look_back_length, mean_offset - current_broad_seed.left_offset);
    }
    // terminate previous trapezoid where the new one will start and push
    if (current_seed_end - look_back_length) > 0 {
      current_broad_seed.ref_end = current_seed_end - look_back_length;
      broad_seeds.push(current_broad_seed);
    } else {
      current_broad_seed.ref_end = 0;
    }

    // generate trapezoid for the gap between seeds and push
    current_broad_seed = TrapezoidDirectParams {
      ref_start: current_broad_seed.ref_end,
      ref_end: next_seed.ref_pos as isize + look_forward_length,
      left_offset: mean_offset - look_back_length - excess_bandwidth,
      right_offset: mean_offset + look_back_length + excess_bandwidth,
    };
    broad_seeds.push(current_broad_seed);

    // generate new current trapezoid for the body of the next seed
    current_broad_seed = TrapezoidDirectParams {
      ref_start: current_broad_seed.ref_end,
      ref_end: (next_seed.ref_pos + next_seed.length) as isize,
      left_offset: next_seed.offset - allowed_mismatches,
      right_offset: next_seed.offset + allowed_mismatches,
    };
    current_seed = next_seed;
    current_seed_end = (current_seed.ref_pos + current_seed.length) as isize;
  }

  look_back_length = max(terminal_bandwidth, look_back_length);
  // rewind the broad seeds until the ref_start of the last one preceeds the one too add
  while current_broad_seed.ref_start > max(0, current_seed_end - look_back_length) {
    current_broad_seed = if let Some(tmp_seed) = broad_seeds.pop() {
      tmp_seed
    } else {
      // we rewound all the way to the beginning, add a new terminal
      TrapezoidDirectParams {
        ref_start: 0,
        ref_end: min(current_seed.ref_pos as isize + look_back_length, ref_len + 1),
        left_offset: current_seed.offset - look_back_length,
        right_offset: current_seed.offset + look_back_length,
      }
    };
    look_back_length = max(look_back_length, current_seed.offset - current_broad_seed.left_offset);
  }
  // terminate previous trapezoid where the new one will start and push
  current_broad_seed.ref_end = max(0, current_seed_end - look_back_length);
  broad_seeds.push(current_broad_seed);

  current_broad_seed = TrapezoidDirectParams {
    ref_start: current_broad_seed.ref_end,
    ref_end: ref_len + 1,
    left_offset: current_seed.offset - look_back_length,
    right_offset: current_seed.offset + look_back_length,
  };
  broad_seeds.push(current_broad_seed);

  let mut stripes = Vec::<Stripe>::with_capacity(ref_len as usize + 1);
  for broad_seed in broad_seeds {
    for ref_pos in broad_seed.ref_start..broad_seed.ref_end {
      stripes.push(Stripe {
        begin: min(qry_len - allowed_mismatches, max(0, ref_pos + broad_seed.left_offset)) as usize,
        end: min(qry_len + 1, ref_pos + broad_seed.right_offset) as usize,
      });
    }
  }
  // write_stripes_to_file(&stripes, "stripes.csv");

  let regularized_stripes = regularize_stripes(stripes, qry_len as usize);

  // For debugging of stripes and matches:
  // write_stripes_to_file(&regularized_stripes, "regularized_stripes.csv");
  // write_matches_to_file(seed_matches, "matches.csv");
  // Usefully visualized using `python scripts/visualize-stripes.py`
  //
  trace_stripe_stats(&regularized_stripes);

  Ok(regularized_stripes)
}

#[derive(Clone, Copy, Debug)]
// Default qry_len is qry_len
struct TrapezoidOffsetParams {
  ref_start: isize,
  ref_end: isize,
  offset: isize,
  bandwidth: usize,
}

#[derive(Clone, Copy, Debug)]
struct TrapezoidDirectParams {
  ref_start: isize,
  ref_end: isize,
  left_offset: isize,
  right_offset: isize,
}

// Implement a function on TrapezoidOffsetParams to convert to TrapezoidDirectParams
impl TrapezoidOffsetParams {
  const fn to_direct_params(self) -> TrapezoidDirectParams {
    TrapezoidDirectParams {
      ref_start: self.ref_start,
      ref_end: self.ref_end,
      left_offset: self.offset - self.bandwidth as isize,
      right_offset: self.offset + self.bandwidth as isize,
    }
  }
}

/// Chop off unreachable parts of the stripes.
/// Overhanging parts are pruned
fn regularize_stripes(mut stripes: Vec<Stripe>, qry_len: usize) -> Vec<Stripe> {
  // assure stripe begin are non-decreasing -- such states would be unreachable in the alignment
  let stripes_len = stripes.len();
  stripes[0].begin = 0;
  for i in 1..stripes_len {
    stripes[i].begin = clamp(stripes[i].begin, stripes[i - 1].begin, qry_len);
  }

  // analogously, assure that strip ends are non-decreasing. this needs to be done in reverse.
  stripes[stripes_len - 1].end = qry_len + 1;
  for i in (0..(stripes_len - 1)).rev() {
    stripes[i].end = clamp(stripes[i].end, 1, stripes[i + 1].end);
  }

  stripes
}

fn trace_stripe_stats(stripes: &[Stripe]) {
  let mut stripe_lengths = Vec::new();
  for stripe in stripes {
    assert!(
      stripe.begin <= stripe.end,
      "Stripe begin must be <= stripe end for stripe {stripe:?}",
    );
    stripe_lengths.push(stripe.end - stripe.begin);
  }
  stripe_lengths.sort_unstable();
  let median = stripe_lengths[stripe_lengths.len() / 2];
  let mean = stripe_lengths.iter().sum::<usize>() as f32 / stripe_lengths.len() as f32;
  let max = stripe_lengths[stripe_lengths.len() - 1];
  let min = stripe_lengths[0];
  trace!("Stripe width stats: min: {min}, max: {max}, mean: {mean:.1}, median: {median}",);
}

fn trace_matches(matches: &[SeedMatch2]) {
  for (i, seed) in matches.iter().enumerate() {
    trace!(
      "Match {}: ref_pos: {}, qry_offset: {}, length: {}",
      i,
      seed.ref_pos,
      -seed.offset,
      seed.length,
    );
  }
}

fn write_stripes_to_file(stripes: &[Stripe], filename: &str) {
  use std::io::Write;
  let mut file = std::fs::File::create(filename).unwrap();
  writeln!(file, "ref,begin,end").unwrap();
  for (i, stripe) in stripes.iter().enumerate() {
    writeln!(file, "{i},{begin},{end}", begin = stripe.begin, end = stripe.end).unwrap();
  }
}

pub fn write_matches_to_file(matches: &[SeedMatch2], filename: &str) {
  use std::io::Write;
  let mut file = std::fs::File::create(filename).unwrap();
  writeln!(file, "ref_pos,qry_pos,length").unwrap();
  for match_ in matches {
    writeln!(file, "{},{},{}", match_.ref_pos, match_.qry_pos, match_.length).unwrap();
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Report;
  use rstest::rstest;

  #[rstest]
  fn test_create_stripes_basic() -> Result<(), Report> {
    let seed_matches = vec![
      SeedMatch2 {
        qry_pos: 5,
        ref_pos: 10,
        length: 0,
        offset: 0,
      },
      SeedMatch2 {
        qry_pos: 20,
        ref_pos: 30,
        length: 0,
        offset: 0,
      },
    ];

    let terminal_bandwidth = 5;
    let excess_bandwidth = 2;
    let allowed_mismatches = 2;
    let max_indel = 100;
    let qry_len = 30;
    let ref_len = 40;

    let result = create_stripes(
      &seed_matches,
      qry_len,
      ref_len,
      terminal_bandwidth,
      excess_bandwidth,
      max_indel,
      allowed_mismatches,
    );

    Ok(())
  }
}
