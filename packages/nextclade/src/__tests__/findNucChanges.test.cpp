#include "../analyze/findNucChanges.h"

#include <common/safe_vector.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <string>

#include "../../include/nextclade/nextclade.h"
#include "../../include/nextclade/private/nextclade_private.h"

#define EXPECT_ARR_EQ(expected, actual) ASSERT_THAT(actual, ::testing::ElementsAreArray(expected));

using Nextclade::findNucChanges;
using Nextclade::NucleotideDeletion;
using Nextclade::NucleotideSubstitution;

TEST(FindNucChanges, ReportsAlignmentStartAndEnd) {
  std::stringstream input;

  // clang-format off
  const auto ref   = toNucleotideSequence("ACGTCAGTG");
  const auto query = toNucleotideSequence("--CTC-GT-");
  // clang-format on                       012345678

  const auto results = findNucChanges(ref, query);

  EXPECT_EQ(2, results.alignmentStart);
  EXPECT_EQ(8, results.alignmentEnd);
}

TEST(FindNucChanges, ReportsSubstitutions) {
  std::stringstream input;

  // clang-format off
  const auto ref   = toNucleotideSequence("CTA" "ATA" "GTA");
  const auto query = toNucleotideSequence("ATA" "TTA" "GTA");
  // clang-format on                       012   345   678

  const auto results = findNucChanges(ref, query);

  const auto expected = safe_vector<NucleotideSubstitution>({
    {.ref = Nucleotide::C, .pos = 0, .qry = Nucleotide::A},
    {.ref = Nucleotide::A, .pos = 3, .qry = Nucleotide::T},
  });

  EXPECT_ARR_EQ(expected, results.substitutions)
}

TEST(FindNucChanges, ReportsDeletions) {
  std::stringstream input;

  // clang-format off
  const auto ref   = toNucleotideSequence("CTA" "ATA" "GTA");
  const auto query = toNucleotideSequence("CTA" "---" "GTA");
  // clang-format on                       012   345   678

  const auto results = findNucChanges(ref, query);

  const auto expected = safe_vector<NucleotideDeletion>({
    {.start = 3, .length = 3}
  });

  EXPECT_ARR_EQ(expected, results.deletions)
}
