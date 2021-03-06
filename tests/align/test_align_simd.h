// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================

#ifndef TESTS_ALIGN_TEST_ALIGN_SIMD_H_
#define TESTS_ALIGN_TEST_ALIGN_SIMD_H_

#include <tuple>

#include <seqan/basic.h>
#include <seqan/align.h>

#include "test_mock.h"

namespace impl
{
namespace test_align_simd
{

struct LocalAlignTester_
{
    template <typename TAlign,
              typename TScoreValue, typename TScoreSpec,
              typename TConfig>
    static auto
    run(TAlign & align,
        seqan::Score<TScoreValue, TScoreSpec> const & score,
        TConfig const &,
        int const lDiag,
        int const uDiag)
    {
        if (lDiag == std::numeric_limits<int>::min() && uDiag == std::numeric_limits<int>::max())
            return localAlignment(align, score);
        else
            return localAlignment(align, score, lDiag, uDiag);
    }
};

struct GlobalAlignTester_
{
    template <typename TAlign,
              typename TScoreValue, typename TScoreSpec,
              typename TConfig>
    static auto
    run(TAlign & align,
        seqan::Score<TScoreValue, TScoreSpec> const & score,
        TConfig const & config,
        int const lDiag,
        int const uDiag)
    {
        if (lDiag == std::numeric_limits<int>::min() && uDiag == std::numeric_limits<int>::max())
            return globalAlignment(align, score, config);
        else
            return globalAlignment(align, score, config, lDiag, uDiag);
    }
};

struct GlobalAlignScoreTester_
{
    template <typename TStringsH,
              typename TStringsV,
              typename TScoreValue, typename TScoreSpec,
              typename TConfig>
    static auto
    run(TStringsH const & strH,
        TStringsV const & strV,
        seqan::Score<TScoreValue, TScoreSpec> const & score,
        TConfig const & config,
        int const lDiag,
        int const uDiag)
    {
        if (lDiag == std::numeric_limits<int>::min() && uDiag == std::numeric_limits<int>::max())
            return globalAlignmentScore(strH, strV, score, config);
        else
            return globalAlignmentScore(strH, strV, score, config, lDiag, uDiag);
    }
};

struct LocalScoreTester_
{
    template <typename TStringsH,
              typename TStringsV,
              typename TScoreValue, typename TScoreSpec,
              typename TConfig>
    static auto
    run(TStringsH const & strH,
        TStringsV const & strV,
        seqan::Score<TScoreValue, TScoreSpec> const & score,
        TConfig const & /*config*/,
        int const lDiag,
        int const uDiag)
    {
        if (lDiag == seqan::MinValue<int>::VALUE && uDiag == seqan::MaxValue<int>::VALUE)
            return localAlignmentScore(strH, strV, score);
        else
            return localAlignmentScore(strH, strV, score, lDiag, uDiag);
    }
};
}  // namespace test_align_simd
}  // namespace impl

// ----------------------------------------------------------------------------
// Class SimdAlignTest
// ----------------------------------------------------------------------------

// Common test class instance, which stores the types to be accessed.
template <typename TTuple>
class SimdAlignTest : public seqan::Test
{
public:
    using TAlignConfig = std::tuple_element_t<0, TTuple>;
    using TLengthParam = std::tuple_element_t<1, TTuple>;
    using TBandSwitch = std::tuple_element_t<2, TTuple>;
};

// ----------------------------------------------------------------------------
// Configuration of typed tests for global alignment.
// ----------------------------------------------------------------------------

template <typename T>
class SimdAlignTestCommon : public SimdAlignTest<T>
{};

typedef
        seqan::TagList<std::tuple<seqan::AlignConfig<>,                         impl::test_align_mock::EqualLengthSimd,     seqan::BandOff>,
        seqan::TagList<std::tuple<seqan::AlignConfig<>,                         impl::test_align_mock::VariableLengthSimd,  seqan::BandOff>,
        seqan::TagList<std::tuple<seqan::AlignConfig<true, true, true, true>,   impl::test_align_mock::EqualLengthSimd,     seqan::BandOff>,
        seqan::TagList<std::tuple<seqan::AlignConfig<true, true, true, true>,   impl::test_align_mock::VariableLengthSimd,  seqan::BandOff>,
        seqan::TagList<std::tuple<seqan::AlignConfig<true, false, false, true>, impl::test_align_mock::EqualLengthSimd,     seqan::BandOff>,
        seqan::TagList<std::tuple<seqan::AlignConfig<true, false, false, true>, impl::test_align_mock::VariableLengthSimd,  seqan::BandOff>,
        seqan::TagList<std::tuple<seqan::AlignConfig<false, true, true, false>, impl::test_align_mock::EqualLengthSimd,     seqan::BandOff>,
        seqan::TagList<std::tuple<seqan::AlignConfig<false, true, true, false>, impl::test_align_mock::VariableLengthSimd,  seqan::BandOff>,
        seqan::TagList<std::tuple<seqan::AlignConfig<>,                         impl::test_align_mock::EqualLengthSimd,     seqan::BandOn>,
        seqan::TagList<std::tuple<seqan::AlignConfig<true, true, true, true>,   impl::test_align_mock::EqualLengthSimd,     seqan::BandOn>,
        seqan::TagList<std::tuple<seqan::AlignConfig<true, false, false, true>, impl::test_align_mock::EqualLengthSimd,     seqan::BandOn>,
        seqan::TagList<std::tuple<seqan::AlignConfig<false, true, true, false>, impl::test_align_mock::EqualLengthSimd,     seqan::BandOn>
        > > > > > > > > > > > > SimdAlignTestCommonCommonTypes;

SEQAN_TYPED_TEST_CASE(SimdAlignTestCommon, SimdAlignTestCommonCommonTypes);

// ----------------------------------------------------------------------------
// Configuration of typed tests for local alignment.
// ----------------------------------------------------------------------------

template <typename T>
class SimdAlignLocalTestCommon : public SimdAlignTest<T>
{};

typedef
        seqan::TagList<std::tuple<seqan::AlignConfig<>, impl::test_align_mock::EqualLengthSimd,    seqan::BandOff>,
        seqan::TagList<std::tuple<seqan::AlignConfig<>, impl::test_align_mock::VariableLengthSimd, seqan::BandOff>,
        seqan::TagList<std::tuple<seqan::AlignConfig<>, impl::test_align_mock::EqualLengthSimd,    seqan::BandOn>
        > > > SimdAlignLocalTestCommonCommonTypes;

SEQAN_TYPED_TEST_CASE(SimdAlignLocalTestCommon, SimdAlignLocalTestCommonCommonTypes);

// ----------------------------------------------------------------------------
// Function testAlignSimd()
// ----------------------------------------------------------------------------

template <typename TAlphabet,
          typename TFunctor,
          typename TScoreValue, typename TScoreSpec,
          typename TAlignConfig,
          typename TSimdLength>
void testAlignSimd(TFunctor const &,
                   seqan::Score<TScoreValue, TScoreSpec> const & score,
                   TAlignConfig const & config,
                   TSimdLength const & /*tag*/,
                   int const lDiag = std::numeric_limits<int>::min(),
                   int const uDiag = std::numeric_limits<int>::max())
{
    auto sets = impl::test_align_mock::TestSequences_<TAlphabet, TSimdLength>::getSequences();

    // Prepare an align object with the sequences.
    seqan::StringSet<seqan::Align<seqan::String<TAlphabet> > > alignments;
    resize(alignments, length(std::get<0>(sets)));
    auto zipCont = makeZipView(alignments, std::get<0>(sets), std::get<1>(sets));

    for(auto tuple : zipCont)
    {
        resize(rows(std::get<0>(tuple)), 2);
        assignSource(row(std::get<0>(tuple), 0), std::get<1>(tuple));
        assignSource(row(std::get<0>(tuple), 1), std::get<2>(tuple));
    }

    // Run the SIMD accelerated alignment.
    seqan::String<TScoreValue> scores = TFunctor::run(alignments, score, config, lDiag, uDiag);
    SEQAN_ASSERT_EQ(length(scores), length(alignments));

    // Check correctness of alignments using sequential alignment.
    auto zipRes = makeZipView(scores, alignments);
    for(auto res : zipRes)
    {
        typename std::decay<decltype(std::get<1>(res))>::type goldAlign;
        resize(rows(goldAlign), 2);
        assignSource(row(goldAlign, 0), source(row(std::get<1>(res), 0)));
        assignSource(row(goldAlign, 1), source(row(std::get<1>(res), 1)));

        TScoreValue goldScore = TFunctor::run(goldAlign, score, config, lDiag, uDiag);

        SEQAN_ASSERT_EQ(std::get<0>(res), goldScore);
        SEQAN_ASSERT(row(std::get<1>(res), 0) == row(goldAlign, 0));
        SEQAN_ASSERT(row(std::get<1>(res), 1) == row(goldAlign, 1));
    }
}

// Helper function to set band parameters.
template <typename TAlphabet,
          typename TFunctor,
          typename TScoreValue, typename TScoreSpec,
          typename TAlignConfig,
          typename TSimdLength,
          typename TBandFlag>
void testAlignSimd(TFunctor const &,
                   seqan::Score<TScoreValue, TScoreSpec> const & score,
                   TAlignConfig const & config,
                   TSimdLength const & /*tag*/,
                   TBandFlag const &)
{
    if (seqan::IsSameType<TBandFlag, seqan::BandOff>::VALUE)
        testAlignSimd<TAlphabet>(TFunctor(), score, config, TSimdLength());
    else
        testAlignSimd<TAlphabet>(TFunctor(), score, config, TSimdLength(), -4, 6);
}

// ----------------------------------------------------------------------------
// Function testAlignScoreSimd()
// ----------------------------------------------------------------------------

template <typename TAlphabet,
          typename TTester,
          typename TScoreValue, typename TScoreSpec,
          typename TAlignConfig,
          typename TSimdLength>
void testAlignSimdScore(TTester const &,
                        seqan::Score<TScoreValue, TScoreSpec> const & score,
                        TAlignConfig const & config,
                        TSimdLength const & /*tag*/,
                        int const lDiag = std::numeric_limits<int>::min(),
                        int const uDiag = std::numeric_limits<int>::max())
{
    auto sets = impl::test_align_mock::TestSequences_<TAlphabet, TSimdLength>::getSequences();

    seqan::String<TScoreValue> scores = TTester::run(std::get<0>(sets), std::get<1>(sets), score, config, lDiag, uDiag);

    SEQAN_ASSERT_EQ(length(scores), length(std::get<0>(sets)));

    auto zipRes = makeZipView(scores, std::get<0>(sets), std::get<1>(sets));
    for(auto res : zipRes)
    {
        TScoreValue goldScore = TTester::run(std::get<1>(res), std::get<2>(res), score, config, lDiag, uDiag);
        SEQAN_ASSERT_EQ(std::get<0>(res), goldScore);
    }
}

// Helper function to set band parameters.
template <typename TAlphabet,
          typename TFunctor,
          typename TScoreValue, typename TScoreSpec,
          typename TAlignConfig,
          typename TSimdLength,
          typename TBandFlag>
void testAlignSimdScore(TFunctor const &,
                        seqan::Score<TScoreValue, TScoreSpec> const & score,
                        TAlignConfig const & config,
                        TSimdLength const & /*tag*/,
                        TBandFlag const &)
{
    if (seqan::IsSameType<TBandFlag, seqan::BandOff>::VALUE)
        testAlignSimdScore<TAlphabet>(TFunctor(), score, config, TSimdLength());
    else
        testAlignSimdScore<TAlphabet>(TFunctor(), score, config, TSimdLength(), -4, 6);
}

// ----------------------------------------------------------------------------
// Global Alignments.
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(SimdAlignTestCommon, Linear_Align)
{
    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    testAlignSimd<seqan::Dna>(impl::test_align_simd::GlobalAlignTester_(), seqan::Score<int>(2, -1, -1),
                              TAlignConf(), TLengthParam(), TBandSwitch());
    testAlignSimd<seqan::AminoAcid>(impl::test_align_simd::GlobalAlignTester_(), seqan::Blosum62(-2),
                                    TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignTestCommon, Linear_Score)
{
    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    testAlignSimdScore<seqan::Dna>(impl::test_align_simd::GlobalAlignScoreTester_(), seqan::Score<int>(2, -1, -1),
                                   TAlignConf(), TLengthParam(), TBandSwitch());
    testAlignSimdScore<seqan::AminoAcid>(impl::test_align_simd::GlobalAlignScoreTester_(), seqan::Blosum62(-2),
                                         TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignTestCommon, Affine_Align)
{
    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    testAlignSimd<seqan::Dna>(impl::test_align_simd::GlobalAlignTester_(), seqan::Score<int>(2, -1, -1, -3),
                              TAlignConf(), TLengthParam(), TBandSwitch());
    testAlignSimd<seqan::AminoAcid>(impl::test_align_simd::GlobalAlignTester_(), seqan::Blosum62(-2, -4),
                                    TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignTestCommon, Dynamic_Score_Matrix_Align)
{
    using namespace seqan;

    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable> > dynScore{-2, -4};
    setScoreMatrixById(dynScore, AminoAcidScoreMatrixID::PAM120);

    testAlignSimd<seqan::AminoAcid>(::impl::test_align_simd::GlobalAlignTester_(), dynScore,
                                    TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignTestCommon, Affine_Score)
{
    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    testAlignSimdScore<seqan::Dna>(impl::test_align_simd::GlobalAlignScoreTester_(), seqan::Score<int>(2, -1, -1, -3),
                                   TAlignConf(), TLengthParam(), TBandSwitch());
    testAlignSimdScore<seqan::AminoAcid>(impl::test_align_simd::GlobalAlignScoreTester_(), seqan::Blosum62(-2, -4),
                                         TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignTestCommon, Dynamic_Score_Matrix)
{
    using namespace seqan;

    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable> > dynScore{-2, -4};
    setScoreMatrixById(dynScore, AminoAcidScoreMatrixID::PAM120);

    testAlignSimdScore<seqan::AminoAcid>(::impl::test_align_simd::GlobalAlignScoreTester_(), dynScore,
                                         TAlignConf(), TLengthParam(), TBandSwitch());
}

// ----------------------------------------------------------------------------
// Local Alignments.
// ----------------------------------------------------------------------------

SEQAN_TYPED_TEST(SimdAlignLocalTestCommon, Linear_Align)
{
    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    testAlignSimd<seqan::Dna>(impl::test_align_simd::LocalAlignTester_(), seqan::Score<int>(2, -1, -1),
                              TAlignConf(), TLengthParam(), TBandSwitch());
    testAlignSimd<seqan::AminoAcid>(impl::test_align_simd::LocalAlignTester_(), seqan::Blosum62(-2),
                                    TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignLocalTestCommon, Linear_Score)
{
    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    testAlignSimdScore<seqan::Dna>(impl::test_align_simd::LocalScoreTester_(), seqan::Score<int>(2, -1, -1),
                                   TAlignConf(), TLengthParam(), TBandSwitch());
    testAlignSimdScore<seqan::AminoAcid>(impl::test_align_simd::LocalScoreTester_(), seqan::Blosum62(-2),
                                         TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignLocalTestCommon, Affine_Align)
{
    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    testAlignSimd<seqan::Dna>(impl::test_align_simd::LocalAlignTester_(), seqan::Score<int>(2, -1, -1, -3),
                              TAlignConf(), TLengthParam(), TBandSwitch());
    testAlignSimd<seqan::AminoAcid>(impl::test_align_simd::LocalAlignTester_(), seqan::Blosum62(-2, -4),
                                    TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignLocalTestCommon, Dynamic_Score_Matrix_Align)
{
    using namespace seqan;

    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable> > dynScore{-2, -4};
    setScoreMatrixById(dynScore, AminoAcidScoreMatrixID::PAM120);

    testAlignSimd<seqan::AminoAcid>(::impl::test_align_simd::LocalAlignTester_(), dynScore,
                                    TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignLocalTestCommon, Affine_Score)
{
    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    testAlignSimdScore<seqan::Dna>(impl::test_align_simd::LocalScoreTester_(), seqan::Score<int>(2, -1, -1, -3),
                                   TAlignConf(), TLengthParam(), TBandSwitch());
    testAlignSimdScore<seqan::AminoAcid>(impl::test_align_simd::LocalScoreTester_(), seqan::Blosum62(-2, -4),
                                         TAlignConf(), TLengthParam(), TBandSwitch());
}

SEQAN_TYPED_TEST(SimdAlignLocalTestCommon, Dynamic_Score_Matrix)
{
    using namespace seqan;

    using TAlignConf = typename TestFixture::TAlignConfig;
    using TLengthParam = typename TestFixture::TLengthParam;
    using TBandSwitch = typename TestFixture::TBandSwitch;

    Score<int, ScoreMatrix<AminoAcid, ScoreSpecSelectable> > dynScore{-2, -4};
    setScoreMatrixById(dynScore, AminoAcidScoreMatrixID::PAM120);

    testAlignSimdScore<seqan::AminoAcid>(::impl::test_align_simd::LocalScoreTester_(), dynScore,
                                         TAlignConf(), TLengthParam(), TBandSwitch());
}

#endif  // #ifndef TESTS_ALIGN_TEST_ALIGN_SIMD_H_
