//  smithWaterman.cpp
//  Skippy _search
//
// VERSION CONTROL:
// Edit scoring_config.h to switch between scoring methods

#include "smithWaterman.hpp"

#include "encode_kmer_2bit.h"
#include "scoring_config.h" // Single source of truth for version control
#include "ssw_cpp.h"
// #include "stripe_pattern.h"
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <sstream>
#include <string>
#include <string_view>

#include "stripe_pattern.hpp"

#define USE_ABSOLUTE_MATCH_COUNT

// Forward declaration for edlib NW alignment (implemented in
// alignSequences_edlib.cpp)
float alignNW_identity(const std::string_view& refWindow, const std::string& query);

// NW alignment via edlib
float alignNW_identity(const std::string_view& refWindow, const std::string& query);
// Forward declaration — implemented in alignSequences_edlib.cpp
float alignNW_identity(const std::string_view& refWindow, const std::string& query);

// Globals defined in main.cpp
extern bool MATCH_ALL;
extern bool GET_SW;

// You already have this global elsewhere; we only *use* it here.
extern uint32_t GATE_COVERAGE;
extern StripePattern gPattern;

// 'r' is an alias for 'Result' (defined in SafeQueue.hpp)
using r = Result;

// Configurable thresholds for density-based scoring
#ifdef USE_ALIGNMENT_DENSITY
extern uint32_t MIN_ALIGNMENT_SPAN;
extern float IDENTITY_THRESHOLD;
extern bool BALANCED_MODE;
extern bool BALANCED_MODE;
extern uint32_t SUSPICIOUS_SPAN_THRESHOLD;
extern float MIN_ALIGNMENT_DENSITY;
extern uint32_t MIN_ABSOLUTE_MATCHES;
#endif

// -----------------------------------------------------------------------------
// Utilities
// -----------------------------------------------------------------------------

std::string toUpperCase(const std::string& str)
{
    std::string upper = str;
    std::transform(upper.begin(), upper.end(), upper.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return upper;
}

std::string reverseComplement(const std::string& sequence)
{
    std::string rev;
    rev.reserve(sequence.size());

    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        switch (*it) {
        case 'A':
        case 'a':
            rev += 'T';
            break;
        case 'T':
        case 't':
            rev += 'A';
            break;
        case 'C':
        case 'c':
            rev += 'G';
            break;
        case 'G':
        case 'g':
            rev += 'C';
            break;
        default:
            rev += *it;
            break;
        }
    }
    return rev;
}

static inline void rstrip_inplace(std::string& s)
{
    while (!s.empty()) {
        const char c = s.back();
        if (c == '\n' || c == '\r' || c == ' ' || c == '\t')
            s.pop_back();
        else
            break;
    }
}

// Count exact matches as '=' in CIGAR. 'M' is ignored (ambiguous).
static inline uint32_t countEqualsInCigar(const std::string& cigar)
{
    uint32_t eq = 0;
    uint32_t num = 0;

    for (size_t i = 0; i < cigar.size(); ++i) {
        const char c = cigar[i];
        if (c >= '0' && c <= '9') {
            num = num * 10u + static_cast<uint32_t>(c - '0');
            continue;
        }

        if (c == '=')
            eq += num;

        num = 0;
    }

    return eq;
}

uint32_t aligned_len_from_cigar(const std::string& cigar)
{
    uint32_t aligned = 0;
    uint32_t num = 0;

    for (char c : cigar) {
        if (c >= '0' && c <= '9') {
            num = num * 10 + (c - '0');
        }
        else {
            switch (c) {
            case 'M': // match or mismatch
            case '=': // match
            case 'X': // mismatch
            case 'I': // insertion consumes query
                aligned += num;
                break;

            case 'D': // deletion consumes reference only
            case 'S': // soft clip
            case 'H': // hard clip
            default:
                break;
            }
            num = 0;
        }
    }

    return aligned;
}

#ifdef USE_ALIGNMENT_DENSITY
// -----------------------------------------------------------------------------
// NEW: Compute alignment span from CIGAR
// Span = total reference positions covered (matches + mismatches + deletions)
// -----------------------------------------------------------------------------
static inline uint32_t computeAlignmentSpan(const std::string& cigar)
{
    uint32_t span = 0;
    uint32_t num = 0;

    for (char c : cigar) {
        if (c >= '0' && c <= '9') {
            num = num * 10 + (c - '0');
        }
        else {
            // Operations that consume reference positions
            switch (c) {
            case 'M': // match or mismatch
            case '=': // match
            case 'X': // mismatch
            case 'D': // deletion (consumes reference)
                span += num;
                break;

            case 'I': // insertion (consumes query only)
            case 'S': // soft clip
            case 'H': // hard clip
            default:
                break;
            }
            num = 0;
        }
    }

    return span;
}
#endif

// -----------------------------------------------------------------------------
// reverse_complement - for k-mer operations
// -----------------------------------------------------------------------------
uint64_t reverse_complement(uint64_t kmer, uint64_t bases)
{
    uint64_t rev = 0;
    for (uint64_t i = 0; i < bases; ++i) {
        const uint64_t base = (kmer >> (2 * i)) & 3ULL;
        rev |= (base ^ 3ULL) << (2 * (bases - 1 - i));
    }
    return rev;
}

// -----------------------------------------------------------------------------
// alignSequences
// -----------------------------------------------------------------------------
void alignSequences(const std::string_view& refWindow, const std::string& query, std::string& cigar,
                    int& sw_score, uint32_t& ref_begin_local, uint32_t& ref_end_local)
{
    StripedSmithWaterman::Aligner aligner(2, 2, 3, 1);
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;

    std::string refStr(refWindow);

    // Correct API: Align(query, ref, ref_len, filter, alignment, maskLen)
    aligner.Align(query.c_str(), refStr.c_str(), static_cast<int>(refStr.size()), filter,
                  &alignment, 15);

    cigar = alignment.cigar_string;
    sw_score = alignment.sw_score;
    ref_begin_local = alignment.ref_begin;
    ref_end_local = alignment.ref_end;
}

// -----------------------------------------------------------------------------
// getSW - Main alignment function
//
// This function has TWO implementations controlled by macro:
//   1. USE_ABSOLUTE_MATCH_COUNT: Original scoring (count of '=' in CIGAR)
//   2. USE_ALIGNMENT_DENSITY: New density-based scoring
//
// candidate_results tuple layout (REQUIRED):
//   (localIdx, refBegin, coverage, span, revCompFlag)
//
// Mode A: GET_SW=true, cutoff>0  => SW always, gate on score >= cutoff, output
// SW. Mode B: GET_SW=true, cutoff==0 => SW always, no gate, output SW. Mode C:
// GET_SW=false           => SW only in twilight zone (low cov + low span),
//                                  gate on score >= cutoff (if cutoff>0),
//                                  output coverage-based (cigar="*").
// -----------------------------------------------------------------------------
bool getSW(const char* blob, const std::vector<ReadInfo>& readInfos, size_t startIndex,
           const std::span<char>& genomeStr, SafeQueue<r>& queue,
           const std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, uint32_t>>&
               candidate_results,
           int matchScore, int mismatchScore, int gapScore, int cutoff)
{
    (void)matchScore;
    (void)mismatchScore;
    (void)gapScore;

    if (candidate_results.empty())
        return false;

    // ------------------------------------------------------------
    // Twilight-zone policy (read-length independent)
    //
    // FP signature (your model):
    //   - marginal coverage
    //   - localised (small span)
    //
    // If coverage is marginal BUT span is wide, we do NOT try to "fix" it
    // because that can be a noisy human read signature.
    // ------------------------------------------------------------

    // "Marginal coverage": allow a little headroom above GATE_COVERAGE.
    // One extra seed-worth is a reasonable first cut; tune later.
    const uint32_t TWILIGHT_COV_MAX = GATE_COVERAGE + (uint32_t)gPattern.weight;

    // "Localised": span near the minimum implied by coverage.
    // Since coverage counts unique striped bases, span should not be much bigger
    // than coverage. Allow some slack for INDEL_SLACK / jitter / seed selection.
    const uint32_t TWILIGHT_SPAN_MAX = GATE_COVERAGE + 10;

    // IMPORTANT: Even when cutoff==0 we still want SW cleanup in twilight.
    // So twilight veto uses a deterministic bp threshold not tied to cutoff.
    // Start with GATE_COVERAGE (same scale as your coverage gate; read-length
    // independent).
    const uint32_t TWILIGHT_EQ_MIN = GATE_COVERAGE;

    const bool want_sw_global = GET_SW;

#ifdef USE_ABSOLUTE_MATCH_COUNT
    // ORIGINAL: cutoff is absolute match count threshold
    const bool gate_sw_global = (cutoff > 0);
#else
    // DENSITY: cutoff disabled in favor of density thresholds
    const bool gate_sw_global = false; // Use density gate instead
#endif

    for (const auto& cand : candidate_results) {
        const uint32_t localIdx = std::get<0>(cand); // chunk-local read index
        const uint32_t refBegin = std::get<1>(cand);
        const uint32_t coverage = std::get<2>(cand);
        const uint32_t span = std::get<3>(cand);
        const uint32_t revCompFlag = std::get<4>(cand);

        if (localIdx >= readInfos.size())
            continue;

        const uint32_t globalSampleIndex = static_cast<uint32_t>(startIndex + localIdx);
        const ReadInfo& ri = readInfos[localIdx];

        std::string readSeq(blob + ri.start, blob + ri.start + ri.length);
        rstrip_inplace(readSeq);
        if (readSeq.empty())
            continue;
        uint32_t read_len = static_cast<uint32_t>(readSeq.size());

        // ------------------------------------------------------------
        // MODE C: GET_SW == false
        //   - SW only in twilight zone (low coverage AND low span)
        //   - SW is a veto gate (cleanup), output remains coverage-based.
        //   - SW must run even if cutoff==0 (use TWILIGHT_EQ_MIN).
        // ------------------------------------------------------------
        if (!want_sw_global) {
            const bool twilight = (coverage <= TWILIGHT_COV_MAX) && (span <= TWILIGHT_SPAN_MAX);

            //            if (!twilight)
            //	    if (span > 79)
            // Balanced mode: length-dependent secondary gate in getSW
            if (BALANCED_MODE) {
                const uint32_t balancedGate =
                    std::max(33u, (uint32_t)(7.2f * std::log((float)read_len) + -3.4413f + 0.5f));
                if (coverage < balancedGate)
                    continue;
            }

            // Fast-accept only if coverage is good AND span is not suspiciously
            // short. Reads with span below threshold fall through to SW regardless of
            // cutoff.
            const uint32_t SUSPICIOUS_SPAN = SUSPICIOUS_SPAN_THRESHOLD;
            if (cutoff == 0) {
                Result r;
                r.sampleIndex = globalSampleIndex;
                r.refBegin = refBegin;
                r.score = static_cast<float>(coverage);
                r.cigar = "*";
                r.revComp = revCompFlag;
                queue.enqueue(r);
                return true;
            }
            // else return true; // get outa here

            // DUBUG
            // std::cerr << "in the twilight zone, coverage " << coverage << ", span "
            // << span << "\n";

            // Twilight: run SW and veto by cut-off

            std::string query = (revCompFlag ? reverseComplement(readSeq) : readSeq);
            const uint32_t PAD_LEFT = 2;
            const uint32_t PAD_RIGHT = 2;

            const uint32_t genomeLen = static_cast<uint32_t>(genomeStr.size());
            uint32_t winStart = (refBegin > PAD_LEFT) ? (refBegin - PAD_LEFT) : 0;
            uint64_t winEnd64 =
                static_cast<uint64_t>(refBegin) + static_cast<uint64_t>(query.size()) + PAD_RIGHT;
            uint32_t winEnd = (winEnd64 > genomeLen) ? genomeLen : static_cast<uint32_t>(winEnd64);

            if (winEnd <= winStart) {
                std::cerr << "in the twilight zone, winEnd " << winEnd << " <= winStart "
                          << winStart << "\n";
                continue;
            }

            std::string_view refWindow(genomeStr.data() + winStart, winEnd - winStart);

            //            std::string cigar;
            //            int sw_score = 0;
            //            uint32_t ref_begin_local = 0;
            //            uint32_t ref_end_local   = 0;
            //            alignSequences(refWindow, query, cigar, sw_score,
            //            ref_begin_local, ref_end_local);
            //	      const float eq = (float) countEqualsInCigar(cigar);
            //            const float al = (float) aligned_len_from_cigar(cigar);

#ifdef USE_EDLIB
            // Global alignment (NW) via edlib — forces full read alignment
            std::string cigar = "*";
            const float eq_frac = alignNW_identity(refWindow, query);
            const float eq = eq_frac * read_len;
            const float al = (float)read_len; // global alignment covers full read
#else
            // Original SW local alignment
            std::string cigar;
            int sw_score = 0;
            uint32_t ref_begin_local = 0;
            uint32_t ref_end_local = 0;
            alignSequences(refWindow, query, cigar, sw_score, ref_begin_local, ref_end_local);
            const float eq = (float)countEqualsInCigar(cigar);
            const float al = (float)aligned_len_from_cigar(cigar);
#endif

            //	    std::cerr << "@" << globalSampleIndex << " " << read_len <<
            //" "<< coverage << " " << span
            //		    << " " << al << " "  << eq/read_len << " " << eq/al
            //<< "\n";

            // Twilight veto threshold:
            // - independent of cutoff (so SW still works when cutoff==0)
            //	    std::cerr << "read_len " << read_len << ", alignment len "
            //<< al << ", identity " << eq << "\n";
            // Span + identity filter: reject if chain span is suspiciously short
            // AND full-read identity is low — smoking gun for bacterial contamination
            // fragments.
            if (span < SUSPICIOUS_SPAN && eq < read_len * IDENTITY_THRESHOLD)
                continue;

            if (eq / read_len < cutoff / 100.0) {
                continue;
            }

            // Accepted: still output coverage-based (SW was cleanup veto only)

            //            std::cerr << "in the twilight zone, coverage " << coverage
            //            << ", span " << span << "\n"; std::cerr << "Accept: Identity
            //            frac from read " << eq/read_len << " >= cutoff " <<
            //            cutoff/100.0 << "\n";

            Result r;
            r.sampleIndex = globalSampleIndex;
            r.refBegin = refBegin;
            r.score = static_cast<float>(coverage);
            r.cigar = cigar;
            r.revComp = revCompFlag;
            queue.enqueue(r);
            return true;
        }

        // ------------------------------------------------------------
        // MODE A/B: GET_SW == true
        //   - Always SW
        //   - Optional gate: if cutoff>0 then require score >= cutoff
        //   - Output SW-derived result
        // ------------------------------------------------------------
        {
            const uint32_t genomeLen = static_cast<uint32_t>(genomeStr.size());
            if (genomeLen == 0)
                continue;

            std::string query = (revCompFlag ? reverseComplement(readSeq) : readSeq);
            const uint32_t PAD_LEFT = 80;
            const uint32_t PAD_RIGHT = 80;

            uint32_t winStart = (refBegin > PAD_LEFT) ? (refBegin - PAD_LEFT) : 0;
            uint64_t winEnd64 =
                static_cast<uint64_t>(refBegin) + static_cast<uint64_t>(query.size()) + PAD_RIGHT;
            uint32_t winEnd = (winEnd64 > genomeLen) ? genomeLen : static_cast<uint32_t>(winEnd64);

            if (winEnd <= winStart)
                continue;

            std::string_view refWindow(genomeStr.data() + winStart, winEnd - winStart);

            std::string cigar;
            int sw_score = 0;
            uint32_t ref_begin_local = 0;
            uint32_t ref_end_local = 0;

            alignSequences(refWindow, query, cigar, sw_score, ref_begin_local, ref_end_local);
            if (cigar.empty())
                continue;

            const uint32_t eq = countEqualsInCigar(cigar);

#ifdef USE_ABSOLUTE_MATCH_COUNT
            // ============================================================
            // ORIGINAL VERSION: Absolute match count scoring
            // ============================================================

            // Gate: require at least 'cutoff' exact matching bases
            if (gate_sw_global && eq < static_cast<uint32_t>(cutoff))
                continue;

            Result r;
            r.sampleIndex = globalSampleIndex;
            r.refBegin = winStart + ref_begin_local;
            r.score = static_cast<float>(eq); // bp identity count
            r.cigar = cigar;
            r.revComp = revCompFlag;

            queue.enqueue(r);
            return true;

#else // USE_ALIGNMENT_DENSITY
      // ============================================================
      // NEW VERSION: Alignment density scoring
      // ============================================================

            // Compute alignment span (reference positions covered)
            const uint32_t alignment_span = computeAlignmentSpan(cigar);

            // Avoid division by zero
            if (alignment_span == 0)
                continue;

            // Compute density: fraction of aligned region that matches exactly
            const float density = static_cast<float>(eq) / static_cast<float>(alignment_span);

            // Three-way gate:
            // 1. Alignment must span substantial region
            // 2. Density must be high enough (not random matches)
            // 3. Absolute matches must exceed minimum (not just 20bp at 100%)
            const bool passes_span_gate = (alignment_span >= MIN_ALIGNMENT_SPAN);
            const bool passes_density_gate = (density >= MIN_ALIGNMENT_DENSITY);
            const bool passes_absolute_gate = (eq >= MIN_ABSOLUTE_MATCHES);

            if (!passes_span_gate || !passes_density_gate || !passes_absolute_gate)
                continue;

            Result r;
            r.sampleIndex = globalSampleIndex;
            r.refBegin = winStart + ref_begin_local;
            r.score = density * 100.0f; // Report as percentage (0-100)
            r.cigar = cigar;
            r.revComp = revCompFlag;

            queue.enqueue(r);
            return true;
#endif
        }
    }

    return false;
}
