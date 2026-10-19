// findSpans.cpp — skip covered seeds + coverage-based admission (and optional span reporting)
// CURRENT-RIG VERSION: bitmap union for same start + optional adjacent-start merge
// Implements EVERYTHING discussed:
//  1) "Dead means dead": once a chain is accepted, ALL participating seeds are removed from ANY future consideration
//     (not anchors, not members).
//  2) Safe feasibility pruning using GATE_COVERAGE: if even the upper bound can't reach the gate, abort this anchor-hit.
//  3) Stop early on a "very good" chain (>= 50% read coverage, and >= GATE_COVERAGE), since chains are disjoint sets.
//  4) Hybrid search: linear for short tails, binary for long tails (already present).
//  5) Cursor discipline: cursors are NOT reset per anchor occurrence; only after seed i is consumed or a chain is found.
//  6) Removes the "seed_covered over best_spans" loop (logically unnecessary under dead-means-dead and can block fragments).

#include "encode_kmer_2bit.h"
#include "findSpans.hpp"
#include "stripe_pattern.hpp"

#include <algorithm>
#include <bit>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <span>
#include <tuple>
#include <vector>

extern uint32_t GATE_COVERAGE;
extern StripePattern gPattern;


namespace
{
    struct CandidateWithBits
    {
        uint32_t coverage;
        uint32_t minBit;
        uint32_t maxBit;
        std::vector<uint64_t> words; // read-space coverage bitmap
    };

    static inline uint32_t abs_diff_u32(uint32_t a, uint32_t b)
    {
        return (a >= b) ? (a - b) : (b - a);
    }

    static std::vector<std::tuple<uint32_t, uint32_t, uint32_t>>
    findSpansPtrCore(std::vector<ListStructPtr> &seedLists,
                     std::vector<uint32_t> &pos,
                     uint32_t minMatches,
                     uint32_t range,
                     uint32_t slack)
    {
        const size_t n = seedLists.size();
        if (n == 0 || range == 0)
            return {};

        // ----------------------------------------------------
        // In-situ reordering: shortest occurrence lists first
        // ----------------------------------------------------
        std::vector<size_t> order(n);
        std::iota(order.begin(), order.end(), 0);

        std::sort(order.begin(), order.end(),
                  [&](size_t a, size_t b)
                  { return seedLists[a].length < seedLists[b].length; });

        for (size_t i = 0; i < n; ++i)
        {
            if (order[i] == i)
                continue;

            size_t curr = i;
            while (order[curr] != i)
            {
                size_t next = order[curr];
                std::swap(seedLists[curr], seedLists[next]);
                std::swap(pos[curr], pos[next]);
                order[curr] = curr;
                curr = next;
            }
            std::swap(seedLists[curr], seedLists[i]);
            std::swap(pos[curr], pos[i]);
            order[curr] = curr;
        }

        // Stripe-aware coverage bitmap in read-space.
        const uint32_t effectiveRange = range;
        const size_t nWords = (static_cast<size_t>(effectiveRange) + 63) >> 6;

        auto set_bit = [&](uint32_t bit,
                           std::vector<uint64_t> &words,
                           uint32_t &coverage,
                           uint32_t &minBit,
                           uint32_t &maxBit)
        {
            if (bit >= effectiveRange)
                return;
            const size_t wi = static_cast<size_t>(bit >> 6);
            const uint64_t mask = 1ULL << (bit & 63);
            if ((words[wi] & mask) == 0)
            {
                words[wi] |= mask;
                ++coverage;
                if (bit < minBit)
                    minBit = bit;
                if (bit > maxBit)
                    maxBit = bit;
            }
        };

        // Candidate starts: genome_start -> bitmap union candidate
        std::map<uint32_t, CandidateWithBits> best_spans;
        std::vector<size_t> cursor(n, 0);

        // DEAD MEANS DEAD: if used[s]==1, that seed is excluded from ALL future work.
        std::vector<uint8_t> used(n, 0);

        // === Anchor bound tightened by coverage feasibility ===
        const uint32_t minSeedsForCoverage =
            (GATE_COVERAGE + STRIPED_WEIGHT - 1) / STRIPED_WEIGHT; // ceil
        const uint32_t requiredSeeds = std::max<uint32_t>(minMatches, minSeedsForCoverage);

        const size_t max_anchor_i_exclusive =
            (requiredSeeds <= n) ? (n - requiredSeeds + 1) : 0;

        // Minimum partners needed beyond the anchor to reach GATE_COVERAGE.
        const uint32_t minPartnersNeeded =
            (GATE_COVERAGE > STRIPED_WEIGHT)
            ? ((GATE_COVERAGE - STRIPED_WEIGHT + STRIPED_WEIGHT - 1) / STRIPED_WEIGHT)
            : 0u;

        bool stop_all = false;

        for (size_t i = 0; i < max_anchor_i_exclusive; ++i)
        {
            if (used[i])
                continue;

            const auto &list = seedLists[i];
            if (list.start == nullptr || list.length == 0)
                continue;
            if (cursor[i] >= list.length)
                continue;

            // Dynamic anchor skip: count unused seeds after i.
            // If fewer than minPartnersNeeded remain, this anchor can never
            // reach GATE_COVERAGE regardless of what it matches.
            {
                uint32_t availableAfter = 0;
                for (size_t t = i + 1; t < n; ++t)
                    if (!used[t]) ++availableAfter;
                if (availableAfter < minPartnersNeeded)
                    break;  // no subsequent anchor can do better either
            }

            const uint32_t *list_start = list.start;

            // Anchor scan over occurrences in list i.
            for (size_t k = cursor[i]; k < list.length; ++k)
            {
                const uint32_t g0 = *(list_start + k);

                // Track which seeds participate in THIS accepted chain attempt.
                std::vector<size_t> chainSeeds;
                std::vector<uint32_t> chainGenomePos;
                chainSeeds.reserve(n);
                chainGenomePos.reserve(n);
                chainSeeds.push_back(i);
                chainGenomePos.push_back(g0);

                uint32_t hits = 1;

                // Remaining *eligible* seeds in suffix (upper bound for feasibility).
                uint32_t remainingSeeds = 0;
                for (size_t t = i + 1; t < n; ++t)
                    if (!used[t])
                        ++remainingSeeds;

                // Extend diagonal by scanning seeds j>i (sorted-by-length order).
                for (size_t j = i + 1; j < n; ++j)
                {
                    if (used[j])
                        continue;

                    // This seed is now being "considered" (either matches or not).
                    if (remainingSeeds > 0)
                        --remainingSeeds;

                    // ---- SAFE EARLY IMPOSSIBILITY PRUNE ----
                    // Cannot reach minMatches even if all remaining seeds matched.
                    if (hits + remainingSeeds < minMatches)
                        break;
                    // ----------------------------------------

                    const auto &next = seedLists[j];
                    if (next.start == nullptr || next.length == 0)
                        continue;

                    const uint32_t *base = next.start;
                    const uint32_t *end = base + next.length;

                    const uint32_t target = g0 + pos[j] - pos[i];
                    const uint32_t lo = (slack == 0) ? target : ((target >= slack) ? (target - slack) : 0u);
                    const uint32_t hi = (slack == 0) ? target : (target + slack);

                    const uint32_t *ptr = base + cursor[j];
                    if (ptr >= end)
                        continue;

                    const size_t remaining = static_cast<size_t>(end - ptr);

                    bool matched = false;

                    if (slack == 0)
                    {
                        // Exact match
                        if (remaining <= 16)
                        {
                            for (const uint32_t *p = ptr; p < end; ++p)
                            {
                                const uint32_t g = *p;
                                if (g > target)
                                {
                                    cursor[j] = static_cast<size_t>(p - base);
                                    break;
                                }
                                if (g == target)
                                {
                                    matched = true;
                                    cursor[j] = static_cast<size_t>(p - base) + 1;
                                    break;
                                }
                            }
                        }
                        else
                        {
                            const uint32_t *it = std::lower_bound(ptr, end, target);
                            cursor[j] = static_cast<size_t>(it - base);
                            if (it != end && *it == target)
                            {
                                matched = true;
                                cursor[j] += 1;
                            }
                        }
                    }
                    else
                    {
                        // Slack-aware: choose closest within [lo,hi]
                        if (remaining <= 16)
                        {
                            const uint32_t *p = ptr;
                            for (; p < end; ++p)
                            {
                                const uint32_t g = *p;
                                if (g >= lo)
                                    break;
                            }

                            if (p < end)
                            {
                                const uint32_t g = *p;
                                if (g <= hi)
                                {
                                    matched = true;
                                    cursor[j] = static_cast<size_t>(p - base) + 1;
                                }
                                // on failure: conservative advance to lo-slack
                                // (safe for future anchors whose target may be up to slack lower)
                            }
                            else
                            {
                                // p >= end and no match: conservative advance
                                const uint32_t safe_lo = (lo >= slack) ? (lo - slack) : 0u;
                                const uint32_t *safe_it = std::lower_bound(base + cursor[j], end, safe_lo);
                                cursor[j] = static_cast<size_t>(safe_it - base);
                            }
                        }
                        else
                        {
                            const uint32_t *it = std::lower_bound(ptr, end, target);

                            const uint32_t *best = nullptr;
                            uint32_t best_dist = 0;

                            if (it != end)
                            {
                                const uint32_t g = *it;
                                if (g >= lo && g <= hi)
                                {
                                    best = it;
                                    best_dist = abs_diff_u32(g, target);
                                }
                            }

                            if (it != ptr)
                            {
                                const uint32_t *pit = it - 1;
                                const uint32_t g = *pit;
                                if (g >= lo && g <= hi)
                                {
                                    const uint32_t d = abs_diff_u32(g, target);
                                    if (best == nullptr || d < best_dist)
                                    {
                                        best = pit;
                                        best_dist = d;
                                    }
                                }
                            }

                            if (best != nullptr)
                            {
                                matched = true;
                                cursor[j] = static_cast<size_t>(best - base) + 1;
                            }
                            else
                            {
                                // conservative advance: lo - slack
                                // safe for future anchors whose target may be up to slack lower
                                const uint32_t safe_lo = (lo >= slack) ? (lo - slack) : 0u;
                                const uint32_t *safe_it = std::lower_bound(ptr, end, safe_lo);
                                cursor[j] = static_cast<size_t>(safe_it - base);
                            }
                        }
                    }

                    if (matched)
                    {
                        ++hits;
                        chainSeeds.push_back(j);
                        const uint32_t gj = *(base + cursor[j] - 1);
                        chainGenomePos.push_back(gj);
                        // Coverage will be recomputed in position order at acceptance
                    }
                }

                if (hits >= minMatches)
                {
                    // --------------------------------------------------------
                    // SECOND PASS: recompute coverage in read-position order
                    // to correctly propagate cumulative indel shifts.
                    // chainSeeds were collected in inverted-list-length order
                    // (not position order), so without this pass coverage is
                    // incorrect for indel reads.
                    // --------------------------------------------------------

                    // Sort chain members by read position
                    std::vector<size_t> sortedIdx(chainSeeds.size());
                    std::iota(sortedIdx.begin(), sortedIdx.end(), 0);
                    std::sort(sortedIdx.begin(), sortedIdx.end(),
                              [&](size_t a, size_t b2) {
                                  return pos[chainSeeds[a]] < pos[chainSeeds[b2]];
                              });

                    // Compute coverage bitmap in position order, tracking cumulative shift
                    std::vector<uint64_t> coveredWords(nWords, 0);
                    uint32_t coverage = 0;
                    uint32_t minBit = std::numeric_limits<uint32_t>::max();
                    uint32_t maxBit = 0;
                    int32_t cumShift = 0;

                    for (size_t si = 0; si < sortedIdx.size(); ++si)
                    {
                        const size_t ci   = sortedIdx[si];
                        const size_t seed = chainSeeds[ci];
                        const uint32_t gj = chainGenomePos[ci];
                        const uint32_t pj = pos[seed];

                        // Cumulative shift at this position
                        const uint32_t expected = g0 + pj - pos[i];
                        cumShift = (int32_t)gj - (int32_t)expected;

                        // Coverage in anchor coordinate system
                        for (uint32_t b = 0; b < STRIPED_WEIGHT; ++b)
                        {
                            const int32_t bit = (int32_t)(pj + STRIPED_USED_INDICES[b]) - cumShift;
                            if (bit >= 0 && bit < (int32_t)effectiveRange)
                                set_bit((uint32_t)bit, coveredWords, coverage, minBit, maxBit);
                        }
                    }

                    const uint32_t genome_start = (g0 >= pos[i]) ? (g0 - pos[i]) : 0;

                    auto it = best_spans.find(genome_start);
                    if (it == best_spans.end())
                    {
                        best_spans.emplace(genome_start,
                                           CandidateWithBits{coverage, minBit, maxBit, std::move(coveredWords)});
                    }
                    else
                    {
                        CandidateWithBits &cand = it->second;
                        if (cand.words.size() != coveredWords.size())
                            cand.words.resize(coveredWords.size(), 0);
                        uint32_t added = 0;
                        for (size_t w = 0; w < coveredWords.size(); ++w)
                        {
                            const uint64_t new_bits = coveredWords[w] & ~cand.words[w];
                            if (new_bits)
                                added += static_cast<uint32_t>(__builtin_popcountll(new_bits));
                            cand.words[w] |= coveredWords[w];
                        }
                        cand.coverage += added;
                        if (minBit < cand.minBit) cand.minBit = minBit;
                        if (maxBit > cand.maxBit) cand.maxBit = maxBit;
                    }

                    // DEAD MEANS DEAD
                    for (size_t s : chainSeeds)
                    {
                        used[s] = 1;
                        cursor[s] = seedLists[s].length;
                    }
#if 0
                    // Stop early on a very good chain
                    if (coverage >= GATE_COVERAGE + STRIPED_WEIGHT)
                    {
                        stop_all = true;
                        goto done_search;
                    }
#endif
                    // Stop early only if the chain is both well-covered AND
                    // wide enough to clear the SW twilight zone fast-accept
                    // threshold (span > 79 in smithWaterman.cpp line 330).
                    // GATE_COVERAGE + STRIPED_WEIGHT = 33 + 16 = 49 bases
                    // ensures at least 4 surviving seeds contributed.
                    // span > 79 ensures SW will fast-accept without veto.
                    // Both conditions must hold to avoid stopping on a tight
                    // chain that SW would reject in the twilight zone.
                    if (coverage >= GATE_COVERAGE + STRIPED_WEIGHT &&
                        (maxBit - minBit) > 79)
                    {
                        stop_all = true;
                        goto done_search;
                    }

                    break;
                }
            }

            // Reset cursors for the suffix after seed i attempt, but only for seeds still alive.
            std::fill(cursor.begin() + i + 1, cursor.end(), 0);

            // Anchor i is exhausted either way.
            cursor[i] = seedLists[i].length;

            if (stop_all)
                break;
        }

    done_search:

        // ------------------------------------------------------------

        const bool MERGE_ADJ_STARTS = true;
        const uint32_t MERGE_START_DIST = 15; // 20

        if (MERGE_ADJ_STARTS && best_spans.size() >= 2)
        {
            auto it = best_spans.begin();
            auto prev = it++;
            while (it != best_spans.end())
            {
                const uint32_t s0 = prev->first;
                const uint32_t s1 = it->first;

                if (s1 - s0 <= MERGE_START_DIST)
                {
                    CandidateWithBits &A = prev->second;
                    CandidateWithBits &B = it->second;

                    if (A.words.size() != B.words.size())
                    {
                        const size_t m = std::max(A.words.size(), B.words.size());
                        A.words.resize(m, 0);
                        B.words.resize(m, 0);
                    }

                    uint32_t added = 0;
                    for (size_t w = 0; w < A.words.size(); ++w)
                    {
                        const uint64_t new_bits = B.words[w] & ~A.words[w];
                        if (new_bits)
                            added += static_cast<uint32_t>(__builtin_popcountll(new_bits));
                        A.words[w] |= B.words[w];
                    }

                    A.coverage += added;
                    if (B.minBit < A.minBit)
                        A.minBit = B.minBit;
                    if (B.maxBit > A.maxBit)
                        A.maxBit = B.maxBit;

                    it = best_spans.erase(it);
                    continue;
                }

                prev = it++;
            }
        }

        // === COVERAGE ADMISSION (still no span gate here) ===
        std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> spans;
        spans.reserve(best_spans.size());

#ifdef FAKE_COVERAGE_PASSTHROUGH
        // DEBUG: push ALL candidates through with fake coverage=999
        // so SW sees every read regardless of actual chain quality.
        if (best_spans.empty())
        {
            // No chain at all — emit a dummy entry at position 0
            spans.emplace_back(0u, 999u, 999u);
        }
        else
        {
            for (const auto &[start, cand] : best_spans)
            {
                const uint32_t span =
                    (cand.minBit <= cand.maxBit) ? (cand.maxBit - cand.minBit + 1) : 0u;
                spans.emplace_back(start, 999u, span);
            }
        }
#else
        for (const auto &[start, cand] : best_spans)
        {
            if (cand.coverage < GATE_COVERAGE)
                continue;

            const uint32_t span =
                (cand.minBit <= cand.maxBit) ? (cand.maxBit - cand.minBit + 1) : 0u;

            spans.emplace_back(start, cand.coverage, span);
        }
#endif

        std::sort(spans.begin(), spans.end(),
                  [](const auto &a, const auto &b)
                  {
                      if (std::get<1>(a) != std::get<1>(b))
                          return std::get<1>(a) > std::get<1>(b);
                      return std::get<0>(a) < std::get<0>(b);
                  });

        return spans;
    }
} // namespace

// ----------------------------------------------------------------------------
// Existing API (unchanged): returns (refBegin, coverage)
// ----------------------------------------------------------------------------
std::vector<std::tuple<uint32_t, uint32_t>>
findSpansPtr(std::vector<ListStructPtr> &seedLists,
             std::vector<uint32_t> &pos,
             uint32_t minMatches,
             uint32_t range,
             uint32_t slack)
{
    auto spans3 = findSpansPtrCore(seedLists, pos, minMatches, range, slack);

    std::vector<std::tuple<uint32_t, uint32_t>> spans2;
    spans2.reserve(spans3.size());

    for (const auto &t : spans3)
        spans2.emplace_back(std::get<0>(t), std::get<1>(t));

    return spans2;
}

// ----------------------------------------------------------------------------
// New API: returns (refBegin, coverage, span)
// ----------------------------------------------------------------------------
std::vector<std::tuple<uint32_t, uint32_t, uint32_t>>
findSpansPtrWithSpan(std::vector<ListStructPtr> &seedLists,
                     std::vector<uint32_t> &pos,
                     uint32_t minMatches,
                     uint32_t range,
                     uint32_t slack)
{
    return findSpansPtrCore(seedLists, pos, minMatches, range, slack);
}
