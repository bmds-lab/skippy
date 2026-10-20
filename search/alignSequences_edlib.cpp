// alignSequences_edlib.cpp
//
// Drop-in overload of alignSequences() using edlib global alignment (NW).
//
// Edlib global mode forces alignment across the full query length against
// the reference window — dragging down the score for bacterial reads that
// only align locally, unlike SW which finds the best local hit.
//
// The function signature matches the existing alignSequences() in
// smithWaterman.cpp exactly, so it can be called alongside it.
//
// To use: add -DUSE_EDLIB to CFLAGS and include this file in the build.
// The caller in smithWaterman.cpp selects between SW and edlib via #ifdef.
//
// Build addition in GNUmakefile:
//   Add edlib-1.2.7/edlib/src/edlib.cpp to the object list
//   Add -Iedlib-1.2.7/edlib/include to CFLAGS

#include <cstdint>
#include <string>
#include <string_view>

#include "edlib-1.2.7/edlib/include/edlib.h"

// ---------------------------------------------------------------------------
// alignSequencesNW
//
// Global alignment (Needleman-Wunsch style) using edlib.
// Forces alignment of the full query against the reference window.
//
// Parameters match alignSequences() in smithWaterman.cpp:
//   refWindow     — reference substring to align against
//   query         — read sequence (forward or RC as needed)
//   eq            — output: number of exact matching bases (= operations)
//   read_len      — output: length of query (for eq/read_len calculation)
//   edit_distance — output: edit distance (substitutions + indels)
// ---------------------------------------------------------------------------

void alignSequencesNW(const std::string_view& refWindow, const std::string& query, uint32_t& eq,
                      uint32_t& read_len, int& edit_distance)
{
    eq = 0;
    edit_distance = -1;
    read_len = static_cast<uint32_t>(query.size());

    if (query.empty() || refWindow.empty())
        return;

    // edlib global mode: align full query against reference window
    EdlibAlignResult result =
        edlibAlign(query.c_str(), static_cast<int>(query.size()), refWindow.data(),
                   static_cast<int>(refWindow.size()),
                   edlibNewAlignConfig(-1,              // no edit distance limit
                                       EDLIB_MODE_NW,   // global alignment (Needleman-Wunsch)
                                       EDLIB_TASK_PATH, // compute full alignment path
                                       nullptr, 0));    // no custom alphabet

    if (result.status != EDLIB_STATUS_OK || result.editDistance < 0) {
        edlibFreeAlignResult(result);
        return;
    }

    edit_distance = result.editDistance;

    // Count exact matches from the alignment path
    // edlib path operations:
    //   0 = match
    //   1 = insert into target (deletion from query)
    //   2 = insert into query  (insertion in query)
    //   3 = mismatch
    if (result.alignment != nullptr) {
        for (int i = 0; i < result.alignmentLength; i++) {
            if (result.alignment[i] == 0) // exact match
                eq++;
        }
    }

    edlibFreeAlignResult(result);
}

// ---------------------------------------------------------------------------
// Convenience wrapper matching the SW interface more closely.
// Returns eq/read_len as a float — same quantity used in the SW veto.
// ---------------------------------------------------------------------------

float alignNW_identity(const std::string_view& refWindow, const std::string& query)
{
    uint32_t eq = 0;
    uint32_t read_len = 0;
    int edit_dist = -1;

    alignSequencesNW(refWindow, query, eq, read_len, edit_dist);

    if (read_len == 0)
        return 0.0f;
    return static_cast<float>(eq) / static_cast<float>(read_len);
}
