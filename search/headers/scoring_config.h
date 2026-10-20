// scoring_config.h
// =============================================================================
// SKIPPY SCORING VERSION CONTROL
// =============================================================================
// Edit this ONE file to switch between scoring methods.
// Both main.cpp and smithWaterman.cpp include this header.
//
// After changing, recompile:
//   make clean && make
// =============================================================================

#ifndef SCORING_CONFIG_H
#define SCORING_CONFIG_H

// -----------------------------------------------------------------------------
// UNCOMMENT ONE OF THE FOLLOWING:
// -----------------------------------------------------------------------------

// #define USE_ABSOLUTE_MATCH_COUNT    // Version 1: Count of exact matches
#define USE_ALIGNMENT_DENSITY // Version 2: Density-based scoring

// -----------------------------------------------------------------------------
// DO NOT EDIT BELOW THIS LINE
// -----------------------------------------------------------------------------

// Sanity check: ensure exactly one version is defined
#if defined(USE_ABSOLUTE_MATCH_COUNT) && defined(USE_ALIGNMENT_DENSITY)
#error "ERROR: Both scoring methods defined! Uncomment only ONE."
#endif

#if !defined(USE_ABSOLUTE_MATCH_COUNT) && !defined(USE_ALIGNMENT_DENSITY)
#warning "No scoring method defined, defaulting to USE_ABSOLUTE_MATCH_COUNT"
#define USE_ABSOLUTE_MATCH_COUNT
#endif

// Version identifiers for runtime display
#ifdef USE_ALIGNMENT_DENSITY
#define SCORING_VERSION "Density-Based (v2)"
#else
#define SCORING_VERSION "Absolute Count (v1)"
#endif

#endif // SCORING_CONFIG_H
