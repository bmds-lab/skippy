# Baseline results

Captured from the current working tree on 2026-10-19. Both output files were
identical across two searches using the same plain index. Search uses
`-stripe plain`, `illumina` for short reads, and `ont` for long reads.
Smith-Waterman is off by default, so scores represent coverage and CIGAR is `*`.

| Queries | Detected |
| --- | ---: |
| Short exact | 19/20 |
| Short reverse complement | 20/20 |
| Short substitutions | 19/20 |
| Random controls | 0/20 false positives |
| Long exact | 9/10 |
| Long indels | 9/10 |

This is a known-failing baseline, not the desired result. Forward reads starting
at reference offset 0 have positive scores but reference `*` and position 0.
This affects `short_exact_01`, `short_substitutions_01`, `long_exact_01`, and
`long_indels_01`. The test still requires every exact query to be detected.

Keep these files when changing code. An output diff shows behavior changes,
not whether those changes are correct. Repeatability after rebuilding the
Monte Carlo seed matrices has not been checked.
