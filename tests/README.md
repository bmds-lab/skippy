# Mitochondrial smoke test

Requires Python 3, make, and the C++ compiler. No pip packages or downloads.

`make test` includes index boundary checks on small references and failed
seed-matrix writes. The SAM test also searches 120 files with a low open-file limit.

```bash
python3 -m tests.mito        # build, index, search, check
python3 -m tests.mito search # rebuild skippy, reuse index, check
python3 -m tests.mito check  # check existing results only
python3 -m tests.test_mito   # fast checker self-test
python3 -m tests.test_sam    # SAM check after the plain index exists
```

The included human mitochondrial reference is NC_012920.1, 16,569 bases.
The 100 fixed queries cover exact matches, reverse complements, substitutions,
indels, and random controls. All exact queries must match. The test reports
mutated-read detection and false positives without pass/fail thresholds.
The `.out` checker does not check alignment coordinates, strand accuracy, or
circular boundaries. `tests.test_sam` checks SAM coordinates, strand, CIGAR,
qualities, and unmapped records on three reads. Pass a genome directory as its
argument to reuse an index outside this repository.

Logs and search output go in `tests/results/`. Seed-matrix generation can be slow.
The runner refuses an existing `tests/data/mito/plain` index. Move it aside to
rebuild. Source details and query coordinates are in `data/provenance.json` and
`data/truth.json`. Truth coordinates start at 0; output read numbers start at 1.

Saved per-read outputs are in `baseline/`. The current baseline has known
failures described in `baseline/README.md`. After rerunning search, compare:

```bash
diff -u tests/baseline/short_reads.out tests/results/short/short_reads.out
diff -u tests/baseline/long_reads.out tests/results/long/long_reads.out
```
