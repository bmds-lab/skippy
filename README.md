# Skippy

A Statistically Guaranteed High-Sensitivity Read Classifier

K-mer index and search over large reference genomes (indexer + search tools).

## Build

**Indexer** (builds a reference index from FASTA):

```bash
cd indexer
make
```

Produces `indexReference`.

**Search** (maps reads against a built index):

```bash
cd search
make -f GNUmakefile
```

Produces `skippy`.

## Usage

Both tools use a genome-directory layout. A genome root holds
`<name>/<name>.fasta` and one subdirectory per stripe pattern:

```text
GRCh38/
  GRCh38.fasta
  plain/              # index artefacts for the plain 16-mer stripe
  illumina -> plain   # technology symlinks (created by indexReference)
  ont -> plain
```

Index a reference:

```bash
./indexer/indexReference -genome /path/to/GRCh38 -stripe plain
```

Search:

```bash
./search/skippy -genome /path/to/GRCh38 -technology illumina \
  -data reads.fasta -out_dir results
```

`-technology` selects the contract (`illumina`/`hifi` → subs;
`pacbio`/`nanopore`/`ont` → indels) and a default error rate. Override with
`-error_target`. See `skippy --help all` for the full flag set.

## Requirements

- g++ with C++17 (indexer) and C++20 (search)
- pthread

## Small end-to-end test

See [tests/README.md](tests/README.md) for a human mitochondrial reference,
reproducible query sequences, and an automated smoke test:

```bash
python3 -m tests.mito
```
