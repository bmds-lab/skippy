# Skippy

A Statistically Guaranteed High-Sensitivity Read Classifier

K-mer index and search over large reference genomes.

## Build

```bash
make
```

Produces one executable, `./skippy`, with `index` and `search` subcommands.

## Usage

Both tools use a genome-directory layout. A genome root holds
`<name>/<name>.fasta` and one subdirectory per stripe pattern:

```text
GRCh38/
  GRCh38.fasta
  plain/              # index artefacts for the plain 16-mer stripe
  illumina -> plain   # technology symlinks (created by skippy index)
  ont -> plain
```

Index a reference:

```bash
./skippy index -genome /path/to/GRCh38 -stripe plain
```

Search:

```bash
./skippy search -genome /path/to/GRCh38 -technology illumina \
  -data reads.fasta -out_dir results
```

`-technology` selects the contract (`illumina`/`hifi` → subs;
`pacbio`/`nanopore`/`ont` → indels) and a default error rate. Override with
`-error_target`. See `./skippy search --help all` for the full flag set.

## Requirements

- g++ with C++20
- pthread

## Small end-to-end test

See [tests/README.md](tests/README.md) for a human mitochondrial reference,
reproducible query sequences, and an automated smoke test:

```bash
python3 -m tests.mito
```
