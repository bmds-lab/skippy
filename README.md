# Skippy

A Statistically Guaranteed High-Sensitivity Read Classifier

K-mer index and search over large reference genomes.

## Build

```bash
make
```

Produces one executable, `./skippy`, with `index` and `search` subcommands.

Use `./skippy --help`, `./skippy index --help`, or `./skippy search --help`
for command options. `--version` works at each level. Commands exit with 0 on
success, 2 for invalid command-line arguments, and 1 for runtime failures.

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

## Formatting

The root `.clang-format` follows the C++ Core Guidelines' suggested Stroustrup
brace layout and C++ pointer/reference layout. Four-space indentation and
100-column wrapping are project choices. Use clang-format 23.1.0 (the version
checked on GitHub pushes and pull requests). Format a file from the repository
root with:

```bash
clang-format -i indexer/main.cpp
```

To format all tracked C and C++ files:

```bash
git ls-files -z '*.c' '*.cpp' '*.h' '*.hpp' | xargs -0 clang-format -i
```

To check formatting without changing files, replace `-i` with
`--dry-run --Werror`.

## Small end-to-end test

See [tests/README.md](tests/README.md) for a human mitochondrial reference,
reproducible query sequences, and an automated smoke test:

```bash
python3 -m tests.mito
```
