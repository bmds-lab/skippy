#!/usr/bin/env python3
"""Run the fixed mitochondrial smoke test. No third-party packages."""

import argparse
import collections
import pathlib
import re
import subprocess
import sys
from typing import Sequence, Union

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / 'tests/data'
RESULTS = ROOT / 'tests/results'


def check(queries: pathlib.Path, output: pathlib.Path) -> None:
    """Checks exact-read detection and prints counts for each query category.

    Args:
        queries: Fixed FASTA query file, in search output read-number order.
        output: Search output containing one row per reported alignment.

    Raises:
        OSError: A file cannot be read.
        ValueError: Categories or results are invalid, or exact reads are
            missing.
    """
    names = [
        line[1:] for line in queries.read_text().splitlines()
        if line.startswith('>')
    ]
    kinds = [name.split('_')[1] if '_' in name else '' for name in names]
    if 'exact' not in kinds or set(kinds) - {
        'exact', 'reverse', 'substitutions', 'negative', 'indels'
    }:
        raise ValueError(f'Invalid query categories in {queries}')
    seen, hits = set(), set()
    for line in output.read_text().splitlines():
        match = re.fullmatch(
            r'([><=])\s+(\d+)\s+(\S+)\s+(\d+)\s+'
            r'(-?\d+(?:\.\d+)?)\s+(\S+)\s*',
            line,
        )
        if not match:
            raise ValueError(f'Malformed result: {line}')
        direction, read_number, reference, position, score, _ = match.groups()
        number = int(read_number)
        if not 1 <= number <= len(names) or number in seen:
            raise ValueError(f'Duplicate or unknown read number: {number}')
        seen.add(number)
        if (
            direction != '=' and reference != '*'
            and int(position) > 0 and float(score) > 0
        ):
            hits.add(number)
    counts = collections.Counter(kinds)
    found = collections.Counter(kinds[number - 1] for number in hits)
    for kind, total in counts.items():
        label = 'false positives' if kind == 'negative' else 'detected'
        print(f'{queries.stem}/{kind}: {found[kind]}/{total} {label}')
    missing = [
        name for number, (name, kind) in enumerate(zip(names, kinds), 1)
        if kind in ('exact', 'reverse') and number not in hits
    ]
    if missing:
        raise ValueError('Missing exact queries: ' + ', '.join(missing))
    print('PASS: all exact and reverse-complement queries detected')


def _logged(
    command: Sequence[Union[str, pathlib.Path]], log: pathlib.Path
) -> None:
    """Runs a command, logging output and reporting its tail on failure."""
    print(f'Running {command[0]}; log: {log}', flush=True)
    with log.open('w') as handle:
        result = subprocess.run(
            [str(part) for part in command], cwd=ROOT,
            stdout=handle, stderr=subprocess.STDOUT, check=False,
        )
    if result.returncode:
        tail = '\n'.join(log.read_text().splitlines()[-20:])
        raise RuntimeError(f'{command[0]} failed:\n{tail}')


def run() -> None:
    """Builds both tools, creates a fresh index, and searches the fixed queries.

    Raises:
        ValueError: An index directory already exists.
        OSError: A command cannot start or a file cannot be accessed.
        RuntimeError: Indexing or search fails.
        subprocess.CalledProcessError: Compilation fails.
    """
    genome = DATA / 'mito'
    if (genome / 'plain').exists():
        raise ValueError(
            'Index already exists. Move tests/data/mito/plain aside '
            'before rebuilding.'
        )
    subprocess.run(['make', '-C', str(ROOT / 'indexer')], check=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    _logged(
        [ROOT / 'indexer/indexReference', '-genome', genome,
         '-stripe', 'plain'],
        RESULTS / 'index.log',
    )
    search()


def search() -> None:
    """Builds search and runs the queries using a completed plain index."""
    subprocess.run(
        ['make', '-C', str(ROOT / 'search'), '-f', 'GNUmakefile'], check=True
    )
    for group, technology in [('short', 'illumina'), ('long', 'ont')]:
        output = RESULTS / group
        output.mkdir(parents=True, exist_ok=True)
        _logged(
            [ROOT / 'search/skippy', '-genome', DATA / 'mito',
             '-stripe', 'plain', '-technology', technology,
             '-data', DATA / f'{group}_reads.fasta',
             '-out_dir', output],
            RESULTS / f'{group}.log',
        )


def main() -> None:
    """Runs the smoke test, reuses its index, or checks existing results."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'action', nargs='?', default='run', choices=['run', 'search', 'check']
    )
    args = parser.parse_args()
    try:
        if args.action == 'run':
            run()
        elif args.action == 'search':
            search()
        for group in ('short', 'long'):
            check(
                DATA / f'{group}_reads.fasta',
                RESULTS / group / f'{group}_reads.out',
            )
    except (
        OSError, ValueError, RuntimeError, subprocess.CalledProcessError
    ) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    sys.exit(main())
