"""Self-checks for the mitochondrial result checker."""

import contextlib
import io
import pathlib
import sys
import tempfile

from tests import mito


def _expect_error(
    queries: pathlib.Path, output: pathlib.Path, message: str
) -> None:
    """Checks that invalid input raises ValueError with the expected message."""
    try:
        mito.check(queries, output)
    except ValueError as error:
        assert message in str(error), str(error)
    else:
        raise AssertionError(
            f'Checker accepted invalid input: {output.read_text()!r}'
        )


def main() -> None:
    """Checks valid results and rejection of malformed or missing results."""
    with tempfile.TemporaryDirectory() as directory:
        with contextlib.redirect_stdout(io.StringIO()):
            _check_results(pathlib.Path(directory))
    print('PASS: checker self-test')


def _check_results(directory: pathlib.Path) -> None:
    """Exercises the checker using temporary FASTA and output files."""
    queries = directory / 'reads.fasta'
    output = directory / 'reads.out'
    queries.write_text(
        '>short_exact_01\nACGT\n>short_reverse_01\nACGT\n'
        '>short_negative_01\nAAAA\n'
    )
    valid = '> 1 mito 1 100.00 4M\n< 2 mito 1 100.00 4M\n'
    output.write_text(valid)
    mito.check(queries, output)
    for bad in [
        'garbage', '> 4 mito 1 100.00 4M', '> 1 mito 1 100.00 4M',
        '> 01 mito 1 100.00 4M', '> 3 mito bad 100.00 4M',
        '> 3 mito 1 100oops 4M',
    ]:
        output.write_text(valid + bad + '\n')
        _expect_error(queries, output, '')
    for unmapped in [
        '', '> 1 mito 0 100.00 4M\n', '> 1 mito 1 0.00 4M\n',
        '= 1 mito 1 100.00 4M\n', '> 1 * 1 100.00 4M\n',
    ]:
        output.write_text('< 2 mito 1 100.00 4M\n' + unmapped)
        _expect_error(queries, output, 'Missing exact queries')
    output.write_text('')
    _expect_error(queries, output, 'Missing exact queries')
    for fixture in [
        '', '>short_typo_01\nACGT\n', '>short_negative_01\nACGT\n'
    ]:
        queries.write_text(fixture)
        _expect_error(queries, output, 'Invalid query categories')


if __name__ == '__main__':
    sys.exit(main())
