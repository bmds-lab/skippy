"""Fast checks for command-line help and failure statuses."""

import pathlib
import subprocess


SKIPPY = pathlib.Path(__file__).resolve().parents[1] / "skippy"


def check(args: list[str], status: int, message: str, stream: str) -> None:
    result = subprocess.run([str(SKIPPY), *args], capture_output=True, text=True)
    output = getattr(result, stream)
    assert result.returncode == status, (args, result.returncode, result.stderr)
    assert message in output, (args, output)


def main() -> None:
    for args, message in [
        (["--help"], "Usage: skippy"),
        (["index", "--help"], "skippy index -genome"),
        (["search", "--help"], "Usage: skippy search"),
        (["search", "--help", "all"], "Power user flags"),
        (["--version"], "skippy v1"),
        (["index", "--version"], "skippy v1"),
        (["search", "--version"], "skippy v1"),
    ]:
        check(args, 0, message, "stdout")

    for args, status, message in [
        ([], 2, "expected index or search"),
        (["unknown"], 2, "unknown command"),
        (["index"], 2, "specify -genome"),
        (["index", "-genome"], 2, "missing value for -genome"),
        (["index", "-genome", "-stripe", "plain"], 2, "missing value for -genome"),
        (["index", "-bogus", "value"], 2, "unknown index option"),
        (["index", "-stripe", "99999999999999999999"], 2,
         "weight exceeds maximum of 64"),
        (["index", "-stripe", "40X40"], 2,
         "weight exceeds maximum of 64"),
        (["index", "-stripe", "plain"], 2, "specify -genome"),
        (["index", "-reference", "/no-such-skippy-reference.fasta"], 1,
         "reference FASTA not found"),
        (["search"], 2, "missing search options"),
        (["search", "--bogus"], 2, "unknown search option"),
        (["search", "-bogus"], 2, "unknown search option"),
        (["search", "-threads", "oops"], 2, "invalid value for -threads"),
        (["search", "-threads", "0"], 2, "-threads must be at least 1"),
        (["search", "-gate", "-1"], 2, "-gate must be nonnegative"),
        (["search", "-threads", "-data"], 2, "missing value for -threads"),
        (["search", "-data", "reads.fasta"], 2, "missing -genome"),
        (["search", "-technology", "unknown"], 2, "unknown technology"),
        (["search", "-format", "bam"], 2, "-format must be out or sam"),
        (["search", "-data", "reads.fasta", "-genome", "/no-such-genome",
          "-format", "sam", "-paired", "true"], 2,
         "paired SAM output is not supported"),
    ]:
        check(args, status, message, "stderr")

    print("PASS: CLI help, version, and error statuses")


if __name__ == "__main__":
    main()
