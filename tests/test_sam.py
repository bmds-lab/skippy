"""End-to-end SAM check using an existing mitochondrial plain index.

Run after indexing: python3 -m tests.test_sam [path/to/mito]
"""

import pathlib
import re
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[1]


def main() -> None:
    genome = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "tests/data/mito"
    if not (genome / "plain/mito_refID.idx").exists():
        raise SystemExit(f"Build the plain mitochondrial index first: {genome}")

    lines = (ROOT / "tests/data/short_reads.fasta").read_text().splitlines()
    reads = dict(zip((line[1:] for line in lines[::2]), lines[1::2]))
    names = ("short_exact_01", "short_reverse_01", "short_negative_01")
    qualities = {name: "".join(chr(33 + i % 40) for i in range(len(reads[name])))
                 for name in names}

    with tempfile.TemporaryDirectory() as directory:
        work = pathlib.Path(directory)
        fastq = work / "reads.fastq"
        fastq.write_text("".join(
            f"@{name} description\n{reads[name]}\n+\n{qualities[name]}\n"
            for name in names
        ))
        result = subprocess.run(
            [str(ROOT / "skippy"), "search", "-genome", str(genome), "-stripe", "plain",
             "-technology", "illumina", "-data", str(fastq), "-out_dir", str(work),
             "-format", "sam", "-sw", "false", "-threads", "2"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        assert "compute Smith Waterman: TRUE" in result.stdout

        lines = (work / "reads.sam").read_text().splitlines()
        assert lines[:2] == ["@HD\tVN:1.6\tSO:unsorted",
                             "@SQ\tSN:NC_012920.1\tLN:16569"]
        records = [line.split("\t") for line in lines[2:]]
        assert [row[0] for row in records] == list(names)
        assert all(len(row) == 11 for row in records)
        forward, reverse, unmapped = records

        for row in (forward, reverse):
            assert row[2:5] == ["NC_012920.1", "1", "255"]
            operations = re.findall(r"(\d+)([MIDNSHP=X])", row[5])
            assert operations and "".join(f"{n}{op}" for n, op in operations) == row[5]
            assert sum(int(n) for n, op in operations if op in "MIS=X") == len(row[9])

        assert forward[1] == "0"
        assert forward[9:] == [reads[names[0]], qualities[names[0]]]
        assert reverse[1] == "16"
        complements = str.maketrans("ACGT", "TGCA")
        assert reverse[9:] == [reads[names[1]].translate(complements)[::-1],
                               qualities[names[1]][::-1]]
        assert unmapped[1:9] == ["4", "*", "0", "0", "*", "*", "0", "0"]
        assert unmapped[9:] == [reads[names[2]], qualities[names[2]]]

        blocker = work / "not_a_directory"
        blocker.write_text("")
        failed = subprocess.run(
            [str(ROOT / "skippy"), "search", "-genome", str(genome), "-stripe", "plain",
             "-technology", "illumina", "-data", str(fastq), "-out_dir", str(blocker),
             "-format", "sam", "-threads", "2"],
            capture_output=True, text=True, timeout=10,
        )
        assert failed.returncode == 1 and "Error opening file" in failed.stderr

    print("PASS: SAM records and output-error handling")


if __name__ == "__main__":
    main()
