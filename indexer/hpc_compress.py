#!/usr/bin/env python3
"""
hpc_compress.py — Homopolymer compress a FASTA file.
Usage: python3 hpc_compress.py input.fasta output.fasta
"""
import sys

if len(sys.argv) != 3:
    print("Usage: hpc_compress.py input.fasta output.fasta")
    sys.exit(1)

infile, outfile = sys.argv[1], sys.argv[2]

with open(infile) as f, open(outfile, 'w') as out:
    for line in f:
        if line.startswith('>'):
            out.write(line)
        else:
            seq = line.rstrip()
            if not seq:
                continue
            result = [seq[0]]
            for c in seq[1:]:
                if c != result[-1]:
                    result.append(c)
            out.write(''.join(result) + '\n')

print(f"Done: {infile} -> {outfile}")
