#!/usr/bin/env python3
import sys
from Bio import SeqIO

def main(file_name, min_length = None, max_length = None):
    filters = []
    if min_length:
        filters.append(lambda x: len(x) >= min_length)
    if max_length:
        filters.append(lambda x: len(x) <= max_length)
    if not filters:
        raise ValueError("No filters specified")

    with open(file_name, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if all(f(record) for f in filters):
                print(f">{record.id}\n{record.seq}")

if __name__ == '__main__':
    main(*sys.argv[1:])