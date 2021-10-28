#!/usr/bin/env python3
"""
Filter a Fasta file to only include records based on a list of IDs
"""
import sys
import argparse
import gzip
from functools import partial
from Bio import SeqIO

def main(args):
    """
    Main
    """
    args = parse_args()

    with open(args.ids, mode="r") as id_file:
        ids = set(i.strip() for i in id_file)

    _open = partial(gzip.open, mode="rt") if args.gzip else open
    with _open(args.fasta) as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            ident = record.id.split("|")[1] if args.uniprot else record.id
            if ident in ids:
                print(record.format("fasta"), file=sys.stdout, end="")

def parse_args():
    """
    Parse comman line arguments
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('fasta', metavar='F', help="Fasta file")
    parser.add_argument('ids', metavar='I', help="File containing a list of IDs")

    parser.add_argument('--gzip', '-g', action="store_true", help="Fasta file is gzipped")
    parser.add_argument('--uniprot', '-u', action="store_true",
                        help="Fasta file is Uniprot formatted and ID list contains Uniprot IDs")

    return parser.parse_args()

if __name__ == "__main__":
    main()
