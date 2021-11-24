#!/usr/bin/env python3
"""
Filter a Fasta file to only include records based on a list of IDs
"""
import sys
import argparse
import gzip
from functools import partial
from Bio import SeqIO

from sequence_unet.predict import AMINO_ACIDS

def main(args):
    """
    Main
    """
    args = arg_parser().parse_args()

    with open(args.ids, mode="r") as id_file:
        ids = set(i.strip() for i in id_file)

    _open = partial(gzip.open, mode="rt") if args.gzip else open
    with _open(args.fasta) as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            ident = record.id.split("|")[1] if args.uniprot else record.id

            if args.iupac and not all(i in AMINO_ACIDS for i in record):
                continue

            if ident in ids:
                print(record.format("fasta"), file=sys.stdout, end="")

def arg_parser():
    """Argument parser"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('fasta', metavar='F', help="Fasta file")
    parser.add_argument('ids', metavar='I', help="File containing a list of IDs")

    parser.add_argument('--gzip', '-g', action="store_true", help="Fasta file is gzipped")
    parser.add_argument('--uniprot', '-u', action="store_true",
                        help="Fasta file is Uniprot formatted and ID list contains Uniprot IDs")
    parser.add_argument('--iupac', '-i', action="store_true",
                        help="Filter to only include cannonical amino acids")

    return parser

if __name__ == "__main__":
    main()
