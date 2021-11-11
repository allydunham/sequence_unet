#!/usr/bin/env python3
"""
Split a fasta file into single sequence files
"""
import argparse
from Bio import SeqIO

def main(args):
    """Main"""
    fasta = SeqIO.parse(args.fasta, "fasta")

    for record in fasta:
        name = record.id.split("|")[1] if args.uniprot else record.id
        with open(f"{args.outdir}/{name}.fa", "w") as fasta_file:
            print(record.format("fasta"), file=fasta_file, end="")

def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('fasta', metavar="F", help="Fasta file to split")

    parser.add_argument('--outdir', '-o', default=".", help="Output directory")
    parser.add_argument('--uniprot', '-u', action="store_true",
                        help="Extract Uniprot IDs from record metadata")

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())