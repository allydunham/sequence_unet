#!/usr/bin/env python3
"""
Split a fasta file into sub-files
"""
import argparse
import math
from pathlib import Path
from Bio import SeqIO

def main():
    """Main"""
    args = parse_args()

    root_name = Path(args.fasta).stem
    fasta = SeqIO.parse(args.fasta, "fasta")

    if args.files:
        max_seqs = math.ceil(sum(1 for _ in fasta) / args.files)
        fasta = SeqIO.parse(args.fasta, "fasta")

    elif args.seqs:
        max_seqs = args.seqs

    else:
        raise ValueError("Neither files seqs arguments passed")

    seqs_left = max_seqs
    file_number = 0
    fasta_file = open(f"{args.outdir}/{root_name}_{file_number}.fa", "w")
    try:
        for record in fasta:
            if seqs_left == 0:
                fasta_file.close()
                file_number += 1
                fasta_file = open(f"{args.outdir}/{root_name}_{file_number}.fa", "w")
                seqs_left = max_seqs

            print(record.format("fasta"), file=fasta_file, end="")
            seqs_left -= 1
        else:
            fasta_file.close()
    except Exception as err:
        fasta_file.close()
        raise err

def arg_parser():
    """Argument parser"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('fasta', metavar="F", help="Fasta file to split")

    parser.add_argument('--outdir', '-o', default=".", help="Output directory")
    parser.add_argument('--files', '-n', default=0, type=int, help="Number of files to split into")
    parser.add_argument('--seqs', '-s', default=0, type=int, help="Number of sequences per file")

    return parser

def parse_args():
    """Process arguments"""
    parser = arg_parser()
    args = parser.parse_args()

    if bool(args.files) == bool(args.seqs):
        raise ValueError("Use --files or --seqs not both")

    return args

if __name__ == "__main__":
    main()