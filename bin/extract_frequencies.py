#!/usr/bin/env python3
"""
Extract a table of PSSM frequencies from a ProteinNet file
"""
import argparse
import sys

from proteinnetpy.data import ProteinNetDataset
from pn_maps import AMINO_ACIDS

def main(args):
    """
    Run the Frequency UNET model on all records of a ProteinNet file
    """
    data = ProteinNetDataset(path=args.proteinnet, preload=False)

    print("id", "position", "wt", *AMINO_ACIDS, sep="\t", file=sys.stdout)
    for record in data:
        for i, (wt, pssm) in enumerate(zip(record.primary, record.evolutionary.T)):
            print(record.id, i + 1, wt, *pssm, sep="\t", file=sys.stdout)

def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('proteinnet', metavar='P', help="Input ProteinNet file")

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())