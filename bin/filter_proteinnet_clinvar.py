#!/usr/bin/env python3
"""
Filter a ProteinNet file to only include records contained in a ClinVar dataset.
A ClinVar TSV of the format produced by extract_clinvar.R is expected.
"""
import argparse
import sys
import pandas as pd

from proteinnetpy.data import ProteinNetDataset
from models.pssm_top_model import make_id_filter

def main(args):
    """
    Run the Frequency UNET model on all records of a ProteinNet file
    """
    clinvar = pd.read_csv(args.clinvar, sep="\t")
    filter_func = make_id_filter(list(clinvar.pdb_id), list(clinvar.chain))
    data = ProteinNetDataset(path=args.proteinnet, preload=False, filter_func=filter_func)

    for rec in data:
        print(rec, file=sys.stdout)

def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('proteinnet', metavar='P', help="Input ProteinNet file")
    parser.add_argument('clinvar', metavar='C', help="Input ClinVar TSV file")

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())