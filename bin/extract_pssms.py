#!/usr/bin/env python3
"""
Extract PSSM profiles from a ProteinNet file into the same form as SPBuild/
"""
import sys
import argparse
import numpy as np
import pandas as pd
from tensorflow.keras.models import load_model
from proteinnetpy.data import ProteinNetDataset, ProteinNetMap
from models.frequency_unet import UnetPNMapFunction

AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')

def format_profile(record):
    """
    Format a predicted profile into a table
    """
    df = pd.DataFrame(record.evolutionary.T, columns=AMINO_ACIDS).reset_index()
    df = df.rename(columns={'index': 'position'})
    df['position'] = df['position'] + 1
    df['protein'] = record.id
    df['wt'] = record.primary
    return df[['protein', 'position', 'wt'] + AMINO_ACIDS]

def main(args):
    """
    Run the Frequency UNET model on all records of a ProteinNet file
    """
    data = ProteinNetDataset(args.proteinnet, preload=False)
    df = pd.concat([format_profile(r) for r in data], axis=0)
    df.to_csv(sys.stdout, sep='\t', index=False)

def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('proteinnet', metavar='P', help="Input ProteinNet file")

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())
