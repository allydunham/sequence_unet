#!/usr/bin/env python3
"""
Time predictions from Sequence UNET
"""
import argparse
import sys
import time
import numpy as np
from Bio import SeqIO
from tensorflow.keras.models import load_model

from pn_maps import one_hot_sequence
from metrics import CUSTOM_OBJECTS

def main(args):
    """
    Main function
    """
    model = load_model(args.model, custom_objects=CUSTOM_OBJECTS)

    l = args.layers

    print("id", "length", "time", sep="\t")
    for seq in SeqIO.parse(args.fasta, format="fasta"):
        try:
            one_hot = one_hot_sequence(seq)
        except KeyError as err:
            print(f"Unknown amino acid in {seq.id} - {err}", file=sys.stderr)
            continue

        t0 = time.perf_counter()
        # Pad to be divisable
        pad_rows = 2 ** (l - 1) - one_hot.shape[0] % 2 ** (l - 1) if l > 0 else 0
        if pad_rows:
            one_hot = np.pad(one_hot, ((0, pad_rows), (0, 0)), mode='constant')
            preds = model(np.array([one_hot])).numpy()[0, :-pad_rows, :]
        else:
            preds = model(np.array([one_hot])).numpy()[0, :, :]
        t1 = time.perf_counter()

        print(seq.id, len(seq), t1 - t0, sep="\t")


def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('model', metavar="M", help="Model to predict with")
    parser.add_argument('fasta', metavar="F", help="Fasta file")

    parser.add_argument('--layers', "-l", help="Number of layers in bottom UNET model",
                        type=int, default=6)

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())
