#!/usr/bin/env python3
"""
Time predictions from Sequence UNET
"""
import argparse
import sys
import numpy as np
import pandas as pd
from tensorflow.keras.models import load_model

from proteinnetpy.data import ProteinNetDataset
from proteinnetpy.record import AMINO_ACIDS
from metrics import CUSTOM_OBJECTS

PADDING = [
    ("0", 0), ("32", 32), ("64", 32), ("128", 64),
    ("256", 128), ("512", 256), ("1024", 512)
]

def main(args):
    """
    Main function
    """
    model = load_model(args.model, custom_objects=CUSTOM_OBJECTS)
    data = ProteinNetDataset(path=args.proteinnet, preload=False, filter_func=None)

    for i, record in enumerate(data):
        print(f"Processing {record.id}", file=sys.stderr)

        one_hot = record.get_one_hot_sequence().T

        pad_rows = 2 ** 5 - one_hot.shape[0] % 2 ** 5
        if pad_rows:
            one_hot = np.pad(one_hot, ((0, pad_rows), (0, 0)), mode='constant')

        preds = {}
        for name, pad in PADDING:
            one_hot = np.pad(one_hot, ((0, pad), (0, 0)), mode='constant')
            preds[name] = model(np.array([one_hot])).numpy()[0, :len(record), :]

        df = pd.DataFrame(record.evolutionary.T, columns=AMINO_ACIDS).reset_index()
        df = df.rename(columns={'index': 'position'})
        df['position'] = df['position'] + 1
        df['id'] = record.id
        df['wt'] = record.primary
        df = df.melt(id_vars=["id", "position", "wt"], var_name="mut", value_name="freq")

        for k, v in preds.items():
            df[f"pad_{k}"] = v.flatten(order="F")

        df = df.sort_values(by=["id", "position", "wt", "mut"])

        df.to_csv(sys.stdout, sep="\t", index=False, header=i==0)

def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('model', metavar="M", help="Model to predict with")
    parser.add_argument('proteinnet', metavar="F", help="ProteinNet file")

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())
