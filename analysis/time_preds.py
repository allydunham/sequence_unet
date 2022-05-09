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

def time_single(model, fasta, layers):
    """
    Time prediction against individual proteins.
    """
    print("id", "length", "time", sep="\t")
    for seq in fasta:
        try:
            one_hot = one_hot_sequence(seq)
        except KeyError as err:
            print(f"Unknown amino acid in {seq.id} - {err}", file=sys.stderr)
            continue

        t0 = time.perf_counter()
        # Pad to be divisable
        pad_rows = 2 ** (layers - 1) - one_hot.shape[0] % 2 ** (layers - 1) if layers > 0 else 0
        if pad_rows:
            one_hot = np.pad(one_hot, ((0, pad_rows), (0, 0)), mode='constant')
            preds = model(np.array([one_hot])).numpy()[0, :-pad_rows, :]
        else:
            preds = model(np.array([one_hot])).numpy()[0, :, :]
        t1 = time.perf_counter()

        print(seq.id, len(seq), t1 - t0, sep="\t")

def batch_fasta(fasta, batch_size):
    """
    Batch fasta sequences
    """
    batch = []
    for seq in fasta:
        try:
            one_hot = one_hot_sequence(seq)
        except KeyError as err:
            print(f"Unknown amino acid in {seq.id} - {err}", file=sys.stderr)
            continue

        batch.append((seq, one_hot))
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch # Last batch may be smaller

def time_bactched(model, fasta, layers, batch_size):
    """
    Predict prediction time in batches of proteins.
    """
    print("batch", "id", "length", "time", sep="\t")
    for batch_num, batch in enumerate(batch_fasta(fasta, batch_size)):
        t0 = time.perf_counter()

        longest = max(i.shape[0] for _, i in batch)
        halving_padding = 2 ** (layers - 1) - longest % 2 ** (layers - 1) if layers > 0 else 0
        length = longest + halving_padding

        model_input = np.stack(
            [np.pad(i, ((0, length - i.shape[0]), (0, 0)), mode='constant') for _, i in batch]
        )

        preds = model(model_input).numpy()

        t1 = time.perf_counter()
        per_prot_time = (t1 - t0) / len(batch)

        for seq, _ in batch:
            print(batch_num, seq.id, len(seq), per_prot_time, sep="\t")

def main(args):
    """
    Main function
    """
    model = load_model(args.model, custom_objects=CUSTOM_OBJECTS)
    fasta = SeqIO.parse(args.fasta, format="fasta")

    if args.batch_size == 1:
        time_single(model, fasta, args.layers)
    else:
        time_bactched(model, fasta, args.layers, args.batch_size)

def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('model', metavar="M", help="Model to predict with")
    parser.add_argument('fasta', metavar="F", help="Fasta file")

    parser.add_argument('--layers', "-l", help="Number of layers in bottom UNET model",
                        type=int, default=6)

    parser.add_argument('--batch_size', "-b", help="Batch size",
                        type=int, default=1)

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())
