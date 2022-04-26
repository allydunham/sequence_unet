#!/usr/bin/env python3
"""
Time execution of ESM1b model using ProteinNet data
"""
import sys
import argparse
import time
import torch
import proteinnetpy

def main():
    """Main"""
    args = parse_args()

    model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm1b_t33_650M_UR50S")
    batch_converter = alphabet.get_batch_converter()
    model = model.eval()

    use_gpu = torch.cuda.is_available() and not args.no_gpu

    if use_gpu:
        model = model.cuda()
        print("Transferred model to GPU", file=sys.stderr)

    filter_func = proteinnetpy.data.make_length_filter(max_length=1022)
    proteinnet = proteinnetpy.data.ProteinNetDataset(path=args.proteinnet, filter_func=filter_func,
                                                     preload=False)

    print("id", "length", "time", sep="\t", file=sys.stdout)

    for record in proteinnet:
        esm_input = [(record.id, "".join(record.primary))]

        # Time from input generation to returning data to CPU
        t0 = time.perf_counter()

        _, _, batch_tokens = batch_converter(esm_input)

        if use_gpu:
            batch_tokens = batch_tokens.to(device="cuda", non_blocking=True)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)
        token_representations = results["representations"][33].to(device="cpu")

        t1 = time.perf_counter()

        print(record.id, len(record), t1 - t0, sep="\t")

def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('proteinnet', metavar='P', type=str, help="Input ProteinNet data")

    parser.add_argument('--no_gpu', '-n', action="store_true",
                        help="Prevent GPU usage for ESM1b even when available")

    return parser.parse_args()

if __name__ == "__main__":
    main()