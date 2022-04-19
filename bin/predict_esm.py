#!/usr/bin/env python3
"""
Run ESM1b on all sequences in a ProteinNet file. Produces a TSV file with columns for the record ID, position, wt, each PSSM frequency and each representation column.
"""
import sys
import argparse
import torch
import proteinnetpy

def main():
    """
    Run ESM1b on input Fasta and format into a single output file
    """
    args = parse_args()

    model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm1b_t33_650M_UR50S")
    batch_converter = alphabet.get_batch_converter()
    model.eval()

    if torch.cuda.is_available() and not args.nogpu:
        model = model.cuda()
        print("Transferred model to GPU", file=sys.stderr)

    filter_func = proteinnetpy.data.make_length_filter(max_length=1022)
    proteinnet = proteinnetpy.data.ProteinNetDataset(path=args.proteinnet, filter_func=filter_func)

    print("id", "position", "wt", *proteinnetpy.record.AMINO_ACIDS,
          *[f"rep{i}" for i in range(1279)], sep="\t", file=sys.stdout)
    for record in proteinnet:
        print(f"Predicting {record.id}", file=sys.stderr)
        data = [(record.id, "".join(record.primary))]
        _, _, batch_tokens = batch_converter(data)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)
        token_representations = results["representations"][33]







def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('proteinnet', metavar='P', help="Input ProteinNet data")

    parser.add_argument('--no_gpu', '-n', action="store_true",
                        help="Prevent GPU usage even when available")

    return parser.parse_args()

if __name__ == "__main__":
    main()
