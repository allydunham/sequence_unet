#!/usr/bin/env python3
"""
Run ESM1b on all sequences in a ProteinNet file. Produces a TSV file with columns for the record ID, position, wt, each PSSM frequency and each representation column.
"""
import sys
import argparse
import torch
import proteinnetpy

def batch_proteinnet_data(data, batch_size=5):
    """
    Generator yielding batches of ProteinNet records
    """
    records = []
    for record in data:
        records.append(record)
        if len(records) == batch_size:
            yield records
            records = []

def main():
    """
    Run ESM1b on input Fasta and format into a single output file
    """
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

    print("id", "position", "wt", *proteinnetpy.record.AMINO_ACIDS,
          *[f"rep{i}" for i in range(1,1281)], sep="\t", file=sys.stdout)

    for i, records in enumerate(batch_proteinnet_data(proteinnet, batch_size=args.batch_size)):
        print(f"Predicting batch {i}", file=sys.stderr)
        esm_input = [(r.id, "".join(r.primary)) for r in records]
        _, _, batch_tokens = batch_converter(esm_input)

        if use_gpu:
            batch_tokens = batch_tokens.to(device="cuda", non_blocking=True)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=True)
        token_representations = results["representations"][33].to(device="cpu")

        for record, rep in zip(records, token_representations):
            for p in range(len(record.primary)):
                print(record.id, p + 1, record.primary[p], *record.evolutionary[:,p],
                      *rep.numpy()[p,:], sep="\t", file=sys.stdout)

def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('proteinnet', metavar='P', type=str, help="Input ProteinNet data")

    parser.add_argument('--no_gpu', '-n', action="store_true",
                        help="Prevent GPU usage even when available")

    parser.add_argument('--batch_size', '-b', default=5, type=int, help="Batch size")

    return parser.parse_args()

if __name__ == "__main__":
    main()
