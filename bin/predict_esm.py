#!/usr/bin/env python3
"""
Make predictions from ESM1b top models
"""
import sys
import argparse
import torch
import proteinnetpy

class LinearModel(torch.nn.Module):
    """Simple linear top layer"""
    def __init__(self, classifier=False):
        super(LinearModel, self).__init__()
        self.linear = torch.nn.Linear(1280, 20)
        self.activation = torch.nn.Sigmoid() if classifier else torch.nn.LogSoftmax(dim=1)

    def forward(self, x):
        out = self.linear(x)
        out = self.activation(out)
        return out
        
def main():
    """
    Run ESM1b on input Fasta and format into a single output file
    """
    args = parse_args()
    use_gpu = torch.cuda.is_available() and not args.no_gpu

    # Load top model
    top_model = torch.load(args.model)
    top_model = top_model.eval()

    # Load ESM1b
    esm_model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm1b_t33_650M_UR50S")
    batch_converter = alphabet.get_batch_converter()
    esm_model = esm_model.eval()

    if use_gpu:
        esm_model = esm_model.cuda()
        print("Transferred ESM1b to GPU", file=sys.stderr)

    # Load ProteinNet
    print(f"Initialising ProtienNet dataset: {args.proteinnet}", file=sys.stderr)
    filter_func = proteinnetpy.data.make_length_filter(max_length=1022)
    proteinnet = proteinnetpy.data.ProteinNetDataset(path=args.proteinnet,
                                                     filter_func=filter_func,
                                                     preload=False)

    # Make predictions
    print("id", "position", "wt",
          *[f"true_{i}" for i in proteinnetpy.record.AMINO_ACIDS],
          *[f"pred_{i}" for i in proteinnetpy.record.AMINO_ACIDS],
          sep="\t", file=sys.stdout)

    for record in proteinnet:
        esm_input = [(record.id, "".join(record.primary))]
        _, _, batch_tokens = batch_converter(esm_input)

        if use_gpu:
            batch_tokens = batch_tokens.to(device="cuda", non_blocking=True)

        with torch.no_grad():
            results = esm_model(batch_tokens, repr_layers=[33], return_contacts=True)
        reps = results["representations"][33].to(device="cpu")

        for p in range(len(record)):
            with torch.no_grad():
                preds = top_model(reps[:,p,:]).numpy()[0,:]
            true = record.evolutionary[:,p]

            print(record.id, p + 1, record.primary[p], *true, *preds, sep="\t", file=sys.stdout)

def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('model', metavar='M', type=str, help="Model file")
    parser.add_argument('proteinnet', metavar='P', type=str, help="Input ProteinNet data")

    parser.add_argument('--no_gpu', '-n', action="store_true",
                        help="Prevent GPU usage for ESM1b even when available")

    return parser.parse_args()

if __name__ == "__main__":
    main()
