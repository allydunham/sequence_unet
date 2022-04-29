#!/usr/bin/env python3
"""
Make softmaxed logit predictions from ESM models
"""
import sys
import argparse
from typing import Type
import torch
import gc
import proteinnetpy

from Bio import SeqIO

MODEL_STRS = {
    f"esm1v_{n}": f"esm1v_t33_650M_UR90S_{n}" for n in range(1,6)
}
MODEL_STRS["esm1b"] = "esm1b_t33_650M_UR50S"

def load_esm_model(model, use_gpu=True):
    """
    Load an ESM model
    """
    esm_model, alphabet = torch.hub.load("facebookresearch/esm:main", model)
    batch_converter = alphabet.get_batch_converter()
    esm_model = esm_model.eval()

    if torch.cuda.is_available() and use_gpu:
        esm_model = esm_model.cuda()
        print("Transferred ESM model to GPU", file=sys.stderr)

    return esm_model, alphabet, batch_converter

def main():
    """
    Run ESM1b on input Fasta and format into a single output file
    """
    args = parse_args()
    use_gpu = torch.cuda.is_available() and not args.no_gpu

    print("model", "id", "position", "wt", *[f"pred_{i}" for i in proteinnetpy.record.AMINO_ACIDS],
          sep="\t", file=sys.stdout)
    for model_key in args.models:
        # Load ESM1b
        model, alphabet, batch_converter = load_esm_model(MODEL_STRS[model_key], use_gpu=use_gpu)
        mut_idxs = [alphabet.get_idx(i) for i in proteinnetpy.record.AMINO_ACIDS]

        # Load Data
        if args.input_type == "proteinnet":
            print(f"Initialising ProtienNet dataset: {args.input}", file=sys.stderr)
            filter_func = proteinnetpy.data.make_length_filter(max_length=1022)
            data = proteinnetpy.data.ProteinNetDataset(path=args.input,
                                                       filter_func=filter_func,
                                                       preload=False)
        elif args.input_type == "fasta":
            data = SeqIO.parse(args.input, "fasta")
        else:
            # Should be caught by argparse before here
            raise ValueError("Unrecognised input file type")


        # Make predictions
        for record in data:
            if len(record) > 1022:
                # Maximum size expected by ESM-1v
                continue

            if args.input_type == "proteinnet":
                rec_id = record.id
                seq = "".join(record.primary)
            elif args.input_type == "fasta":
                rec_id = record.id
                seq = str(record.seq)
            else:
                # Should be caught by argparse before here
                raise ValueError("Unrecognised input file type")

            esm_input = [(rec_id, seq)]
            _, _, batch_tokens = batch_converter(esm_input)

            if use_gpu:
                batch_tokens = batch_tokens.to(device="cuda", non_blocking=True)

            with torch.no_grad():
                results = torch.log_softmax(model(batch_tokens)["logits"], dim=-1)[0]
            results = results.cpu().numpy()

            for p in range(len(record)):
                wt_idx = alphabet.get_idx(seq[p])
                if args.wt_norm:
                    scores = results[p + 1, mut_idxs] - results[p + 1, wt_idx]
                else:
                    scores = results[p + 1, mut_idxs]

                print(model_key, rec_id, p + 1, seq[p], *scores,
                      sep="\t", file=sys.stdout)

        # Wipe model from GPU to clear space for next model
        model = model.cpu()
        del model
        gc.collect()
        torch.cuda.empty_cache()


def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input', metavar='I', type=str, help="Input file")

    parser.add_argument('--input_type', '-i', default="fasta", choices=["proteinnet", "fasta"],
                        help="Input file type")

    parser.add_argument('--models', '-m', nargs="+", default=["esm1b"], type=str,
                        choices=list(MODEL_STRS.keys()), help="ESM models to predict with")

    parser.add_argument('--wt_norm', '-w', action="store_true",
                        help="Normalise logits to WT amino acid")

    parser.add_argument('--no_gpu', '-n', action="store_true",
                        help="Prevent GPU usage for ESM1b even when available")

    return parser.parse_args()

if __name__ == "__main__":
    main()
