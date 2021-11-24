#!/usr/bin/env python3
"""
Generate model predictions from Sequence UNET like models.

Three modes are supported:

ClinVar: Use --proteinnet path/to/protiennet, --clinvar and --tsv path/to/clinvar_tsv where
         the ClinVar TSV is formated like the output of extract_clinvar.R. The ProteinNet file
         must contain the PDB ids specified in the TSV for predictions to be generated for them.

Fasta: Use --fasta path/to/fasta and optionally --tsv to specify a list of variants
       to keep predictions for, where the TSV has columns gene, position, wt, mut and
       genes match the fasta IDs.

ProteinNet: Use --proteinnet path/to/protiennet and optionally --tsv path/to/tsv to specify
            predictions to keep, where the TSV has columns pdb_id, chain, position, wt, mut.
"""
import argparse
import sys
import pandas as pd
from Bio import SeqIO

from sequence_unet.models import load_trained_model
from sequence_unet.predict import predict_proteinnet, predict_sequence
from proteinnetpy.data import ProteinNetDataset, make_id_filter

def main():
    """
    Main function
    """
    args = parse_args()

    tsv = pd.read_csv(args.tsv, sep="\t") if args.tsv else None
    model = load_trained_model(args.model, args.model_dir, download=args.download)

    if args.proteinnet:
        filter_func = make_id_filter(list(tsv.pdb_id), list(tsv.chain)) if tsv is not None else None
        data = ProteinNetDataset(path=args.proteinnet, preload=False, filter_func=filter_func)
        preds = predict_proteinnet(model, data, layers=args.layers,
								   contact=args.contacs, wide=args.wide,
                                   make_pssm=args.pssm)

    elif args.fasta:
        fasta = SeqIO.parse(args.fasta, format="fasta")
        preds = predict_sequence(model, sequences=fasta, layers=args.layers, tsv=tsv,
                                 wide=args.wide, make_pssm=args.pssm)

    else:
        raise ValueError("One of --fasta or --proteinnet must be passed")

    for i, df in enumerate(preds):
        df.to_csv(sys.stdout, sep="\t", index=False, header=i==0)

def arg_parser():
    """Argument parser"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('model', metavar="M", help="Model to predict with")

    inputs = parser.add_argument_group("Input Data")
    inputs.add_argument('--tsv', "-t", help="TSV file containing variants of interest")
    inputs.add_argument('--proteinnet', "-p", help="ProteinNet file")
    inputs.add_argument('--fasta', "-f", help="Fasta file")

    options = parser.add_argument_group("Options")
    options.add_argument('--contacts', "-c", help="Use contact graph input", action="store_true")
    options.add_argument('--layers', "-l", help="Number of layers in bottom UNET model",
                         type=int, default=6)

    options.add_argument('--wide', "-w", help="Output a wide table", action="store_true")
    options.add_argument('--pssm', "-s", help="Convert output frequency predictions to PSSMs",
	                     action="store_true")

    options.add_argument('--model_dir', "-m", help="Directory to locate/download model files to")
    options.add_argument('--download', "-d", action="store_true",
                         help="Download model if not located")

    return parser

def parse_args():
    """Process arguments"""
    parser = arg_parser()
    args = parser.parse_args()

    if args.fasta is None and args.proteinnet is None:
        raise ValueError("Use one of --proteinnet/-p or --fasta/-f")

    if args.fasta is not None and args.proteinnet is not None:
        raise ValueError("Use either --proteinnet/-p or --fasta/-f, not both")

    if args.fasta is not None and args.contacts:
        raise ValueError("Cannot use --contacts with Fasta input")

    if args.clinvar and args.tsv is None:
        raise ValueError("Must pass an annotated ClinVar --tsv when using --clinvar")

    if args.clinvar and args.proteinnet is None:
        raise ValueError("Must use --proteinnet when using --clinvar")

    return args

if __name__ == "__main__":
    main()
