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
            This may require a large amount of memory if the ProteinNet file is large.
"""
import argparse
import sys
import numpy as np
import pandas as pd
from Bio import SeqIO
from tensorflow.keras.models import load_model

from pn_maps import SequenceUNETMapFunction, ClinVarMapFunction, one_hot_sequence, AMINO_ACIDS
from proteinnetpy.data import ProteinNetDataset, make_id_filter
from metrics import CUSTOM_OBJECTS

def predict_fasta(model, fasta, layers, tsv=None):
    """
    Predict values from a Fasta file
    """
    tsv = tsv[["gene", "position", "wt", "mut"]]
    output = []
    for seq in SeqIO.parse(fasta, format="fasta"):
        one_hot = one_hot_sequence(seq)

        # Pad to be divisable
        pad_rows = 2 ** (layers - 1) - one_hot.shape[0] % 2 ** (layers - 1)
        if pad_rows:
            one_hot = np.pad(one_hot, ((0, pad_rows), (0, 0)), mode='constant')

        preds = model(np.array([one_hot])).numpy()[0, :-pad_rows, :]

        df = pd.DataFrame(preds, columns=AMINO_ACIDS).reset_index()
        df = df.rename(columns={'index': 'position'})
        df['position'] = df['position'] + 1
        df = df[df.index < len(seq.seq)] # Remove padded records
        df['gene'] = seq.id
        df['wt'] = seq.seq
        df = df.melt(id_vars=["gene", "position", "wt"], var_name="mut", value_name="pred")

        output.append(df[["gene", "position", "wt", "mut", "pred"]])

    output = pd.concat(output, axis=0)
    if tsv is not None:
        output = output.merge(tsv, on=["gene", "position", "wt", "mut"], how="right")
    return output

# TODO maybe have to stream this if want to predict on larger datasets
def predict_proteinnet(model, data, func, tsv=None):
    """
    Predict values from a ProteinNetDataset
    """
    tsv = tsv[["pdb_id", "chain", "position", "wt", "mut"]]
    output = []
    for record in data:
        # Format model input
        try:
            mod_input = func(record)
        except (ValueError, IndexError) as err:
            print(err, file=sys.stderr)
            continue

        # Predict
        x = [np.array([i]) for i in mod_input[0]]
        preds = model(x).numpy()[0, :, :]

        df = pd.DataFrame(preds, columns=AMINO_ACIDS).reset_index()
        df = df.rename(columns={'index': 'position'})
        df['position'] = df['position'] + 1
        df = df[df.index < len(record)] # Remove padded records
        df['pdb_id'] = record.pdb_id
        df['chain'] = record.pdb_chain
        df['wt'] = record.primary
        df = df.melt(id_vars=["pdb_id", "chain", "position", "wt"],
                     var_name="mut", value_name="pred")

        output.append(df[["pdb_id", "chain", "position", "wt", "mut", "pred"]])

    output = pd.concat(output, axis=0)
    if tsv is not None:
        output = output.merge(tsv, on=["pdb_id", "chain", "position", "wt", "mut"], how="right")
    return output

def predict_clinvar(model, clinvar, proteinnet, layers, contact, pssm):
    """
    Predict values from a ProteinNetDataset
    """
    filter_func = make_id_filter(list(clinvar.pdb_id), list(clinvar.chain))
    data = ProteinNetDataset(path=proteinnet, preload=False, filter_func=filter_func)
    func = ClinVarMapFunction(clinvar=clinvar, num_layers=layers, contact_graph=contact, pssm=pssm)

    output = {"pdb_id": [], "chain": [], "pdb_pos": [], "wt": [], "mut": [], "pred": []}
    for record in data:
        cln = clinvar[((clinvar.pdb_id == record.pdb_id) &
                       (clinvar.chain == record.pdb_chain))]

        offset = cln.first_pdb_pos.values[0] - np.nonzero(record.mask)[0][0]

        # Format model input
        try:
            mod_input = func(record)
        except (ValueError, IndexError) as err:
            print(err, file=sys.stderr)
            continue

        x = [np.array([i]) for i in mod_input[0]]

        # Predict, normalise and select del preds
        preds = model(x).numpy()[0, :, :]

        # Select relavent predictions
        positions = np.where(mod_input[1] > 0)
        preds = preds[positions[0], positions[1]]
        wt = record.primary[positions[0]]
        mut = AMINO_ACIDS[positions[1]]
        pos = positions[0] + offset

        # Add to output
        output["pdb_id"].extend([record.pdb_id] * len(preds))
        output["chain"].extend([record.pdb_chain] * len(preds))
        output["pdb_pos"].extend(pos)
        output["wt"].extend(wt)
        output["mut"].extend(mut)
        output["pred"].extend(preds)

    return pd.DataFrame(data=output)

def main(args):
    """
    Main function
    """
    tsv = pd.read_csv(args.tsv, sep="\t") if args.tsv else None
    model = load_model(args.model, custom_objects=CUSTOM_OBJECTS)

    if args.clinvar:
        preds = predict_proteinnet(model, clinvar=tsv, proteinnet=args.proteinnet,
                                   layers=args.layers, contact=args.contact,
                                   pssm=args.pssm)
        preds.to_csv(sys.stdout, sep="\t", index=False)

    elif args.proteinnet:
        filter_func = make_id_filter(list(tsv.pdb_id), list(tsv.chain)) if tsv is not None else None
        data = ProteinNetDataset(path=args.proteinnet, preload=False, filter_func=filter_func)
        func = SequenceUNETMapFunction(num_layers=args.layers if args.layers else None,
                                       contact_graph=args.contact,
                                       pssm=args.pssm)
        preds = predict_proteinnet(model, data, func)

    elif args.fasta:
        preds = predict_fasta(model, args.fasta, args.layers, args.tsv)
        preds.to_csv(sys.stdout, sep="\t", index=False)

    else:
        raise ValueError("One of --fasta or --proteinnet must be passed")

def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('model', metavar="M", help="Model to predict with")

    inputs = parser.add_argument_group("Input Data")
    inputs.add_argument('--tsv', "-t", help="TSV file containing variants of interest")
    inputs.add_argument('--clinvar', "-v", action="store_true",
                        help="TSV file contains ClinVar variants")
    inputs.add_argument('--proteinnet', "-p", help="ProteinNet file")
    inputs.add_argument('--fasta', "-f", help="Fasta file")

    options = parser.add_argument_group("Options")
    options.add_argument('--contacts', "-c", help="Use contact graph input", action="store_true")
    options.add_argument('--pssm', "-s", help="Use PSSM input", action="store_true")
    options.add_argument('--layers', "-l", help="Number of layers in bottom UNET model",
                        type=int, default=6)

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
    main(parse_args())
