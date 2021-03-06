#!/usr/bin/env python3
"""
Generate model predictions from Sequence UNET like models.

Three modes are supported:

ClinVar: Use --proteinnet path/to/protiennet, --clinvar and --tsv path/to/clinvar_tsv where
         the ClinVar TSV is formated like the output of extract_clinvar.R. The ProteinNet file
         must contain the PDB ids specified in the TSV for predictions to be generated for them.

Fasta: U    se --fasta path/to/fasta and optionally --tsv to specify a list of variants
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

def predict_fasta(model, fasta, layers, tsv=None, wide=False):
    """
    Predict values from a Fasta file
    """
    ind_cols = ["gene", "position", "wt"]
    if not wide:
        ind_cols.append("mut")

    if tsv is not None:
        tsv = tsv[ind_cols]

    for seq in SeqIO.parse(fasta, format="fasta"):
        try:
            one_hot = one_hot_sequence(seq)
        except KeyError as err:
            print(f"Unknown amino acid in {seq.id} - {err}", file=sys.stderr)
            continue

        # Pad to be divisable
        pad_rows = 2 ** (layers - 1) - one_hot.shape[0] % 2 ** (layers - 1) if layers > 0 else 0
        if pad_rows:
            one_hot = np.pad(one_hot, ((0, pad_rows), (0, 0)), mode='constant')
            preds = model(np.array([one_hot])).numpy()[0, :-pad_rows, :]
        else:
            preds = model(np.array([one_hot])).numpy()[0, :, :]

        df = pd.DataFrame(preds, columns=AMINO_ACIDS).reset_index()
        df = df.rename(columns={'index': 'position'})
        df['position'] = df['position'] + 1
        df = df[df.index < len(seq.seq)] # Remove padded records
        df['gene'] = seq.id
        df['wt'] = list(seq.seq)

        if not wide:
            df = df.melt(id_vars=["gene", "position", "wt"], var_name="mut", value_name="pred")

        if tsv is not None:
            df = df.merge(tsv, on=ind_cols, how="inner")

        df = df.sort_values(by=ind_cols)

        yield df[ind_cols + (list(AMINO_ACIDS) if wide else ["pred"])]

def predict_proteinnet(model, data, func, wide=False):
    """
    Predict values from a ProteinNetDataset
    """
    ind_cols = ["pdb_id", "chain", "position", "wt"]
    if not wide:
        ind_cols.append("mut")

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
        df['pdb_id'] = record.pdb_id if record.pdb_id is not None else record.id
        df['chain'] = record.pdb_chain
        df['wt'] = record.primary

        if not wide:
            df = df.melt(id_vars=ind_cols,
                         var_name="mut", value_name="pred")

        df = df.sort_values(by=ind_cols)

        yield df[ind_cols + (list(AMINO_ACIDS) if wide else ["pred"])]

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

    yield pd.DataFrame(data=output)

def main(args):
    """
    Main function
    """
    tsv = pd.read_csv(args.tsv, sep="\t") if args.tsv else None
    model = load_model(args.model, custom_objects=CUSTOM_OBJECTS)

    if args.clinvar:
        preds = predict_clinvar(model, clinvar=tsv, proteinnet=args.proteinnet,
                                layers=args.layers if args.layers else None,
                                contact=args.contacts, pssm=args.pssm)

    elif args.proteinnet:
        filter_func = make_id_filter(list(tsv.pdb_id), list(tsv.chain)) if tsv is not None else None
        data = ProteinNetDataset(path=args.proteinnet, preload=False, filter_func=filter_func)
        func = SequenceUNETMapFunction(num_layers=args.layers if args.layers else None,
                                       contact_graph=args.contacts, pssm=args.pssm, wide=args.wide)
        preds = predict_proteinnet(model, data, func)

    elif args.fasta:
        preds = predict_fasta(model, fasta=args.fasta, layers=args.layers, tsv=tsv, wide=args.wide)

    else:
        raise ValueError("One of --fasta or --proteinnet must be passed")

    for i, df in enumerate(preds):
        df.to_csv(sys.stdout, sep="\t", index=False, header=i==0)

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

    options.add_argument('--wide', "-w", help="Output a wide table", action="store_true")

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
