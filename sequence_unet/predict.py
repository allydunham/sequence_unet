"""
Predict values using SequenceUNET models
"""
import logging
import numpy as np
import pandas as pd
from Bio import SeqIO

from sequence_unet.graph_cnn import contact_graph
from proteinnetpy.data import LabeledFunction

__all__ = ["SequenceUNETMapFunction", "predict_fasta", "predict_proteinnet"]

AMINO_ACIDS = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                        'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
AA_HASH = {aa: index for index, aa in enumerate(AMINO_ACIDS)}

def one_hot_sequence(seq):
    """
    Convert a Biopython AA sequnece to one hot representation
    """
    indeces = np.array([AA_HASH[aa] for aa in seq])
    one_hot = np.zeros((len(indeces), 20), dtype=np.int)
    one_hot[np.arange(len(indeces)), indeces] = 1
    return one_hot

class SequenceUNETMapFunction(LabeledFunction):
    """
    ProteinNetMapFunction returning the required data for the sequence UNET model
    """
    def __init__(self, num_layers=4, threshold=None, contact_graph=False,
                 weights=False, freq_input=False):
        self.num_layers = num_layers
        self.threshold = threshold
        self.contact_graph = contact_graph
        self.weights = weights
        self.freq_input = freq_input

        output_shapes = ([None, 20], [None, None]) if self.contact_graph else ([None, 20],)
        output_types = ('float32', 'float32') if self.contact_graph else ('float32',)

        output_shapes = [output_shapes, [None, 20]]
        output_types = [output_types, 'int32' if self.threshold is not None else 'float32']

        if self.weights:
            output_shapes.append([None])
            output_types.append('float32')

        self.output_shapes = tuple(output_shapes)
        self.output_types = tuple(output_types)

    def __call__(self, record):
        data = record.evolutionary.T if self.freq_input else record.get_one_hot_sequence().T
        labels = record.evolutionary.T
        weights = np.ones(data.shape[0])

        if self.threshold is not None:
            labels = labels < self.threshold
            labels = labels.astype(int) + 1

        # Need to make sure output can be halved a sufficient number of time
        pad_rows = 2 ** (self.num_layers - 1) - data.shape[0] % 2 ** (self.num_layers - 1)
        if pad_rows:
            data = np.pad(data, ((0, pad_rows), (0, 0)), mode='constant')
            labels = np.pad(labels, ((0, pad_rows), (0, 0)), mode='constant')
            weights = np.pad(weights, (0, pad_rows), mode='constant')

        if self.contact_graph:
            contacts = contact_graph(record)
            if pad_rows:
                contacts = np.pad(contacts, ((0, pad_rows), (0, pad_rows)), mode='constant')

        out = (data, contacts) if self.contact_graph else (data,)

        if self.weights:
            return out, labels, weights
        else:
            return out, labels

def predict_fasta(model, fasta, layers, variants=None, wide=False):
    """
    Predict values from a Fasta file
    """
    ind_cols = ["gene", "position", "wt"]
    if not wide:
        ind_cols.append("mut")

    if variants is not None:
        variants = variants[ind_cols]

    for seq in SeqIO.parse(fasta, format="fasta"):
        try:
            one_hot = one_hot_sequence(seq)
        except KeyError as err:
            logging.log(logging.WARN, "Skipping %s: unknown amino acid (%s)", seq.id, err)
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
        df['wt'] = seq.seq

        if not wide:
            df = df.melt(id_vars=["gene", "position", "wt"], var_name="mut", value_name="pred")

        if variants is not None:
            df = df.merge(variants, on=ind_cols, how="inner")

        df = df.sort_values(by=ind_cols)

        yield df[ind_cols + (list(AMINO_ACIDS) if wide else ["pred"])]

def predict_proteinnet(model, data, layers=6, contacts=False, wide=False):
    """
    Generator yielding predictions from
    """
    ind_cols = ["pdb_id", "chain", "position", "wt"]
    if not wide:
        ind_cols.append("mut")

    func = SequenceUNETMapFunction(num_layers=layers if layers else 1, contact_graph=contacts)

    for record in data:
        try:
            mod_input = func(record)
        except (ValueError, IndexError) as err:
            logging.log(logging.WARN, "Skipping %s: %s", record.id, err)
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
