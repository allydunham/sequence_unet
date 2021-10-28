#!/usr/bin/env Python3
"""
ProtienNetMap LabeledFunctions for working with the Sequence UNET model
"""
from dataclasses import dataclass

import numpy as np
import pandas as pd
from Bio import SeqIO

from sequence_unet.graph_cnn import contact_graph
from proteinnetpy.data import LabeledFunction

AMINO_ACIDS = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                        'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
AA_HASH = {aa: index for index, aa in enumerate(AMINO_ACIDS)}

class SequenceUNETMapFunction(LabeledFunction):
    """
    ProteinNetMapFunction returning the required data for the sequence UNET model
    """
    def __init__(self, num_layers=4, threshold=None, contact_graph=False,
                 weights=False, pssm=False):
        self.num_layers = num_layers
        self.threshold = threshold
        self.contact_graph = contact_graph
        self.weights = weights
        self.pssm = pssm

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
        data = record.evolutionary.T if self.pssm else record.get_one_hot_sequence().T
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

class ClinVarMapFunction(LabeledFunction):
    """
    ProteinNetMapFunction returning the required data for the Mutation PSSM Top model

    clinvar: Annotated clinvar variants pandas data frame
    num_layers: Number of layers of UNET model
    contact_graph: Include contact graph in output
    pssm: Return the true PSSM instead of the sequence
    """
    def __init__(self, clinvar, num_layers=None, contact_graph=True, pssm=False):
        self.clinvar = clinvar
        self.clinvar['pdb_id'] = self.clinvar['pdb_id'].str.upper()

        self.num_layers = num_layers
        self.contact_graph = contact_graph
        self.pssm = pssm

        out_shape = ([None, 20], [None, None]) if self.contact_graph else ([None, 20],)
        out_type = ('float32', 'float32') if self.contact_graph else ('float32',)

        # (Seq, struct graph), labels
        self.output_shapes = (out_shape, [None, 20])
        self.output_types = (out_type, 'int32')

    def __call__(self, record):
        out_main = record.evolutionary.T if self.pssm else record.get_one_hot_sequence().T
        contacts = contact_graph(record, contact_distance=800) if self.contact_graph else None

        labels = self._calc_clinvar(record, shape=out_main.shape)

        if self.num_layers is not None:
            # Need to make sure output can be halved a sufficient number of time
            pad_rows = 2 ** (self.num_layers - 1) - out_main.shape[0] % 2 ** (self.num_layers - 1)
            if pad_rows:
                out_main = np.pad(out_main, ((0, pad_rows), (0, 0)), mode='constant')
                labels = np.pad(labels, ((0, pad_rows), (0, 0)), mode='constant')
                if self.contact_graph:
                    contacts = np.pad(contacts, ((0, pad_rows), (0, pad_rows)), mode='constant')

        if self.contact_graph:
            return (out_main, contacts), labels
        else:
            return (out_main,), labels

    def _calc_clinvar(self, record, shape):
        """
        Generate clinvar label matrix for a given record
        """
        try:
            muts = self.clinvar[((self.clinvar.pdb_id == record.pdb_id) &
                                (self.clinvar.chain == record.pdb_chain))]
            # Offset pos by first unmasked index, since that is what first_pdb_pos corresponds to
            pos = (np.nonzero(record.mask)[0][0] + muts.pdb_pos - muts.first_pdb_pos).values
            aa = np.vectorize(AA_HASH.__getitem__)(muts.mut)

            benign = muts.clnsig.isin(["Benign", "Likely_benign"])
            path = muts.clnsig.isin(["Likely_pathogenic", "Pathogenic", "Pathogenic/Likely_pathogenic"])
            labels = np.full(shape, fill_value=0, dtype=int) # 0 means masked
            labels[pos[benign], aa[benign]] = 1 # These are +1 offset compared to model out index
            labels[pos[path], aa[path]] = 2

        except (ValueError, IndexError) as err:
            # Raise errors as Value error, so they are filtered by map (bit of a hack)
            raise ValueError(str(err))

        return labels

def one_hot_sequence(seq):
    """
    Convert a Biopython AA sequnece to one hot representation
    """
    indeces = np.array([AA_HASH[aa] for aa in seq])
    one_hot = np.zeros((len(indeces), 20), dtype=np.int)
    one_hot[np.arange(len(indeces)), indeces] = 1
    return one_hot
