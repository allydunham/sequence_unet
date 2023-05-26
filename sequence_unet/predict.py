"""
Functions for predicting with Sequence UNET models.
"""
import logging
import numpy as np
import pandas as pd

from sequence_unet.graph_cnn import contact_graph
from proteinnetpy.data import LabeledFunction

__all__ = ["SequenceUNETMapFunction", "predict_sequence", "predict_proteinnet",
           "one_hot_sequence", "padding_rows"]

AMINO_ACIDS = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
                        'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
"""
Numpy array of Amino Acid single letter codes.
"""

AA_HASH = {aa: index for index, aa in enumerate(AMINO_ACIDS)}
"""
Dictionary mapping Amino Acids to their alphabetic index.
"""

# from ExPASy Data Portal (https://web.expasy.org/docs/relnotes/relstat.html, 17/08/2020)
# Values are for AAs alphabetically
AA_FREQS = np.array([8.25, 1.38, 5.46, 6.72, 3.86, 7.08, 2.27, 5.91, 5.80, 9.65,
                     2.41, 4.05, 4.73, 3.93, 5.53, 6.63, 5.35, 6.86, 1.09, 2.92]) / 100
"""
Frequency of each Amino Acid in SwissProt. These were collected from ExPASY Data Portal (https://web.expasy.org/docs/relnotes/relstat.html) on 17/08/2020. Amino Acids are ordered alphabetically.
"""

def one_hot_sequence(seq):
    """
    Convert a Biopython AA sequnece to one hot representation.

    Parameters
    ----------
    seq : BioPython Sequence
        Sequence to convert to one-hot matrix
    Returns
    -------
    `numpy.array`
        One-hot encoded sequence array
    """
    indeces = np.array([AA_HASH[aa] for aa in seq])
    one_hot = np.zeros((len(indeces), 20), dtype=int)
    one_hot[np.arange(len(indeces)), indeces] = 1
    return one_hot

def padding_rows(length, layers=6):
    """
    Calculate the number of padding rows for an input sequence matrix.

    Calculate the number of padding rows required by Sequence UNET for an input sequence matrix. This number of 0 rows should be added to the end of the array so it can be evenly halved a sufficient number of times.

    Parameters
    ----------
    length : int
        Sequence length
    layers : int
        Number of layers in the Sequence UNET model. All pre-trained models have 6 layers.

    Returns
    -------
    int
        Number of padding rows required.
    """
    return 2 ** (layers - 1) - length % 2 ** (layers - 1) if layers > 0 else 0

def freqs_to_pssm(mat):
    """
    Convert a frequency matrix to a PSSM.

    Convert a matrix of raw multiple sequence alignment frequency predictions in a position specific scoring matrix (PSSM). This has the log2 ratio of frequency for each amino acid at each position against the baseline occurance of that amino acid:

    log2 (Freq) / (Baseline AA Freq)

    The baseline amino acid frequencies can be found in `AA_FREQS` and were collected from ExPASY Data Portal (https://web.expasy.org/docs/relnotes/relstat.html) on 17/08/2020.

    Parameters
    ----------
    mat : Numpy array
        N x 20 array of MSA frequencies (floats in the range 0 to 1). This can optionally have an additional batch dimension to transform multiple PSSMs at once.

    Returns
    -------
    Numpy array
        (B x) N x 20 PSSM matrix.
    """
    return np.log2((mat + 0.00001) / AA_FREQS[:,None]).astype(int)

class SequenceUNETMapFunction(LabeledFunction):
    """
    ProteinNetPy LabeledFunction returning the required data for the sequence UNET model.

    LabeledFunction taking a ProteinNetPy record and returning the data required to run Sequence UNET on that record. This can be used with a ProteinNetPy map to generate input data for training or predictions. The function outputs a tuple containing model input, N x 20 labels matrix and optionally N long position weights vector . Model input is itself a tuple contaning the N x 20 one-hot encoded sequence and optionally the protein's amino acid contact graph.

    Attributes
    ----------
    num_layers    : int
        Expected number of layers in the target Sequence UNET model.
    threshold     : float or None
        If set, the function returns categorical output, with variants below this threshold classed as deleterious.
    contact_graph : bool
        The function outputs a contact graph in addition to the one-hot encoded sequence in the model input tuple (the first element of the returned tuple).
    weights       : bool
        The function outputs per position sample weights as the third element of the output tuple, with true positions weighted 1 and padded positions weighted 0.
    freq_input    : bool
        The function outputs true PSSM frequencies instead of one-hot encoded sequences.
    output_shapes : tuple
        Potentially nested tuple describing the shape of the arrays output by the function.
    output_types  : tuple
        Potentially nested tuple describing the type of the arrays output by the function.
    """
    def __init__(self, num_layers=4, threshold=None, contact_graph=False,
                 weights=False, freq_input=False):
        """
        Initialise the SequenceUNETMapFunction.

        Parameters
        ----------
        num_layers     : int
            Number of layers in the Sequence UNET model.
        threshold      : float or None
            Return categorical output, with variants below this threshold classed as deleterious.
        contact_graph  : bool
            Return a contact graph in addition to the one-hot encoded sequence in the model input tuple (the first element of the returned tuple).
        weights        : bool
            Return per position sample weights as the third element of the output tuple, with true positions weighted 1 and padded positions weighted 0.
        freq_input     : bool
            Use true PSSM frequencies instead of one-hot encoded sequence as the main model input. This option is not used by any included pre-trained models.
        """
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
        """
        Generate Sequence UNET input data from a ProteinNet record.

        Generate Sequence UNET input data from a ProteinNet record. The function outputs a tuple containing model input, N x 20 labels matrix and optionally N long position weights vector . Model input is itself a tuple contaning the N x 20 one-hot encoded sequence and optionally the protein's amino acid contact graph.

        Parameters
        ----------
        record : ProteinNetPy ProteinNetRecord
            ProteinNet record to process.

        Returns
        -------
        tuple
            Model input, labels and optionally weights to be passed to Sequence UNET. See above for details.
        """
        data = record.evolutionary.T if self.freq_input else record.get_one_hot_sequence().T
        labels = record.evolutionary.T
        weights = np.ones(data.shape[0])

        if self.threshold is not None:
            labels = labels < self.threshold
            labels = labels.astype(int) + 1

        # Need to make sure output can be halved a sufficient number of time
        pad_rows = padding_rows(data.shape[0], self.num_layers)
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

def predict_sequence(model, sequences, layers=6, variants=None, wide=False, make_pssm=False):
    """
    Predict values from a iterable of sequence strings

    Parameters
    ----------
    model     : keras.Model
        Sequence UNET model to predict with (or another with the same input/output signature).
    sequences      : Iterable of str or BioPython sequences
        Sequences to predict from.
    layers    : int
        Number of layers in the Sequence UNET model. All pre-trained models have 6 layers.
    variants  : Pandas DataFrame or None
        Filter variants to those contained in this data frame. It must have the same identifying columns as the output table: gene, position & wt if `wide=true` and gene, position, wt & mut if not.
    wide      : bool
        Return wide format data frames, with one column for each mutant amino acid prediction instead of a long format table with mut and pred columns.
    make_pssm : bool
        Convert predicted frequencies to PSSM scale scores. If used with predictions other than amino acid frequencies in a N x 20 matrix with alphabetical amino acid columns this will create non-sensical data.

    Yields
    ------
    Pandas DataFrame
        Data frame of predictions, including columns for gene, position, wt, mut and prediction. If `wide=True` a column is included for each mutant amino acid prediction instead of mut and pred columns.
    """
    ind_cols = ["gene", "position", "wt"]
    if not wide:
        ind_cols.append("mut")

    if variants is not None:
        variants = variants[ind_cols]

    for seq in sequences:
        try:
            one_hot = one_hot_sequence(seq)
        except KeyError as err:
            logging.log(logging.WARN, "Skipping %s: unknown amino acid (%s)", seq.id, err)
            continue

        # Pad to be divisable
        pad_rows = padding_rows(one_hot.shape[0], layers)
        if pad_rows:
            one_hot = np.pad(one_hot, ((0, pad_rows), (0, 0)), mode='constant')
            preds = model(np.array([one_hot])).numpy()[0, :-pad_rows, :]
        else:
            preds = model(np.array([one_hot])).numpy()[0, :, :]

        if make_pssm:
            preds = freqs_to_pssm(preds)

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

def predict_proteinnet(model, data, layers=6, contacts=False, wide=False, make_pssm=False):
    """
    Generator yielding predictions from a ProteinNet Dataset

    Parameters
    ----------
    model     : keras.Model
        Sequence UNET model to predict with (or another with the same input/output signature).
    data      : ProteinNetPy Dataset
        ProteinNetPy Dataset to predict from. Use filter functions on this dataset to control what predictions are made.
    layers    : int
        Number of layers in the Sequence UNET model. All pre-trained models have 6 layers.
    contacts  : bool
        Model requires structural contact graph input.
    wide      : bool
        Return wide format data frames, with one column for each mutant amino acid prediction instead of a long format table with mut and pred columns.
    make_pssm : bool
        Convert predicted frequencies to PSSM scale scores. If used with predictions other than amino acid frequencies in a N x 20 matrix with alphabetical amino acid columns this will create non-sensical data.

    Yields
    ------
    Pandas DataFrame
        Data frame of predictions, including columns for pdb_id, chain, position, wt, mut and prediction. If `wide=True` a column is included for each mutant amino acid prediction instead of mut and pred columns.
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

        if make_pssm:
            preds = freqs_to_pssm(preds)

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
