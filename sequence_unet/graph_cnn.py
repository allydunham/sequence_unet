"""
Simple GraphCNN TensorFlow Keras layer and supporting functions.
"""
import numpy as np
from tensorflow.keras import layers, activations, backend as K

__all__ = ['contact_graph', 'GraphCNN']

def contact_graph(record, binary=False, contact_distance=1000):
    r"""
    Calculate a residue contact graph from a ``proteinnetpy.record.ProteinNetRecord``.

    Calculate a contact matrix for the amino acids in a ProteinNet Record.
    This defines the edges of a weighted graph showing which residues connect to each other.
    The graph is determined from the proteins distance matrix, filtering to only include cells below the specified contact distance cutoff.
    Additional self connections and any missing connections to the neighbouring two residues in the primary sequence are added.
    Then each row is normalised to weight connections by inverse distance.

    Parameters
    ----------
    record           : ``proteinnetpy.record.ProteinNetRecord``
        ProteinNet Record to calculate the contact graph for.
    binary           : bool
        Use binary contacts instead of similarity scaled.
    contact_distance : numeric
        Maximum distance in nanometers between two residue for them to be considered in contact.

    Returns
    -------
    numpy.array
        A contact graph matrix, where each cell represents the similarity scaled contact between a pair of amino acids.
    """
    distance_matrix = record.distance_matrix()

    # Handle masked values
    # Set non-neighbour masked residues to high values to reflect unknown position
    # Set 2 adjacent neighbours on either side to approximate values
    mask = np.zeros_like(distance_matrix) == 1
    mask[record.mask == 0, :] = True
    mask[:, record.mask == 0] = True
    mask[np.diag_indices(len(record))] = False

    distance_matrix[mask] = 1000000

    # Immediate neighbours
    diag = np.eye(len(record), k=1) + np.eye(len(record), k=-1) == 1
    diag[~mask] = False
    distance_matrix[diag] = 380

    # Next removed neighbour
    diag = np.eye(len(record), k=2) + np.eye(len(record), k=-2) == 1
    diag[~mask] = False
    distance_matrix[diag] = 610

    # Determine neighbours (and masked assumed neighbours)
    inds = distance_matrix <= contact_distance

    # Normalised binary connections
    if binary:
        row_norm = np.linalg.norm(inds, 1, axis=1, keepdims=True)
        contacts = inds / row_norm

    # Graph connections weighted by distance within contact_distance
    else:
        contacts = 1 / (300 + distance_matrix) # Weight self roughly the same as adjacent
        contacts[~inds] = 0
        row_norm = np.linalg.norm(contacts, 1, axis=1, keepdims=True)
        contacts = contacts / row_norm

    return contacts

class GraphCNN(layers.Layer):
    """
    Graph CNN Keras Layer

    A simple GraphCNN Keras layer, which takes a feature matrix and edge matrix to calculate new features for each node of the graph.
    The matrix multiplication makes this opperation able to be performed on any sized input graph.

    Attributes
    ----------
    units               : int
        Number of features calculated in the output layer
    activation          : str
        Name of activation function applied
    activation_function : function
        Activation function applied

    Notes
    -----
    The layer implements a simple GraphCNN, which calculates new features by combining features from neighbouring nodes using the trained weights and weighting by distance to each node.

    Expected input/output:

    Input
        :Node features: (B, N, F) feature arrays with B batches of N x F matrices, representing the F features for each of the N nodes in the graph.
        :Contact graph: (B, N, N) edge weight matrix giving the weighting for each edge in the graph.

    Output
        (B, N, F) array of B batches of new N x F feature arrays.
    """
    def __init__(self, units=32, activation='relu', **kwargs):
        """
        Initialise the GraphCNN layer

        Parameters
        ----------
        record           : `proteinnetpy.record.ProteinNetRecord`
            ProteinNet Record to calculate the contact graph for.
        binary           : bool
            Use binary contacts instead of similarity scaled.
        contact_distance : numeric
            Maximum distance in nanometers between two residue for them to be considered in contact.
        """
        super().__init__(**kwargs)
        self.units = units
        self.activation = activation
        self.activation_function = activations.get(self.activation)

    def get_config(self):
        """
        Generate a configuration dictionary for serialisation

        Create the serialisation dictionary used by Keras ``save_model`` for serialising models including this layer.

        Returns
        -------
        dict
            Serialisation dictionary
        """
        config = super().get_config().copy()
        config.update({
            'units': self.units,
            'activation': self.activation
        })
        return config

    def build(self, input_shape):
        """
        Build this layer.

        Calculate and initialise layer parameters based on the input tensors.

        Parameters
        ----------
        input_shape : tuple
            Tuple of integer dimension defining the input shape.
        """
        self.W = self.add_weight("W", shape=[input_shape[1][-1], self.units])
        self.bias = self.add_weight("bias", shape=[self.units])
        self.built = True

    def call(self, inputs):
        """
        Call the layer on input tensors.

        Parameters
        ----------
        input : list
            List of input tensors, expected to contain B batches of N x F node feature matrices and paired N x N edge matrics.

        Returns
        -------
        tf.Tensor
            B x N x F output tensor with new node features.
        """
        x = K.dot(K.batch_dot(inputs[0], inputs[1]), self.W) + self.bias
        return self.activation_function(x)