"""
Graph CNN layer
"""
import numpy as np
from tensorflow.keras import layers, models, activations, backend as K

def contact_graph(record, binary=False, contact_distance=1000):
    """
    Calculate a contact graph from a ProteinNetPy Record
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
    """Graph CNN Layer"""
    def __init__(self, units=32, activation='relu', **kwargs):
        super().__init__(**kwargs)
        self.units = units
        self.activation = activation
        self.activation_function = activations.get(self.activation)

    def get_config(self):
        config = super().get_config().copy()
        config.update({
            'units': self.units,
            'activation': self.activation
        })
        return config

    def build(self, input_shape):
        self.W = self.add_weight("W", shape=[input_shape[1][-1], self.units])
        self.bias = self.add_weight("bias", shape=[self.units])
        self.built = True

    def call(self, inputs):
        x = K.dot(K.batch_dot(inputs[0], inputs[1]), self.W) + self.bias
        return self.activation_function(x)