#!/usr/bin/env python3
"""
UNET model for predicting protein variant features from sequence and structure.
"""
import numpy as np
from tensorflow.keras import layers, models
from graph_cnn import GraphCNN

def sequence_unet(filters=8, kernel_size=5, num_layers=4, dropout=0,
                   graph_layers=None, graph_activation="relu",
                   conv_activation="relu", pred_activation="sigmoid",
                   kernel_regulariser=None, batch_normalisation=False):
    """
    Keras model predicting multiple alignment frequency matrix from sequences, based on
    a 1D UNET architecture.

    Expected input: (B, N, 20) arrays with B batches of Nx20 matrices, representing
                    sequences with one hot encoding on each row. An additional (B, N, N)
                    array with NxN contact graphs is required when graph_layers is not none.
                    The UNET architecture means N must be divisble by 2^(num_layers - 1).
                    This can be achieved with padding and masking.
    Output:         (B, N, 20) arrays of B batches of Nx20 arrays, with
                    the predicted features of each AA at each position as rows.
    """
    num_layers = num_layers - 1 # 0 index layers

    input_seq = layers.Input(shape=[None, 20], name='input_seq')
    inputs = input_seq
    x = input_seq

    if graph_layers is not None:
        input_graph = layers.Input(shape=[None, None], name='input_graph')
        inputs = [input_seq, input_graph]
        for i, units in enumerate(graph_layers):
            x = GraphCNN(units, activation=graph_activation, name=f'graph_{i}')([input_graph, x])

        x = layers.concatenate([input_seq, x], name="graph_concat")

    # Contraction
    contraction = []
    for i in range(num_layers):
        x = layers.Conv1D(filters * 2 ** i, kernel_size, 1, padding='same',
                          activation=conv_activation, name=f'down_{i}_conv_1',
                          kernel_regularizer=kernel_regulariser)(x)
        x = layers.Conv1D(filters * 2 ** i, kernel_size, 1, padding='same',
                          activation=conv_activation, name=f'down_{i}_conv_2',
                          kernel_regularizer=kernel_regulariser)(x)
        contraction.append(x)
        x = layers.MaxPool1D(pool_size=2, strides=2, padding='same', name=f'down_{i}_max_pool')(x)

        if dropout:
            x = layers.SpatialDropout1D(dropout, name=f'down_{i}_dropout')(x)

        if batch_normalisation:
            x = layers.BatchNormalization(name=f"down_{i}_batch_normalisation")(x)

    # Bottom layer
    x = layers.Conv1D(filters * 2 ** num_layers, kernel_size, 1,
                      padding='same', activation=conv_activation, name=f'bottom_conv_1',
                      kernel_regularizer=kernel_regulariser)(x)
    x = layers.Conv1D(filters * 2 ** num_layers, kernel_size, 1,
                      padding='same', activation=conv_activation, name=f'bottom_conv_2',
                      kernel_regularizer=kernel_regulariser)(x)

    if dropout:
        x = layers.SpatialDropout1D(dropout, name='bottom_dropout')(x)

    if batch_normalisation:
            x = layers.BatchNormalization(name="bottom_batch_normalisation")(x)

    # Expansion
    for i in range(num_layers - 1, -1, -1):
        x = layers.UpSampling1D(2, name=f'up_{i}_upsample')(x)
        x = layers.Conv1D(filters * 2 ** i, kernel_size, 1, padding='same',
                          activation=conv_activation, name=f'up_{i}_conv_1',
                          kernel_regularizer=kernel_regulariser)(x)
        x = layers.concatenate([contraction[i], x], name=f"up_{i}_concat")

        x = layers.Conv1D(filters * 2 ** i, kernel_size, 1, padding='same',
                          activation=conv_activation, name=f'up_{i}_conv_2',
                          kernel_regularizer=kernel_regulariser)(x)
        x = layers.Conv1D(filters * 2 ** i, kernel_size, 1, padding='same',
                          activation=conv_activation, name=f'up_{i}_conv_3',
                          kernel_regularizer=kernel_regulariser)(x)

        if dropout:
            x = layers.SpatialDropout1D(dropout, name=f'up_{i}_dropout')(x)

        if batch_normalisation:
            x = layers.BatchNormalization(name=f"up_{i}_batch_normalisation")(x)

    preds = layers.Conv1D(20, 1, 1, activation=pred_activation, name="predictor",
                          kernel_regularizer=kernel_regulariser)(x)

    return models.Model(inputs=inputs, outputs=preds)

