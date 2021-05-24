#!/usr/bin/env python3
"""
Baseline CNN to compare UNET model to.
"""
from tensorflow.keras import layers, models

def baseline_cnn(filters=32, kernel_size=7, activation="relu", pred_activation="sigmoid"):
    """
    Simple baseline CNN model predicting an Nx20 matrix

    filter: List of filter counts per layer
    kernel_size: List of kernel_sizes per layer
    activation: Activation function
    pred_activation: Activation for final layer (sigmoid for classifier and softmax for PSSM)
    """
    one_hot = layers.Input(shape=[None, 20], name='input')

    x = one_hot
    for f, k in zip(filters, kernel_size):
        x = layers.Conv1D(f, k, 1, padding='same', activation=activation)(x)

    preds = layers.Conv1D(20, 1, 1, padding='same', activation=pred_activation)(x)
    return models.Model(inputs=one_hot, outputs=preds)
