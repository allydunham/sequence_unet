#!/usr/bin/env python3
"""
Module containing model for predicting variant deleteriousness trained on top of
"""
import numpy as np
from tensorflow.keras import layers, models
from metrics import CUSTOM_OBJECTS

def top_model(bottom_model=None, features=False, tune_layers=3, kernel_size=3,
              activation="sigmoid", kernel_regulariser=None,
              activity_regulariser=None, dropout=0):
    """
    Model predicting variant deleteriousness. Trained on top of a PSSM prediction model.
    Predicts a single probability for each variant.

    bottom_model: PSSM model to train on top of. Bottom_model=None creates a simple CNN model.
    features: Start from model features rather than predicted output.
    tune_layers: Number of top layers to set as trainable.
    """
    if bottom_model is not None:
        bottom_model = models.load_model(bottom_model, custom_objects=CUSTOM_OBJECTS)

        bottom_model.trainable = True

        for layer in bottom_model.layers[:-tune_layers]:
            layer.trainable = False

        x = bottom_model.layers[-2 if features else -1].output
    else:
        input_features = layers.Input(shape=[None, 20], name='input')
        x = input_features

    if dropout:
        x = layers.SpatialDropout1D(dropout, name="top_dropout")(x)

    preds = layers.Conv1D(20, kernel_size, 1, padding='same', activation=activation, name='preds',
                          kernel_regularizer=kernel_regulariser,
                          activity_regularizer=activity_regulariser)(x)

    return models.Model(inputs=bottom_model.inputs if bottom_model is not None else input_features,
                        outputs=preds)
