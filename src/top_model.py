#!/usr/bin/env python3
"""
Module containing model for predicting variant deleteriousness trained on top of
"""
import numpy as np
from tensorflow.keras import layers, models
from metrics import CUSTOM_OBJECTS

def top_model(bottom_model=None, pssm=False, tune_layers=3, kernel_size=3,
              activation="sigmoid", kernel_regulariser=None,
              activity_regulariser=None, dropout=0):
    """
    Model predicting variant deleteriousness. Trained on top of a PSSM prediction model.
    Predicts a single probability for each variant.

    bottom_model: PSSM model to train on top of
    pssm: Start from PSSM values. bottom_model=None and pssm=True creates a model that expects
    a PSSM as input. If bottom_model is specified and pssm=True the predicted PSSM is used, and
    if pssm=False the layer below is used.
    tune_layers: Number of top layers to set as trainable.
    """
    if bottom_model is not None:
        bottom_model = models.load_model(bottom_model, custom_objects=CUSTOM_OBJECTS)

        bottom_model.trainable = True

        for layer in bottom_model.layers[:-tune_layers]:
            layer.trainable = False

        x = bottom_model.layers[-1 if pssm else -2].output
    elif pssm:
        input_pssm = layers.Input(shape=[None, 20], name='input_pssm')
        x = input_pssm
    else:
        raise ValueError("Either bottom or pssm must be specified")

    if dropout:
        x = layers.SpatialDropout1D(dropout, name="top_dropout")(x)

    preds = layers.Conv1D(20, kernel_size, 1, padding='same', activation=activation, name='preds',
                          kernel_regularizer=kernel_regulariser,
                          activity_regularizer=activity_regulariser)(x)

    return models.Model(inputs=bottom_model.inputs if bottom_model is not None else input_pssm,
                        outputs=preds)
