#!/usr/bin/env python3
"""
Module containing custom metrics and losses
"""
import tensorflow as tf
from tensorflow.keras import losses, metrics, backend as K

def masked_binary_crossentropy(y_true, y_pred):
    """
    Binary cross entropy where any classes marked 0 are masked and then
    true labels are offset by minus 1 (i.e a label of 1 corresponds to index 0 in y_pred).
    Accepts predictions/true values as matrices.
    """
    y_true = K.flatten(y_true)
    y_pred = K.flatten(y_pred)
    mask = tf.not_equal(y_true, 0)
    y_true_masked = tf.boolean_mask(y_true - 1, mask)
    y_pred_masked = tf.boolean_mask(y_pred, mask)
    cce = losses.binary_crossentropy(y_true_masked, y_pred_masked, from_logits=False)
    return K.mean(cce)

def masked_accuracy(y_true, y_pred):
    """
    Accuracy metric where classes marked 0 are masked and then true labels are offset by
    minus 1 (i.e a label of 1 corresponds to index 0 in y_pred).
    """
    mask = tf.not_equal(y_true, 0)
    y_true_masked = tf.boolean_mask(y_true - 1, mask)
    y_pred_masked = tf.boolean_mask(y_pred, mask)
    return metrics.binary_accuracy(y_true_masked, y_pred_masked)

CUSTOM_OBJECTS = {
    "masked_binary_crossentropy": masked_binary_crossentropy,
    "masked_accuracy": masked_accuracy
}