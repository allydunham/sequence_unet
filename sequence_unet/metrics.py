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
    y_true = tf.cast(y_true, y_pred.dtype)
    y_true_masked = tf.boolean_mask(y_true - 1, mask)
    y_pred_masked = tf.boolean_mask(y_pred, mask)
    ce = K.binary_crossentropy(y_true_masked, y_pred_masked, from_logits=False)
    return K.mean(ce)

def masked_accuracy(y_true, y_pred):
    """
    Accuracy metric where classes marked 0 are masked and then true labels are offset by
    minus 1 (i.e a label of 1 corresponds to index 0 in y_pred).
    """
    mask = tf.not_equal(y_true, 0)
    y_true_masked = tf.boolean_mask(y_true - 1, mask)
    y_pred_masked = tf.boolean_mask(y_pred, mask)
    return metrics.binary_accuracy(y_true_masked, y_pred_masked)

class WeightedMaskedBinaryCrossEntropy(losses.Loss):
    """
    Binary cross entropy where any classes marked 0 are masked and then
    true labels are offset by minus 1 (i.e a label of 1 corresponds to index 0 in y_pred).
    Output is weighted by the given class weights.
    Accepts predictions/true values as matrices.
    """
    def __init__(self, pos_weight, neg_weight, from_logits=False,
                 reduction=losses.Reduction.AUTO,
                 name='weighted_masked_binary_crossentropy'):
        super().__init__(reduction=reduction, name=name)
        self.pos_weight = pos_weight
        self.neg_weight = neg_weight
        self.from_logits = from_logits

    def call(self, y_true, y_pred):
        y_true = K.flatten(y_true)
        y_pred = K.flatten(y_pred)
        mask = tf.not_equal(y_true, 0)
        y_true = tf.cast(y_true, y_pred.dtype)
        y_true_masked = tf.boolean_mask(y_true - 1, mask)
        y_pred_masked = tf.boolean_mask(y_pred, mask)
        ce = K.binary_crossentropy(y_true_masked, y_pred_masked, from_logits=self.from_logits)
        weights = y_pred_masked * self.pos_weight + (1. - y_pred_masked) * self.neg_weight
        return ce * weights

    def get_config(self):
        config = super().get_config()
        config.update({"pos_weight": self.pos_weight,
                       "neg_weight": self.neg_weight,
                       "from_logits": self.from_logits})
        return config

CUSTOM_OBJECTS = {
    "masked_binary_crossentropy": masked_binary_crossentropy,
    "masked_accuracy": masked_accuracy,
    "WeightedMaskedBinaryCrossEntropy": WeightedMaskedBinaryCrossEntropy
}
