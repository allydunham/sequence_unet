#!/usr/bin/env python3
"""
Custom weighted metrics and losses for training and evaluating Sequence UNET models.
"""
import tensorflow as tf
from tensorflow.keras import losses, metrics, backend as K

__all__ = ["masked_binary_crossentropy", "masked_accuracy", "WeightedMaskedBinaryCrossEntropy"]

def masked_binary_crossentropy(y_true, y_pred):
    """
    Zero masked binary cross-entropy.

    A version of the binary cross-entropy loss function with masked labels. Classes labelled 0 in y_true are masked, those labelled 1 correspond to 0 in y_pred and those labelled 2 to 1 in y_pred (i.e. offset by -1). Accepts predictions/true values as matrices.

    Parameters
    ----------
    y_true            : float
        True class labels. 0 < x < 1.
    y_pred            : int
        Predicted class labels. x = 0,1,2 with 0 being masked and 1,2 converted to 0,1.

    Returns
    -------
    float
        Binary cross-entropy from non-masked positions
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
    Zero masked binary accuracy.

    A version of the binary accuracy metric with masked labels. Classes labelled 0 in y_true are masked, those labelled 1 correspond to 0 in y_pred and those labelled 2 to 1 in y_pred (i.e. offset by -1). Accepts predictions/true values as matrices.

    Parameters
    ----------
    y_true            : float
        True class labels. 0 < x < 1.
    y_pred            : int
        Predicted class labels. x = 0,1,2 with 0 being masked and 1,2 converted to 0,1.

    Returns
    -------
    float
        Binary accuracy from non-masked positions
    """
    mask = tf.not_equal(y_true, 0)
    y_true_masked = tf.boolean_mask(y_true - 1, mask)
    y_pred_masked = tf.boolean_mask(y_pred, mask)
    return metrics.binary_accuracy(y_true_masked, y_pred_masked)

class WeightedMaskedBinaryCrossEntropy(losses.Loss):
    """
    Version of `masked_binary_crossentropy` with class weights.

    A version of the binary cross-entropy loss function with masked labels and class weights. Classes labelled 0 in y_true are masked, those labelled 1 correspond to 0 in y_pred and those labelled 2 to 1 in y_pred (i.e. offset by -1). Class weightings are applied after masking.Accepts predictions/true values as matrices.

    Attributes
    ----------
    pos_weight  : float
        Weight applied to positvely labelled (2) items.
    neg_weight  : float
        Weight applied to negatively labelled (1) items.
    from_logits : bool
        Predictions are assumed to be in logit rather than probability form.
    """
    def __init__(self, pos_weight, neg_weight, from_logits=False,
                 reduction=losses.Reduction.AUTO,
                 name='weighted_masked_binary_crossentropy'):
        """
        Initialise a `WeightedMaskedBinaryCrossEntropy` loss

        Parameters
        ----------
        pos_weight   : float
            Weight for positive class (labelled 2 in y_true).
        neg_weight   : float
            Weight for negative class (labelled 1 in y_true).
        from_logits  : bool
            Treat y_pred as logits rather than probabilities.
        reduction    : Keras loss reduction strategy
            Unused, included for keras compatibility.
        name         : str
            Metric name.
        """
        super().__init__(reduction=reduction, name=name)
        self.pos_weight = pos_weight
        self.neg_weight = neg_weight
        self.from_logits = from_logits

    def call(self, y_true, y_pred):
        """
        Calculate Weighted masked binary cross-entropy.

        Parameters
        ----------
        y_true            : float
            True class labels. 0 < x < 1.
        y_pred            : int
            Predicted class labels. x = 0,1,2 with 0 being masked and 1,2 converted to 0,1.

        Returns
        -------
        float
            Binary cross-entropy from non-masked positions
        """
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
        """
        Generate configuration dictionary for serialisation.

        Generate configuration dictionary used by the Keras save_model function for serialisation.

        Returns
        -------
        dict
            Configuration dictionary.
        """
        config = super().get_config()
        config.update({"pos_weight": self.pos_weight,
                       "neg_weight": self.neg_weight,
                       "from_logits": self.from_logits})
        return config
