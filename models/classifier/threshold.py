"""
Experiment testing various deleteriousness thresholds for the Sequence UNET model
"""
from operator import ne
import os
import sys

import utils
from tensorflow.keras import optimizers

from proteinnetpy.data import ProteinNetDataset, ProteinNetMap
from proteinnetpy.data import make_length_filter

import metrics
import pn_maps
from sequence_unet import sequence_unet

def make_load_data(threshold=0.01):
    """
    Generate load_data function
    """
    def f(validation=False):
        """
        Input data for PSSM top model
        """
        if validation:
            pn_path = 'data/proteinnet/casp12/validation'
        else:
            pn_path = 'data/proteinnet/casp12/training_95'

        filter_func = make_length_filter(min_length=32, max_length=2000)
        data = ProteinNetDataset(path=pn_path, preload=False, filter_func=filter_func)
        func = pn_maps.SequenceUNETMapFunction(num_layers=6, threshold=threshold)
        return ProteinNetMap(data, func=func, static=True, filter_errors=True)
    return f

def main():
    """Main script"""
    root = 'models/classifier/threshold'
    if not os.path.isdir(root):
        os.mkdir(root)

    # threshold, pos_weight, neg_weight
    thresholds = [
        (0.1, 0.15, 1),
        (0.01, 1, 0.9),
        (0.001, 1, 0.38),
        (0.0001, 1, 0.25)
    ]

    for threshold, pos_weight, neg_weight in thresholds:
        model_dir = f"{root}/t{threshold}"

        if os.path.isdir(model_dir):
            print(f"Model {model_dir} already exists, skipping", file=sys.stderr)
            continue

        model = sequence_unet(filters=48, kernel_size=9, num_layers=6,
                              batch_normalisation=True, dropout=0.05, conv_activation="relu")

        optimiser = optimizers.Adam(lr=0.01, epsilon=0.01)
        loss = metrics.WeightedMaskedBinaryCrossEntropy(pos_weight=pos_weight,
                                                        neg_weight=neg_weight)
        acc = metrics.masked_accuracy
        model.compile(optimizer=optimiser, loss=loss, metrics=[acc])

        load_data = make_load_data(threshold=threshold)

        # Create sample train script
        command = utils.model_bsub(f"threshold_{threshold}", model_dir,
                                    ram=10000, epochs=150, validation_epochs=1,
                                    checkpoint=None, big_job=True, save_format='tf')

        # Use this to setup a model directory for the experiment(s)
        utils.make_experiment_dir(model, model_dir, load_data, command, save_format='tf')

if __name__ == "__main__":
    # No argparse as these scripts serve as the config for experiments
    main()
