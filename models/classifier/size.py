"""
Experiment testing various Sequence UNET model sizes
"""
import os
import sys

import utils
from tensorflow.keras import optimizers

from proteinnetpy.data import ProteinNetDataset, ProteinNetMap
from proteinnetpy.data import make_length_filter

import metrics
import pn_maps
from sequence_unet import sequence_unet

def make_load_data(num_layers=4):
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
        func = pn_maps.SequenceUNETMapFunction(num_layers=num_layers, threshold=0.01)
        return ProteinNetMap(data, func=func, static=True, filter_errors=True)
    return f

def main():
    """Main script"""
    root = 'models/classifier/size'
    if not os.path.isdir(root):
        os.mkdir(root)

    # filters, kernel, layers
    size = (
        (32, 5, 2),
        (16, 5, 4),
        (32, 5, 4),
        (48, 5, 4),
        (16, 7, 4),
        (16, 9, 4),
        (16, 5, 6),
        (32, 3, 6),
        (32, 5, 6),
        (32, 7, 6),
        (32, 9, 6),
        (48, 9, 6),
        (64, 9, 4),
    )

    for filters, kernel, layers in size:
        model_dir = f"{root}/f{filters}_k{kernel}_l{layers}"

        if os.path.isdir(model_dir):
            print(f"Model {model_dir} already exists, skipping", file=sys.stderr)
            continue

        model = sequence_unet(filters=filters, kernel_size=kernel, num_layers=layers,
                              batch_normalisation=True, dropout=0.05, conv_activation="elu")

        optimiser = optimizers.Adam(lr=0.01, epsilon=0.01)
        loss = metrics.masked_binary_crossentropy
        acc = metrics.masked_accuracy
        model.compile(optimizer=optimiser, loss=loss, metrics=[acc])

        load_data = make_load_data(num_layers=layers)

        # Create sample train script
        command = utils.model_bsub(f"size_f{filters}_k{kernel}_l{layers}", model_dir,
                                    ram=10000, epochs=150, validation_epochs=1,
                                    checkpoint=None, big_job=True, save_format='tf')

        # Use this to setup a model directory for the experiment(s)
        utils.make_experiment_dir(model, model_dir, load_data, command, save_format='tf')

if __name__ == "__main__":
    # No argparse as these scripts serve as the config for experiments
    main()
