"""
Experiment testing various Sequence UNET model sizes for PSSM prediction
"""
import os
import sys

import utils
from tensorflow.keras import optimizers

from proteinnetpy.data import ProteinNetDataset, ProteinNetMap
from proteinnetpy.data import make_length_filter

import metrics
import pn_maps
from seq_unet import sequence_unet

def load_data(validation=False):
    """
    Input data for PSSM top model
    """
    if validation:
        pn_path = 'data/proteinnet/casp12/validation'
    else:
        pn_path = 'data/proteinnet/casp12/training_95'

    filter_func = make_length_filter(min_length=32, max_length=2000)
    data = ProteinNetDataset(path=pn_path, preload=False, filter_func=filter_func)
    func = pn_maps.SequenceUNETMapFunction(num_layers=6, threshold=None,
                                            weights=True)
    return ProteinNetMap(data, func=func, static=True, filter_errors=True)

def main():
    """Main script"""
    root = 'models/pssm/size'
    if not os.path.isdir(root):
        os.mkdir(root)

    # filters, kernel
    size = (
        (32, 5),
        (48, 5),
        (32, 7),
        (32, 9),
        (32, 5),
        (32, 7),
        (48, 9),
        (64, 9),
        (96, 9)
    )

    for filters, kernel in size:
        model_dir = f"{root}/f{filters}_k{kernel}"

        if os.path.isdir(model_dir):
            print(f"Model {model_dir} already exists, skipping", file=sys.stderr)
            continue

        model = sequence_unet(filters=filters, kernel_size=kernel, num_layers=6,
                              batch_normalisation=True, dropout=0.05, conv_activation="relu",
                              pred_activation="softmax")

        optimiser = optimizers.Adam(lr=0.01, epsilon=0.01)
        model.compile(optimizer=optimiser, loss="kullback_leibler_divergence",
                      metrics=["mean_absolute_error"], sample_weight_mode="temporal")

        # Create sample train script
        ram = 10000 if filters < 50 else 15000
        command = utils.model_bsub(f"size_f{filters}_k{kernel}", model_dir,
                                    ram=ram, epochs=150, validation_epochs=1,
                                    checkpoint=None, big_job=True, save_format='tf')

        # Use this to setup a model directory for the experiment(s)
        utils.make_experiment_dir(model, model_dir, load_data, command, save_format='tf')

if __name__ == "__main__":
    # No argparse as these scripts serve as the config for experiments
    main()
