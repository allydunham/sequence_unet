"""
Experiment testing various regularisations on the Sequence UNET model
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
    func = pn_maps.SequenceUNETMapFunction(num_layers=6, threshold=0.01)
    return ProteinNetMap(data, func=func, static=True, filter_errors=True)

def main():
    """Main script"""
    root = 'models/classifier/regularisation'
    if not os.path.isdir(root):
        os.mkdir(root)

    # dropout, kernel, batch
    regularisation = (
        (0, None, False),
        (0.05, None, False),
        (0.1, None, False),
        (0, "l2", False),
        (0, None, True),
        (0.05, "l2", False),
        (0.05, None, True),
        (0.05, "l2", True),
        (0, "l2", True),
    )

    for dropout, kernel, batch in regularisation:
        model_dir = f"{root}/d{dropout}_{kernel}_{batch}"

        if os.path.isdir(model_dir):
            print(f"Model {model_dir} already exists, skipping", file=sys.stderr)
            continue

        model = sequence_unet(filters=32, kernel_size=5, num_layers=6, dropout=dropout,
                              kernel_regulariser=kernel, batch_normalisation=batch)

        optimiser = optimizers.Adam(lr=0.01, epsilon=0.01)
        loss = metrics.masked_binary_crossentropy
        acc = metrics.masked_accuracy
        model.compile(optimizer=optimiser, loss=loss, metrics=[acc])

        # Create sample train script
        command = utils.model_bsub(f"reg_d{dropout}_{kernel}_{batch}", model_dir,
                                    ram=8000, epochs=150, validation_epochs=1,
                                    checkpoint=None, big_job=True, save_format='tf')

        # Use this to setup a model directory for the experiment(s)
        utils.make_experiment_dir(model, model_dir, load_data, command, save_format='tf')

if __name__ == "__main__":
    # No argparse as these scripts serve as the config for experiments
    main()
