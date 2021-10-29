"""
Experiment testing Sequence UNET representation learning models, with and without structure
"""
import os
import sys

import utils
from tensorflow.keras import optimizers

from proteinnetpy.data import ProteinNetDataset, ProteinNetMap
from proteinnetpy.data import make_length_filter

import pn_maps
from seq_unet import sequence_unet

def make_load_data(structure=False):
    """
    Generate load_data function
    """
    def f(validation=False):
        """
        Input data for Position Representation model
        """
        if validation:
            pn_path = 'data/proteinnet/casp12/validation'
        else:
            pn_path = 'data/proteinnet/casp12/training_95'

        filter_func = make_length_filter(min_length=32, max_length=2000)
        data = ProteinNetDataset(path=pn_path, preload=True, filter_func=filter_func)
        func = pn_maps.MaskedSequenceMapFunction(num_layers=6, mask=0.05, contact_graph=structure)
        return ProteinNetMap(data, func=func, static=False, filter_errors=True)
    return f

def main():
    """Main script"""
    root = 'models/representation/structure'
    if not os.path.isdir(root):
        os.mkdir(root)

    # graph_layers
    structure = (
        None,
        # [32],
        # [32,32],
        [64]
        # [64,64]
    )

    for graph_layers in structure:
        layer_str = '_'.join(str(i) for i in graph_layers) if graph_layers is not None else 'none'
        model_dir = f"{root}/{layer_str}"

        if os.path.isdir(model_dir):
            print(f"Model {model_dir} already exists, skipping", file=sys.stderr)
            continue

        model = sequence_unet(filters=64, kernel_size=9, num_layers=6, conv_activation="relu",
                              graph_layers=graph_layers, graph_activation="elu",
                              batch_normalisation=True, dropout=0.05)

        optimiser = optimizers.Adam(lr=0.01, epsilon=0.01)
        model.compile(optimizer=optimiser, loss="categorical_crossentropy",
                      metrics=["accuracy"], sample_weight_mode="temporal")

        load_data = make_load_data(structure=graph_layers is not None)

        # Create sample train script
        command = utils.model_bsub(f"structure_{layer_str}", model_dir,
                                    ram=15000, epochs=150, validation_epochs=5,
                                    checkpoint=None, big_job=True, save_format='tf',
                                    early_stop=20)

        # Use this to setup a model directory for the experiment(s)
        utils.make_experiment_dir(model, model_dir, load_data, command, save_format='tf')

if __name__ == "__main__":
    # No argparse as these scripts serve as the config for experiments
    main()
