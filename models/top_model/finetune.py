"""
Experiment testing various features for a Sequence UNET top model
"""
import os
import sys

import utils
from tensorflow.keras import optimizers, models
import pandas as pd

from proteinnetpy.data import ProteinNetDataset, ProteinNetMap
from proteinnetpy.data import make_length_filter, combine_filters, make_id_filter

import metrics
import pn_maps
from top_model import top_model

def load_data(validation=False):
    """
    Input data for PSSM top model
    """
    if validation:
        clinvar_path = 'data/clinvar/clinvar_val.tsv'
    else:
        clinvar_path = 'data/clinvar/clinvar_train.tsv'
    clinvar = pd.read_csv(clinvar_path, sep="\t")

    pn_path = 'data/proteinnet/casp12/clinvar'
    filter_func = combine_filters(
        make_length_filter(min_length=32, max_length=2000),
        make_id_filter(list(clinvar.pdb_id), list(clinvar.chain))
    )
    data = ProteinNetDataset(path=pn_path, preload=False, filter_func=filter_func)
    func = pn_maps.ClinVarMapFunction(clinvar=clinvar, num_layers=6,
                                      contact_graph=True, pssm=False)
    return ProteinNetMap(data, func=func, static=True, filter_errors=True)

def main():
    """Main script"""
    model_dir = 'models/top_model/finetune'

    if os.path.isdir(model_dir):
        raise FileExistsError(f"Model {model_dir} already exists, skipping")

    if not os.path.isdir(model_dir):
        os.mkdir(model_dir)

    model_path = "models/classifier/structure/elu_32/model.tf"
    model = models.load_model(model_path, custom_objects=metrics.CUSTOM_OBJECTS)

    optimiser = optimizers.Adam(lr=0.001, epsilon=0.01)
    loss = metrics.masked_binary_crossentropy
    acc = metrics.masked_accuracy
    model.compile(optimizer=optimiser, loss=loss, metrics=[acc])

    # Create sample train script
    command = utils.model_bsub("finetune", model_dir, ram=10000, epochs=50,
                               validation_epochs=1, checkpoint=None, big_job=True,
                               save_format='tf', finetune=2, early_stop=10)

    # Use this to setup a model directory for the experiment(s)
    utils.make_experiment_dir(model, model_dir, load_data, command, save_format='tf')

if __name__ == "__main__":
    # No argparse as these scripts serve as the config for experiments
    main()
