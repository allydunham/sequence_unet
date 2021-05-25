"""
Baseline CNN models of various sizes to compare Sequence UNET models to
"""
import os
import sys
import pandas as pd

import utils
from tensorflow.keras import optimizers

from proteinnetpy.data import ProteinNetDataset, ProteinNetMap
from proteinnetpy.data import make_length_filter, make_id_filter, combine_filters

import metrics
import pn_maps
from baseline_cnn import baseline_cnn

def make_load_data_proteinnet(pssm=False):
    """
    Make ProteinNet load_data function
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
        func = pn_maps.SequenceUNETMapFunction(num_layers=1, threshold=None if pssm else 0.01)
        return ProteinNetMap(data, func=func, static=True, filter_errors=True)
    return f

def load_data_clinvar(validation=False):
    """
    Input data for PSSM top model
    """
    if validation:
        clin_path = 'data/clinvar/clinvar_val.tsv'
    else:
        clin_path = 'data/clinvar/clinvar_train.tsv'

    clinvar = pd.read_csv(clin_path, sep="\t")
    pn_path = 'data/protein_net/text/casp12/clinvar'
    filter_func = combine_filters(
        make_length_filter(min_length=32, max_length=2000),
        make_id_filter(list(clinvar.pdb_id), list(clinvar.chain))
    )
    data = ProteinNetDataset(path=pn_path, preload=False, filter_func=filter_func)
    func = pn_maps.ClinVarMapFunction(clinvar=clinvar, num_layers=1,
                                      contact_graph=False, pssm=False)
    return ProteinNetMap(data, func=func, static=True, filter_errors=True)

def main():
    """Main script"""
    root = 'models/baseline'

    # model, function
    models = {
        "clinvar": (False, load_data_clinvar),
        "freq_classifier": (False, make_load_data_proteinnet(pssm=False)),
        "pssm": (True, make_load_data_proteinnet(pssm=True))
    }

    sizes = {
        "single": ([32], [7]),
        "double": ([32]*2, [7]*2),
        "large": ([32, 64], [7, 7]),
    }

    for model_name, (pssm, load_data) in models.items():
        try:
            os.mkdir(f"{root}/{model_name}/")
        except FileExistsError:
            pass

        for size_name, (filters, kernels) in sizes.items():
            model_dir = f"{root}/{model_name}/{size_name}"

            if os.path.isdir(model_dir):
                print(f"Model {model_dir} already exists, skipping", file=sys.stderr)
                continue

            model = baseline_cnn(filters=filters, kernel_size=kernels,
                                 pred_activation="softmax" if pssm else "sigmoid")

            optimiser = optimizers.Adam(lr=0.01, epsilon=0.01)
            loss = 'kullback_leibler_divergence' if pssm else metrics.masked_binary_crossentropy
            acc = 'mean_absolute_error' if pssm else metrics.masked_accuracy
            model.compile(optimizer=optimiser, loss=loss, metrics=[acc])

            # Create sample train script
            command = utils.model_bsub(f"baseline_{model_name}_{size_name}", model_dir,
                                        ram=8000, epochs=100, validation_epochs=1,
                                        big_job=False, save_format='h5')

            # Use this to setup a model directory for the experiment(s)
            utils.make_experiment_dir(model, model_dir, load_data, command, save_format='h5')

if __name__ == "__main__":
    # No argparse as these scripts serve as the config for experiments
    main()
