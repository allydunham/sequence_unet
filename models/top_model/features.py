"""
Experiment testing various features for a Sequence UNET top model
"""
import os
import sys

import utils
from tensorflow.keras import optimizers
import pandas as pd

from proteinnetpy.data import ProteinNetDataset, ProteinNetMap
from proteinnetpy.data import make_length_filter, combine_filters, make_id_filter

import metrics
import pn_maps
from top_model import top_model

def make_load_data(structure=False, pssm=False, num_layers=6):
    """
    Generate a load_data function
    """
    def f(validation=False):
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
        func = pn_maps.ClinVarMapFunction(clinvar=clinvar, num_layers=num_layers,
                                          contact_graph=structure, pssm=pssm)
        return ProteinNetMap(data, func=func, static=True, filter_errors=True)
    return f

def main():
    """Main script"""
    root = 'models/top_model/features'
    if not os.path.isdir(root):
        os.mkdir(root)

    freq_model = "models/classifier/size/f64_k9_l6/model.tf"
    freq_struct_model = "models/classifier/structure/elu_32/model.tf"

    pssm_model = "models/pssm/size/f64_k9/model.tf"
    pssm_struct_model = "models/pssm/structure/32/model.tf"

    # name: args
    features = {
        "true_pssm": {"bottom_model": None, "kernel_size": 3},
        "true_pssm_small": {"bottom_model": None, "kernel_size": 1},
        "true_seq": {"bottom_model": None, "kernel_size": 3},
        "true_seq_small": {"bottom_model": None, "kernel_size": 1},

        "pred_freq": {"bottom_model": freq_model, "tune_layers": 1,
                      "kernel_size": 3, "features": False},
        "pred_freq_features": {"bottom_model": freq_model, "tune_layers": 1,
                               "kernel_size": 3, "features": True},
        "pred_freq_features_small": {"bottom_model": freq_model, "tune_layers": 0,
                                     "kernel_size": 1, "features": True},

        "pred_freq_structure": {"bottom_model": freq_struct_model, "tune_layers": 1,
                                "kernel_size": 3, "features": False},
        "pred_freq_structure_features": {"bottom_model": freq_struct_model, "tune_layers": 1,
                                         "kernel_size": 3, "features": True},
        "pred_freq_structure_features_small": {"bottom_model": freq_struct_model, "tune_layers": 0,
                                               "kernel_size": 1, "features": True},

        "pred_pssm": {"bottom_model": pssm_model, "tune_layers": 1,
                      "kernel_size": 3, "features": False},
        "pred_pssm_features": {"bottom_model": pssm_model, "tune_layers": 1,
                               "kernel_size": 3, "features": True},
        "pred_pssm_features_small": {"bottom_model": pssm_model, "tune_layers": 0,
                                     "kernel_size": 1, "features": True},

        "pred_pssm_structure": {"bottom_model": pssm_struct_model, "tune_layers": 1,
                                "kernel_size": 3, "features": False},
        "pred_pssm_structure_features": {"bottom_model": pssm_struct_model, "tune_layers": 1,
                                         "kernel_size": 3, "features": True},
        "pred_pssm_structure_features_small": {"bottom_model": pssm_struct_model, "tune_layers": 0,
                                               "kernel_size": 1, "features": True}
    }

    for name, args in features.items():
        model_dir = f"{root}/{name}"

        if os.path.isdir(model_dir):
            print(f"Model {model_dir} already exists, skipping", file=sys.stderr)
            continue

        model = top_model(**args, kernel_regulariser="l2", dropout=0.2)

        optimiser = optimizers.Adam(lr=0.01, epsilon=0.01)
        loss = metrics.masked_binary_crossentropy
        acc = metrics.masked_accuracy
        model.compile(optimizer=optimiser, loss=loss, metrics=[acc])

        load_data = make_load_data(structure="structure" in name, pssm=name == "true_pssm")

        # Create sample train script
        finetune = None if args["bottom_model"] is None else 2 + args["tune_layers"]
        command = utils.model_bsub(f"features_{name}", model_dir,
                                    ram=15000, epochs=150, validation_epochs=1,
                                    checkpoint=None, big_job=True, save_format='tf',
                                    finetune=finetune, early_stop=10)

        # Use this to setup a model directory for the experiment(s)
        utils.make_experiment_dir(model, model_dir, load_data, command, save_format='tf')

if __name__ == "__main__":
    # No argparse as these scripts serve as the config for experiments
    main()
