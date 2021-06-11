#!/usr/bin/bash
# Make predictions for each model against clinvar
# Predictions are generated for all clinvar variants for all models,
# meaning models that were trained on clinvar also predict training data

python bin/predict.py --layers 6 --proteinnet data/proteinnet/casp12/clinvar --clinvar --tsv data/clinvar/clinvar_variants.tsv models/classifier/size/f48_k9_l6/model.tf > data/clinvar/preds/unet_freq.tsv

python bin/predict.py --layers 6 --contacts --proteinnet data/proteinnet/casp12/clinvar --clinvar --tsv data/clinvar/clinvar_variants.tsv models/classifier/structure/elu_32/model.tf > data/clinvar/preds/unet_freq_structure.tsv

python bin/predict.py --layers 6 --proteinnet data/proteinnet/casp12/clinvar --clinvar --tsv data/clinvar/clinvar_variants.tsv models/top_model/features/pred_freq_features_small/model.tf > data/clinvar/preds/unet_freq_features_top.tsv

python bin/predict.py --layers 6 --contacts --proteinnet data/proteinnet/casp12/clinvar --clinvar --tsv data/clinvar/clinvar_variants.tsv models/top_model/features/pred_freq_structure_features_small/model.tf > data/clinvar/preds/unet_freq_structure_features_top.tsv

python bin/predict.py --layers 0 --proteinnet data/proteinnet/casp12/clinvar --clinvar --tsv data/clinvar/clinvar_variants.tsv models/baseline/freq_classifier/single/model.h5 > data/clinvar/preds/baseline_freq.tsv

python bin/predict.py --layers 0 --proteinnet data/proteinnet/casp12/clinvar --clinvar --tsv data/clinvar/clinvar_variants.tsv models/baseline/clinvar/single/model.h5 > data/clinvar/preds/baseline_clinvar.tsv