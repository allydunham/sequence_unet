#!/usr/bin/bash
# Make predictions for each model against clinvar
# Predictions are generated for all clinvar variants for all models,
# meaning models that were trained on clinvar also predict training data

# UNET Model
python bin/predict.py --layers 6 --fasta data/clinvar/clinvar.fa models/final/freq_classifier.tf > data/clinvar/preds/unet_freq_fa.tsv

# Top Model
python bin/predict.py --layers 6 --fasta data/clinvar/clinvar.fa models/final/patho_top.tf > data/clinvar/preds/unet_freq_features_top_fa.tsv

# Finetune Model
python bin/predict.py --layers 6 --fasta data/clinvar/clinvar.fa models/final/patho_finetune.tf > data/clinvar/preds/unet_freq_finetune_fa.tsv
