#!/usr/bin/bash
# Make predictions for models against ProteinNet datasets

# Sequence UNET Frequency Classification Predictions
python bin/predict.py --layers 6 --proteinnet data/proteinnet/casp12/testing models/classifier/size/f48_k9_l6/model.tf > data/freq/unet_sequence_testing.tsv
python bin/predict.py --layers 6 --proteinnet data/proteinnet/casp12/validation models/classifier/size/f48_k9_l6/model.tf > data/freq/unet_sequence_validation.tsv

python bin/predict.py --layers 6 --contacts --proteinnet data/proteinnet/casp12/testing models/classifier/structure/elu_32/model.tf > data/freq/unet_structure_testing.tsv
python bin/predict.py --layers 6 --contacts --proteinnet data/proteinnet/casp12/validation models/classifier/structure/elu_32/model.tf > data/freq/unet_structure_validation.tsv

# Sequence UNET PSSM Predictions
python bin/predict.py --layers 6 --proteinnet data/proteinnet/casp12/testing models/pssm/size/f48_k9/model.tf > data/pssm/unet_sequence_testing.tsv
python bin/predict.py --layers 6 --proteinnet data/proteinnet/casp12/validation models/pssm/size/f48_k9/model.tf > data/pssm/unet_sequence_validation.tsv

python bin/predict.py --layers 6 --contacts --proteinnet data/proteinnet/casp12/testing models/pssm/structure/32/model.tf > data/pssm/unet_structure_testing.tsv
python bin/predict.py --layers 6 --contacts --proteinnet data/proteinnet/casp12/validation models/pssm/structure/32/model.tf > data/pssm/unet_structure_validation.tsv

# ClinVar Model Predictions
python bin/predict.py --layers 6 --proteinnet data/proteinnet/casp12/testing models/top_model/features/pred_freq_features_small/model.tf > data/clinvar/preds/pn_testing_features.tsv
python bin/predict.py --layers 6 --proteinnet data/proteinnet/casp12/testing models/top_model/finetune/classifier_sequence/model.tf > data/clinvar/preds/pn_testing_finetune.tsv