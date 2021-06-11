#!/usr/bin/bash
# Make predictions for models against DMS and Jelier results

# DMS Predictions
python bin/predict.py --layers 6 --fasta data/dms/dms.fa models/classifier/size/f48_k9_l6/model.tf > data/dms/preds/unet_freq.tsv

python bin/predict.py --layers 6 --fasta data/dms/dms.fa models/top_model/features/pred_freq_features_small/model.tf > data/dms/preds/clinvar_top.tsv

python bin/predict.py --layers 0 --fasta data/dms/dms.fa models/baseline/freq_classifier/single/model.h5 > data/dms/preds/baseline_freq.tsv

python bin/predict.py --layers 0 --fasta data/dms/dms.fa models/baseline/clinvar/single/model.h5 > data/dms/preds/baseline_clinvar.tsv


# Jelier Predictions
python bin/predict.py --layers 6 --fasta data/jelier/variant_genes.fa --tsv data/jelier/jelier_variants.tsv models/classifier/size/f48_k9_l6/model.tf > data/jelier/preds/unet_freq.tsv

python bin/predict.py --layers 6 --fasta data/jelier/variant_genes.fa --tsv data/jelier/jelier_variants.tsv models/top_model/features/pred_freq_features_small/model.tf > data/jelier/preds/clinvar_top.tsv

python bin/predict.py --layers 0 --fasta data/jelier/variant_genes.fa --tsv data/jelier/jelier_variants.tsv models/baseline/freq_classifier/single/model.h5 > data/jelier/preds/baseline_freq.tsv

python bin/predict.py --layers 0 --fasta data/jelier/variant_genes.fa --tsv data/jelier/jelier_variants.tsv models/baseline/clinvar/single/model.h5 > data/jelier/preds/baseline_clinvar.tsv