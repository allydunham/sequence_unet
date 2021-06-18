#!/usr/bin/bash
# Extract tensorboard logs to tsv files
python bin/extract_metric.py models/pssm/*/*/*/events* > data/pssm/training_logs.tsv
python bin/extract_metric.py models/classifier/*/*/*/events* > data/freq/training_logs.tsv