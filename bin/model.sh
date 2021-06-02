#!/usr/bin/env bash
# Reset an initialised model directory

model_dir=$1
rm "$model_dir/train/*" "$model_dir/validation/*"

if test -f "$model_dir/model.tf"; then
    rm -r "$model_dir/model.tf"
    cp -r "$model_dir/initial_model.tf" "$model_dir/model.tf"
elif test -f "$model_dir/model.h5"; then
    rm "$model_dir/model.h5"
    cp -r "$model_dir/initial_model.h5" "$model_dir/model.h5"
else
    echo "Exiting: neither model.h5 nor model.tf exist."
    exit 1
fi

