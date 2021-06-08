#!/usr/bin/env bash
# Utility script to run common commands on models
# Path is expected to be:
# reset: model directory(s)
# init: init script(s)
# train: model directory(s)
# check: model directory(s)

command=$1
paths=( "${@:2}" )

check_dir () {
    if ! test -f "$1/train.sh"; then
        echo "$1 does not look like a model dir"
        exit 1
    fi
}

check_init () {
    if ! [[ "$1" == *.py ]]; then
        echo "$1 does not look like an init script"
        exit 1
    fi
}

# Reset model dir(s)
if [ "$command" = "reset" ]; then
    for path in "${paths[@]}"; do
        check_dir "$path"

        rm "$path"/train/* "$path"/validation/*
        rm "$path"/training_log.*
        if test -d "$path/model.tf"; then
            rm -r "$path/model.tf"
            cp -r "$path/initial_model.tf" "$path/model.tf"
        elif test -f "$path/model.h5"; then
            rm "$path/model.h5"
            cp "$path/initial_model.h5" "$path/model.h5"
        else
            echo "Exiting: neither model.h5 nor model.tf exist."
            exit 1
        fi
    done

# Init bsub
elif [ "$command" = "init" ]; then
    for path in "${paths[@]}"; do
        check_init "$path"
        name=${path##*/}
        name=${name%.py}
        bsub -J "${name}_init" -o "logs/${name}_init.%J" -e "logs/${name}_init.%J.err" -M 4000 -R "rusage[mem=4000]" -P "gpu" -q "research-rh74" -m "gpu-009 gpu-011" -gpu - "python $path"
    done

# Run train scripts
elif [ "$command" = "train" ]; then
    for path in "${paths[@]}"; do
        check_dir "$path"
        bash "$path"/train.sh
    done

# Check training log(s)
elif [ "$command" = "check" ]; then
    for path in "${paths[@]}"; do
        check_dir "$path"
        grep -E "Successfully completed|Exited with" "${path}"/training_log.*
    done
fi

