#!/usr/bin/env python3
"""
Module with utility functions for managing experiments using the various models
"""
import os
from inspect import cleandoc

import dill
import tensorflow as tf

from lsf import bsub

def make_experiment_dir(model, path, data_func, overwrite=False, save_format='h5'):
    """
    Setup an experiment directory
    """
    path = path.rstrip('/')
    try:
        os.mkdir(path)
    except FileExistsError:
        if (not overwrite and
                any(os.path.exists(x) for x in (f'{path}/model.{save_format}',
                                                f'{path}/initial_model.{save_format}',
                                                f'{path}/model.yaml',
                                                f'{path}/training_log'))):
            raise FileExistsError(f'Model already exists in {path}')

    # Save model
    model.save(f'{path}/model.{save_format}', save_format=save_format)
    model.save(f'{path}/initial_model.{save_format}', save_format=save_format)

    # Dump data function
    with open(f'{path}/data_loader.dill', 'wb') as data_func_file:
        dill.dump(data_func, data_func_file, recurse=True)

# TODO better way of generating commands
def model_bsub(model_name, model_dir, big_job=False, ram=None,
               prefetch=400, shuffle_buffer=200, tensorboard=True,
               epochs=600, checkpoint=100, batch_size=100, validation_epochs=10,
               finetune=None, early_stop=None, gpu=1, save_format='h5', **kwargs):
    """
    Generate default bsub command string for training a model on the EBI farm
    """
    command = ['python bin/models/train_keras.py',
               f'-d {model_dir}/data_loader.dill',
               f'-l {model_dir}/training_log',
               f'-f {finetune}' if finetune is not None else '',
               f'-p {prefetch}',
               f'-s {shuffle_buffer}',
               f'-b {batch_size}',
               f'-e {epochs}',
               f'-v {validation_epochs}',
               '-t' if tensorboard else '',
               f'-c {checkpoint}' if checkpoint else '',
               f'-s {early_stop}' if early_stop else '',
               f'{model_dir}/model.{save_format}']
    command = [x for x in command if x]

    ram = ram if ram is not None else 30000 if big_job else 20000

    bsub_args = {'name': model_name, 'stdout': f'logs/{model_name}_train.%J',
                 'stderr': f'logs/{model_name}_train.%J.err', 'gpu': gpu,
                 'gpu_exclusive': big_job, 'ram': ram,
                 'hosts': "gpu-009 gpu-011", 'project': 'gpu',
                 'queue': 'research-rh74'}

    # kwargs can override bsub defaults
    bsub_args = {**bsub_args, **kwargs}

    return bsub(' '.join(command), **bsub_args)
