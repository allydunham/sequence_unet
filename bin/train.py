#!/usr/bin/env python3
"""
Train a model using the inbuilt keras training functionality
"""
import argparse
import math
import logging
from pathlib import Path
from datetime import date

import dill
import tensorflow.keras as keras

from proteinnetpy import tfdataset
from metrics import CUSTOM_OBJECTS

def main(args):
    """Main script"""
    # import model
    root = Path(args.model).parent

    save_format = Path(args.model).suffix
    model = keras.models.load_model(args.model, custom_objects=CUSTOM_OBJECTS)

    # Set trainable if applicable (this is not loaded with the model due to a bug)
    if args.finetune is not None:
        model.trainable = True

        for layer in model.layers[:-args.finetune]:
            layer.trainable = False

    # Setup log
    if args.log:
        logging.basicConfig(filename=args.log, format='{asctime} {name} {message}',
                            style='{', level=logging.INFO)

    else:
        logging.basicConfig(format='{asctime} {name} {message}', style='{', level=logging.INFO)

    # Setup data
    if args.data:
        with open(args.data, 'rb') as dill_file:
            pn_map_func = dill.load(dill_file)
    else:
        pn_map_func = dill.load(f'{root}/data_loader.dill')
    train_map = pn_map_func(validation=False)
    val_map = pn_map_func(validation=True)
    train_len = len(train_map)
    val_len = len(val_map)

    train_data = tfdataset.proteinnet_tf_dataset(train_map, batch_size=args.batch_size,
                                                 prefetch=args.prefetch,
                                                 shuffle_buffer=args.shuffle_buffer)


    val_data = tfdataset.proteinnet_tf_dataset(val_map, batch_size=val_len,
                                               prefetch=2*val_len,
                                               shuffle_buffer=val_len)

    # setup checkpoints
    callbacks = []
    if args.tensorboard:
        callbacks.append(keras.callbacks.TensorBoard(log_dir=root,
                                                     histogram_freq=50,
                                                     write_graph=False,
                                                     profile_batch=0))

    if args.checkpoint:
        chkpt_path = f"{root}/chkpt_{date.today().isoformat()}_{{epoch:02d}}"
        samples_per_save = args.checkpoint * train_len
        callbacks.append(keras.callbacks.ModelCheckpoint(filepath=chkpt_path,
                                                         save_freq=samples_per_save))

    if args.early_stop != -1:
        callbacks.append(keras.callbacks.EarlyStopping(monitor="val_loss",
                                                       patience=args.early_stop,
                                                       restore_best_weights=True))

    # fit
    train_steps = math.ceil(train_len / args.batch_size)

    logging.info('Training for %s epochs; batch size: %s, shuffle buffer: %s, prefetch: %s',
                 args.epochs, args.batch_size, args.shuffle_buffer, args.prefetch)
    model.fit(train_data,
              epochs=args.epochs,
              steps_per_epoch=train_steps,
              validation_data=val_data,
              callbacks=callbacks,
              validation_steps=args.val_epochs)

    model.save(args.model, save_format=save_format)
    model.save(f'{root}/{date.today().isoformat()}_model{save_format}', save_format=save_format)

    metrics = model.evaluate(val_data, steps=args.val_epochs)
    if len(model.metrics_names) == 1:
        metrics = [metrics]

    metric_strings = [f'{m}: {v}' for m, v in zip(model.metrics_names, metrics)]
    logging.info('After training: %s', ', '.join(metric_strings))

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('model', metavar='M', help="Model path")

    parser.add_argument('--finetune', '-f', type=int, help="Number of layers to train")

    parser.add_argument('--data', '-d', help="Pickled data loading function")
    parser.add_argument('--log', '-l', help="Path to log file")

    parser.add_argument('--batch_size', '-b', default=50, type=int, help="Batch size")
    parser.add_argument('--prefetch', '-p', default=300, type=int, help="Prefetch size")
    parser.add_argument('--shuffle_buffer', '-s', default=200, type=int, help="Prefetch size")
    parser.add_argument('--epochs', '-e', default=100, type=int, help="Number of epochs")
    parser.add_argument('--val_epochs', '-v', default=5, type=int,
                        help="Number of passes through validation data per epoch")

    parser.add_argument('--tensorboard', '-t', action='store_true',
                        help="Track training using tensorboard")
    parser.add_argument('--checkpoint', '-c', default=0, type=int,
                        help="Checkpoint frequency in epochs")
    parser.add_argument('--early_stop', '-s', default=-1, type=int,
                        help="Enable early stopping on val_loss with a patience of N epochs")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
