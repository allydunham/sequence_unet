#!/usr/bin/env python3
"""
Extract data from TensorBoard logs to
"""
import argparse
import sys
from os import path

from tensorboard.backend.event_processing.event_accumulator import EventAccumulator

def main(args):
    """
    Main function
    """
    print("model", "metric", "step", "value", sep="\t", file=sys.stdout)
    for log_path in args.events:
        model = path.dirname(log_path)
        acc = EventAccumulator(log_path).Reload()
        for tag in acc.Tags()['scalars']:
            for event in acc.Scalars(tag):
                print(model, tag, event.step, event.value, sep="\t", file=sys.stdout)

def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('events', metavar="E", nargs="+", help="TensorBoard event logs")

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())
