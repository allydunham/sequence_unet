#!/usr/bin/env python3
"""
Run SPBuild on all sequences in a fasta files
"""
import sys
import os
import multiprocessing
import subprocess
import itertools
import tempfile
import argparse
import pandas as pd
from Bio import SeqIO

def run_spbuild(seq, tempdir):
    """
    Run SPBuild
    """
    node_name = multiprocessing.current_process().name
    print(f'Processing {seq.id} on {node_name}', file=sys.stderr)

    # Write required Fasta
    fasta_path = f'{tempdir}/{node_name}.fa'
    SeqIO.write(seq, fasta_path, 'fasta')

    # Run SPBuild and cleanup
    mtx_path = f'{tempdir}/{node_name}.mtx'
    spbuild = subprocess.run(['spbuild', '-i', fasta_path, '-m', mtx_path])

    if not spbuild.returncode == 0:
        print(f'Error processing {seq.id}:', spbuild.std, sep='\n', file=sys.stderr)
        return None

    # Process Output
    mtx = pd.read_csv(mtx_path, skiprows=2, sep='\s+')
    mtx = mtx.reset_index().rename(columns={'level_0': 'position', 'level_1': 'wt'})
    mtx['protein'] = seq.id
    cols = mtx.columns.to_list()
    mtx = mtx[['protein'] + cols[:-1]]
    os.remove(fasta_path)
    os.remove(mtx_path)
    return mtx

def main(args):
    """
    Run SPBuild on input Fasta and format into a single output file
    """
    seqs = SeqIO.parse(args.fasta, 'fasta')
    with multiprocessing.Pool(processes=args.processes) as pool,\
         tempfile.TemporaryDirectory(dir=args.temp) as tempdir:
        profiles = pool.starmap(run_spbuild, zip(seqs, itertools.cycle([tempdir])))

    profiles = pd.concat(profiles, axis=0)
    profiles.to_csv(sys.stdout, sep='\t', index=False)

def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('fasta', metavar='F', help="Input Fasta")
    parser.add_argument('--processes', '-p', default=1, type=int,
                        help="Number of processes available")
    parser.add_argument('--temp', '-t', default='.', type=str,
                        help="Root location for tempory storage")

    return parser.parse_args()

if __name__ == "__main__":
    main(parse_args())
