#!/usr/bin/env python
# coding: utf-8

import argparse


def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Retrieve outbreak lineages information.')
    parser.add_argument('--lineage', type=str, help="The lineage of interest eg. 'B.1.617.2'")
    parser.add_argument('--f', type=str, default="0.75", help="Characteristic mutations for a lineage are defined as nonsynonymous substitutions or deletions that occur of sequences within that lineage. For lineages with few sequences, the 75%% threshold may not identify all the mutations specific to that lineage, and as more sequences are found, the characteristic mutations may change. default=0.75")
    return parser

if __name__ == '__main__':
    
    parser = init_argparse()
    args = parser.parse_args()

    lineage = args.lineage
    freq = args.f