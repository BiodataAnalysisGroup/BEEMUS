#!/usr/bin/env python
# coding: utf-8

import json
import urllib
import urllib.request
import pandas as pd
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

    # request url
    urlreq = "https://api.outbreak.info/genomics/lineage-mutations?pangolin_lineage={lineage}&frequency={freq}".format(lineage = lineage, freq = freq)

    # get response
    response = urllib.request.urlopen(urlreq)
    # load as json
    jresponse = json.load(response)

    if jresponse['success'] and jresponse['results'].keys():
        data = pd.DataFrame.from_dict(jresponse['results'][lineage])
        data['amino acid'] = data['mutation'].str.split(':',1).str[1]
        file_name = "{lineage}_f{freq}.csv".format(lineage = lineage, freq = freq)
        data.to_csv(file_name, index = False)
        print("Data has been saved to the {} file in the current directory.".format(file_name))
    else:
        print('Something went wrong!')
