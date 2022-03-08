#!/usr/bin/env python
# coding: utf-8

import json
import pathlib
import urllib
import urllib.request
import pandas as pd
import argparse

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Retrieve outbreak lineages information.')
    parser.add_argument('--lineage', type=str, help="The lineage of interest eg. 'B.1.617.2'")
    parser.add_argument('--f', type=str, default="0.75", help="Characteristic mutations for a lineage are defined as nonsynonymous substitutions or deletions that occur of sequences within that lineage. For lineages with few sequences, the 75%% threshold may not identify all the mutations specific to that lineage, and as more sequences are found, the characteristic mutations may change. default=0.75")
    parser.add_argument('--o', type=pathlib.Path, required=True, help="Specify the output directory of the program.")
    
    return parser

if __name__ == '__main__':
    
    parser = init_argparse()
    args = parser.parse_args()

    lineage = args.lineage
    freq = args.f
    output_path = args.o

    # request url
    urlreq = "https://api.outbreak.info/genomics/lineage-mutations?pangolin_lineage={lineage}&frequency={freq}".format(lineage = lineage, freq = freq)

    req = urllib.request.Request(urlreq)
    req.add_header('authorization', 'Bearer 0ed52bbfb6c79d1fd8e9c6f267f9b6311c885a4c4c6f037d6ab7b3a40d586ad0')
    # get response
    response = urllib.request.urlopen(req)
    # load as json
    jresponse = json.load(response)

    if jresponse['success'] and jresponse['results'].keys():
        data = pd.DataFrame.from_dict(jresponse['results'][lineage])
        data['amino acid'] = data['mutation'].str.split(':',1).str[1]
        file_name = output_path / "{lineage}_f{freq}.csv".format(lineage = lineage, freq = freq)
        data.to_csv(file_name, index = False)
        print("Data has been saved to the {} file in the {} directory.".format(file_name, output_path))
    else:
        print('Something went wrong!')
