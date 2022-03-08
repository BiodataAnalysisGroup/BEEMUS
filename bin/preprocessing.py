#!/usr/bin/env python
# coding: utf-8

###### Up to now command in order to run this file
#### python preprocessing.py --datafolder ../data --vcfs ../data/clinical_variant_files --lineagesfiles ../data/lineages --lineagesClasification "../data/SARS-CoV-2 lineage meta data.csv" --genes_coordinates_path ../data/ref/NC_045512.2_annot.gff3

import argparse
from os import listdir
from os.path import isfile, join
import os
import pathlib
import subprocess
from tqdm.notebook import trange
import re

import numpy as np
import pandas as pd

from src.vcfs_parser import parser as vcf_parser
from src.helpers import *


def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Retrieve outbreak lineages information.')
    parser.add_argument('--datafolder', type=pathlib.Path, required=True)
    parser.add_argument('--vcfs', type=pathlib.Path, required=True)
    parser.add_argument('--lineagesfiles', type=pathlib.Path, required=True)
    parser.add_argument('--lineagesClasification', type=pathlib.Path, required=True)
    parser.add_argument('--genes_coordinates_path', type=pathlib.Path, required=True)





    return parser

def genes_loading():
    orf1a_sub_names = ['NSP'+str(i) for i in range(1, 11)]
    orf1b_sub_names = ['NSP12a', 'NSP12b'] + ['NSP'+str(i) for i in range(13,17)]

    genes_coordinates = pd.read_csv(genes_coordinates_path, delimiter='\t', comment = '#', header = None)


    genes_coordinates = genes_coordinates.append({0: 'NC_045512.2', \
                                                1: 'Giorgi', \
                                                2: 'CDS', \
                                                3: genes_coordinates.loc[genes_coordinates[8].isin(orf1a_sub_names)][3].min(), \
                                                4: genes_coordinates.loc[genes_coordinates[8].isin(orf1a_sub_names)][4].max(), \
                                                5: '.', \
                                                6: '+', \
                                                7: '.', \
                                                8: 'ORF1a', \
                                                9: 'desc'} \
                                                ,ignore_index=True)

    genes_coordinates = genes_coordinates.append({0: 'NC_045512.2', \
                                                1: 'Giorgi', \
                                                2: 'CDS', \
                                                3: genes_coordinates.loc[genes_coordinates[8].isin(orf1b_sub_names)][3].min(), \
                                                4: genes_coordinates.loc[genes_coordinates[8].isin(orf1b_sub_names)][4].max(), \
                                                5: '.', \
                                                6: '+', \
                                                7: '.', \
                                                8: 'ORF1b', \
                                                9: 'desc'} \
                                                ,ignore_index=True)

    genes_coordinates.drop(labels = genes_coordinates[genes_coordinates[8].isin(orf1a_sub_names)].index, inplace = True)
    genes_coordinates.drop(labels = genes_coordinates[genes_coordinates[8].isin(orf1b_sub_names)].index, inplace = True)

    genes_coordinates.rename(columns={3: 'gene_start', 4: 'gene_end', 8: 'gene'}, inplace = True)
    genes_coordinates = genes_coordinates[['gene_start', 'gene_end', 'gene']].copy().set_index('gene')

    print(len(genes_coordinates), 'genes coordinates have been loaded.')

    genes_coordinates.to_csv(genes_coordinates_path.parent / (genes_coordinates_path.stem + "_simplified.csv"))
    return genes_coordinates



if __name__ == '__main__':
    
    parser = init_argparse()
    args = parser.parse_args()

    # Store arguments to variables
    data_folder = args.datafolder
    clinical_data_path = args.vcfs
    lineages_files = args.lineagesfiles
    lineages_clasification_path = args.lineagesClasification
    genes_coordinates_path = args.genes_coordinates_path

    # Create necessary directories
    if not os.path.exists(lineages_files):
        os.makedirs(lineages_files)

    tic('Loading genes...')
    genes_coordinates = genes_loading()
    toc('Genes loading completed.')
    
    tic('Loading lineages clasification file...')
    lineage_metadata = pd.read_csv(lineages_clasification_path, index_col = 'id', skiprows = [1,2])
    print(len(lineage_metadata), 'samples with known lineages have been loaded.')

    lineage_metadata['lineage'] = lineage_metadata['lineage'].apply(lambda x: re.sub(r" ?\([^)]+\)", "", str(x)))
    lineage_metadata['lineage'] = lineage_metadata['lineage'].str.strip()
    lineage_metadata['lineage'] = lineage_metadata['lineage'].apply(lambda x: np.nan if (x in ['nan', '']) else x)
    toc('Lineages clasification loading completed.')
    
    tic('Samples\' metadata loading...')
    metadata = pd.read_csv(data_folder / "biosample_result.csv", index_col = 'Title')
    print(len(metadata), 'metadata samples have been loaded.')
    toc('Samples\' metadata loading completed.')

    metadata = pd.concat([lineage_metadata, metadata], axis = 1)

    tic('Dowloading and merging lineages characteristics...')

    lineages = metadata['lineage'].unique().tolist()
    lineages = [i for i in lineages if pd.notna(i)]

    # comment out to download the files
    for idx, lineage in enumerate(lineages):
        cmd = 'scraper.py --o {} --lineage={}'.format(lineages_files, lineage)
        normal = subprocess.run(cmd,
            shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            check=True,
            text=True)

        print(normal.stdout)
        print(normal.stderr)
        print(lineage, ' have been downloaded')

    # comment out to download the files
    # for idx, lineage in enumerate(lineages):
    #     os.system('python scraper.py --o ../data/lineages_test --lineage={}'.format(lineage))
    # print(idx+1, 'files have been downloaded')
        
    lineages_files = [join(lineages_files, f) for f in listdir(lineages_files) if isfile(lineages_files / f) and \
            any(f.startswith(substring+'_') for substring in lineages)]

    lineages_data = pd.concat([pd.read_csv(f, usecols = ['lineage','gene', 'ref_aa', 'alt_aa', 'codon_num', 'codon_end']) for f in lineages_files])
    lineages_data['codon_end'].replace({"None": np.nan}, inplace=True)
    lineages_data['codon_end'] = pd.to_numeric(lineages_data['codon_end'])

    # lineages_data
    lineages_data = pd.merge(lineages_data, genes_coordinates, how='left', on = 'gene', validate = 'many_to_one')

    # Calculate mutation start-end coordinates
    lineages_data = lineages_data.assign(mut_start = lambda x: ((x['codon_num'] * 3) + x['gene_start'] - 3))

    lineages_data['mut_end'] = lineages_data.apply(lambda x: ((x['codon_num'] * 3) + x['gene_start'] - 1) \
                                                if (pd.isna(x['codon_end'])) \
                                                else ((float(x['codon_end']) * 3) + x['gene_start'] - 1), \
                                                axis=1).astype('int64')

    toc('Dowloading and merging lineages characteristics completed.')


    tic('Modeling...')
    p = vcf_parser(clinical_data_path, metadata, lineages_data)
    p.convert_to_bin()
    data = p.data

    result = pd.merge(data, metadata, how='left', left_index = True, right_index = True)
    msk = (result == 0).all() # get rid of columns containing only zeros
    # result = result.loc[:,~msk].copy()
    result.to_csv(data_folder / "dataset.csv")
    pd.DataFrame(result['lineage'].value_counts()).reset_index().rename(columns={"lineage": "counts", "index": "linage"}).to_csv(data_folder / "dataset_summary.csv", index=False)
    toc('Modeling completed.')