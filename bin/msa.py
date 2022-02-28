#!/usr/bin/env python
# coding: utf-8

###### Up to now command in order to run this file
#### python msa.py --datafolder ../data --datasetname dataset.csv --msafilename combined-MSA-keep_length_ep0.1.fasta --graphsfolder ../graphs --nodesfilename nodes_all --treefilepath ../../tree-of-life/files/life.json --lineages B.1.1.7 B.1.617.2 AY.12 AY.9 AY.4 B.1.1.318 AY.7 --nodes_of_interest 23 42 32 14 15 3 16 25 24 2 4 5 33 38 35 61


import sys
# sys.path.append('../')

import os
import argparse
import pathlib
import pandas as pd

from sklearn.cluster import AgglomerativeClustering

import matplotlib.pyplot as plt
import pickle
import itertools
import json

from scipy.stats import norm


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from src.helpers import *
from src.vcfs_parser import parser
from src.msa_analysis import msa_parser


def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Graph analysis')
    # parser.add_argument('--lineage', type=str, help="The lineage of interest eg. 'B.1.617.2'")
    # parser.add_argument('--f', type=str, default="0.75", help="Characteristic mutations for a lineage are defined as nonsynonymous substitutions or deletions that occur of sequences within that lineage. For lineages with few sequences, the 75%% threshold may not identify all the mutations specific to that lineage, and as more sequences are found, the characteristic mutations may change. default=0.75")

    parser.add_argument('--datafolder', type=pathlib.Path, required=True)
    parser.add_argument('--datasetname', type=pathlib.Path, required=True)
    parser.add_argument('--msafilename', type=pathlib.Path, required=True)
    parser.add_argument('--graphsfolder', type=pathlib.Path, required=True)
    parser.add_argument('--nodesfilename', type=pathlib.Path, required=False, default=False)
    parser.add_argument('--treefilepath', type=pathlib.Path, required=True)
    parser.add_argument('--nodes_of_interest', nargs='+', help='<Required> One ore more', required=True)
    parser.add_argument('--lineages', nargs='+', help='<Required> One ore more', required=True)


    return parser

if __name__ == '__main__':

    start_time = time.time()
    parser = init_argparse()
    args = parser.parse_args()

    # Init argument variables

    data_folder = args.datafolder
    dataset_name = args.datasetname
    msa_name = args.msafilename
    graphs_folder = args.graphsfolder
    nodes_file_name = args.nodesfilename
    tree_file_path = args.treefilepath
    nodes_of_interest = args.nodes_of_interest    # Nodes of interest
    lineages_of_interest = args.lineages    # Lineages of interest
    
    nodes_of_interest = [int(k) for k in nodes_of_interest]
    # print(lineages_of_interest)
    # Init programm constants

    # Create necessary directories
    if not os.path.exists(graphs_folder):
        os.makedirs(graphs_folder)

    if not os.path.exists(graphs_folder / "mutual_info-entropy"):
        os.makedirs(graphs_folder / "mutual_info-entropy")

    tic('Loading the MSA, graph, edges and dataset...')
    # Load the MSA
    names = []
    seqs = []
    lengths = []
    uniques = []
    for seq_record in SeqIO.parse(data_folder / msa_name, "fasta"):
        names.append(seq_record.id)
        seqs.append(str(seq_record.seq))
        uniques.append(list(set(str(seq_record.seq))))
        lengths.append(len(seq_record))

    # Load the existing nodes file
    nodes = pickle.load(open(data_folder / nodes_file_name, "rb"))
    
    # Load the datset
    df = pd.read_csv(data_folder / dataset_name, index_col = 0, low_memory=False)
    
    # Load the arrow edges
    edges = pickle.load(open(data_folder / 'arrow_edges', "rb"))

    # edges = pd.read_csv(data_folder / "arrow_edges.csv")
    # edges = list(edges.itertuples(index=False, name=None))

    toc('Loading completed.')


    permutations_of_nodes_of_interest = list(itertools.permutations(nodes_of_interest, 2))
    tmp = []
    for i in nodes.loc[nodes_of_interest, 'samples']:
        tmp.append(list(i.values()))

    samples = list(set(flatten_list_of_lists(flatten_list_of_lists(tmp))))
    prefix_samples = add_prefix_to_list_of_strings(samples, 'New|1_')
    msa_indices = flatten_list_of_lists(search_elements_in_list(names, prefix_samples))

    sequences = get_list_items_from_idx_list(seqs, msa_indices)
    sequences_names = get_list_items_from_idx_list(names, msa_indices)

    record_list = []
    for (name, seq_string) in zip(sequences_names, sequences):
        record = SeqRecord(Seq(seq_string), id=name, description="")
        record_list.append(record)

    SeqIO.write(record_list, data_folder / "intermidiate_msa.fasta", "fasta")
    msa = msa_parser(data_folder / "intermidiate_msa.fasta")

    #-----------------------------------
    tic('Calculation of entropies and mutual info gain...')
    msa.calculate_probabilities()
    indices_of_pos_with_length_of_probs_larger_than_1 = [idx for idx, value in enumerate(msa.probabilities) if len(value) > 1]

    # retrive edges mutations of the nodes of interest
    edges_of_interest = search_in_edges(edges, permutations_of_nodes_of_interest)
    edges_data = list(set([int(i) for i in flatten_list_of_lists([list(k[2]['data'].keys()) for k in edges_of_interest])]))

    msa.probabilities = get_list_items_from_idx_list(msa.probabilities, edges_data)
    ### ara theloume mono tis theses poy de;ixnoyn ta arrows na ta psaxnoume sto msa kai na ypologizoume tis entropies giati
    ### giati endiaferomaste mono gia ekeines tis theses pou aforoun tha katallhla zeygaria clusters kai OXi gia oles tis theseis
    ### pou fainontai metalagmenes sto MSA metaksy twn samples of interest


    # ## sort the probabilities based on the positions (edge_data) in an asceding order
    # sorted_probabilities = [x for _, x in sorted(zip(edges_data, msa.probabilities), key=lambda pair: pair[0])]
    # sorted_edge_data = sorted(edges_data)

    # ([(ind,prob) for ind,prob in zip(sorted(edges_data), sorted_probabilities)])
    msa.calculate_entropies()

    ## sort the entropies based on the positions (edge_data) in an asceding order
    sorted_entropies = [x for _, x in sorted(zip(edges_data, msa.entropies), key=lambda pair: pair[0])]
    sorted_edge_data = sorted(edges_data)
    sorted_edge_data_frequencies = {}

    for i in sorted_edge_data:
        d = df.loc[samples, str(i)].value_counts(normalize = True).to_dict()
        dd = {}
        for key, value in d.items():
            if key == 2:
                dd['chr'] = value * 100
            elif key == 1:
                dd['non-chr'] = value * 100
                
        sorted_edge_data_frequencies[i] = dd

    plt.bar([str(k) for k in sorted_edge_data], sorted_entropies)
    plt.xticks(rotation=90, fontsize=6)
    plt.savefig(graphs_folder / 'mutual_info-entropy/entropies-nodes_{}.pdf'.format('-'.join([str(k) for k in nodes_of_interest])))
    plt.clf()
    plt.close() 

    msa.calculate_mutual_info(array = msa.transposed_array_of_sequences[sorted_edge_data, :])
    toc('Calculation of entropies and mutual info gain done.')

    tic('Compute heatmap...')
    fig, ax = plt.subplots()
    # ax.set_title(' | '.join([k.split('_')[1] for k in sequences_names]))

    im, cbar = heatmap(msa.mutual_info, sorted_edge_data, sorted_edge_data, ax=ax,
                    cmap="YlGn", cbarlabel="mutual info score", interpolation='none')
    # texts = annotate_heatmap(im, valfmt=matplotlib.ticker.FuncFormatter(func), size=6)

    fig.tight_layout()
    plt.savefig(graphs_folder / 'mutual_info-entropy/mutual_info-nodes_{}.pdf'.format('-'.join([str(k) for k in nodes_of_interest])))
    plt.clf()
    plt.close()
    toc('Heatmap computation completed.')


    tic('Agglomerative clustering...')
    similarity = 1 - msa.mutual_info

    # setting distance_threshold=0 ensures we compute the full tree.
    model = AgglomerativeClustering(distance_threshold=0, n_clusters=None, affinity='precomputed', linkage='complete')

    model = model.fit(similarity)
    # plt.title("Hierarchical Clustering Dendrogram")
    # plot the top three levels of the dendrogram
    plot_dendrogram(model, truncate_mode="level", p=100)#, labels=sorted_edge_data)
    # plt.xlabel("Number of points in node (or index of point if no parenthesis).")
    plt.savefig(graphs_folder / 'mutual_info-entropy/h-clustering-nodes_{}.pdf'.format('-'.join([str(k) for k in nodes_of_interest])))
    plt.clf()
    plt.close()
    toc('Agglomerative clustering completed')


    tic('Tree rendering...')
    genes_coordinates = pd.read_csv(data_folder / 'ref/NC_045512.2_annot_simplified.csv').sort_values(by='gene_start')

    # gene_colors = [(1, 100, 50),
    #                (38, 100, 50),
    #                (58, 100, 50),
    #                (85, 100, 50),
    #                (143, 100, 50),
    #                (177, 100, 50),
    #                (202, 100, 50),
    #                (236, 100, 50),
    #                (263, 100, 50),
    #                (285, 100, 50),
    #                (316, 100, 50),
    #                (1, 100, 70)]
    # gene_faded_color = [(gene_color[0], 30, gene_color[2]) for gene_color in gene_colors]

    gene_colors = [
                    '#BAD1E9',
                    '#E8CBB7',
                    '#7BC146',
                    '#43C3D9',
                    '#3678BE',
                    '#9AA5D1',
                    '#6260AB',
                    '#A67CB8',
                    '#ED5588',
                    '#FDD540',
                    '#F4842F',
                    '#A8CF58',
    ]

    genes_coordinates['gene_color'] = gene_colors
    # genes_coordinates['gene_color_faded'] = gene_faded_color

    wierd_chars = ['†', 'ψ', '⛥', '✠', 'ɸ', '♣', '▲']
    lineages_to_chars = {k : wierd_chars[i] for i, k in enumerate(lineages_of_interest)}
    sorted_edge_metadata = []

    for i in sorted_edge_data:   
        lins = df.loc[df[str(i)] == 2, 'lineage'].unique().tolist()
        if lins:
            tmp = []
            for j in lins:
                try:
                    tmp.append(lineages_to_chars[j])
                except:
                    pass
            sorted_edge_metadata.append(" ".join(tmp))
        else:
            sorted_edge_metadata.append("")

    with open(tree_file_path, "w") as write_file:
        labels = []
        for i, k in enumerate(sorted_edge_data):
            chr_f = sorted_edge_data_frequencies[k].get('chr', '')
            non_chr_f = sorted_edge_data_frequencies[k].get('non-chr', '')
            
            pos_color = coordinate_belongs_to_gene_color(k, genes_coordinates, fade = False)
            non_color = 'blue'
            if non_chr_f != '' and chr_f == '' and non_chr_f <= 5:
                pos_color = '#BABABA'
                non_color = '#BFC7FF'
            
            labels.append(
                {'pos' : {'value' : str(k), 'color' : pos_color},
                'symbols' : {'value' : sorted_edge_metadata[i], 'color' : 'green'},
                'freqs' : {'chr' : {'value' : check_if_frequency_is_empty(chr_f), 'color' : 'red'},
                            'non_chr' : {'value' : check_if_frequency_is_empty(non_chr_f), 'color' : non_color}
                        }
                }
            )
        json.dump(model_tree_to_nested_dics(model, labels), write_file, indent = 4)

    toc('Tree rendering completed.')