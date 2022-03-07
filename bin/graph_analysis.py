#!/usr/bin/env python
# coding: utf-8

###### Up to now command in order to run this file
#### python graph_analysis.py --datafolder ../data --datasetname dataset.csv --graphsfolder ../graphs --nodesfilename nodes_all --lineages B.1.1.7 B.1.617.2 AY.12 AY.9 AY.4 B.1.1.318 AY.7 --export_arrow_edges --export_nodes_long_format

import os
import subprocess
import argparse
import pathlib
import pandas as pd
import numpy as np
import networkx as nx
from pyvis.network import Network

from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN

import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import pickle
import json
import statistics
from scipy.stats import norm
from joblib import Parallel, delayed

from src.helpers import *
from src.vcfs_parser import parser
from src.counts import count_snps

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Graph analysis')
    # parser.add_argument('--lineage', type=str, help="The lineage of interest eg. 'B.1.617.2'")
    # parser.add_argument('--f', type=str, default="0.75", help="Characteristic mutations for a lineage are defined as nonsynonymous substitutions or deletions that occur of sequences within that lineage. For lineages with few sequences, the 75%% threshold may not identify all the mutations specific to that lineage, and as more sequences are found, the characteristic mutations may change. default=0.75")

    parser.add_argument('--datafolder', type=pathlib.Path, required=True)
    parser.add_argument('--datasetname', type=pathlib.Path, required=True)
    parser.add_argument('--graphsfolder', type=pathlib.Path, required=True)
    parser.add_argument('--nodesfilename', type=pathlib.Path, required=False, default=False)
    parser.add_argument('--lineages', nargs='+', help='<Required> One ore more', required=True)
    parser.add_argument('--export_arrow_edges', action='store_true', help="Export arrow edges to ... format.")
    parser.add_argument('--export_nodes_long_format', action='store_true', help="Export graph's nodes to a proper long format as a csv file.")
    parser.add_argument('--plot_dates_vs_nodes', action='store_true', help="Plot collection dates per node for each lineage.")

    return parser

def export_nodes_long():
    cluster_1 = []
    cluster_2 = []
    position = []
    f = []
    snps = []
    based = []
    node = []
    lineage = []
    mcl_group = []
    central = []
    bet_central = []
    l_central = []


    for bo in based_on:
        for i in nodes.index:
            line = nodes.loc[i, :]
            for key, value in line["based_on_{}_pos".format(bo)].items():
                
                position = extend_list(position, [int(k) for k in value])
                f = extend_list(f, line["based_on_{}_counts".format(bo)][key])
                snps = extend_list(snps, line["snps_on_{}".format(bo)][key])
                
                length = len(line["based_on_{}_pos".format(bo)][key])
                
                lineage = extend_list(lineage, [key] * length)
                based = extend_list(based, [bo] * length)
                
                node = extend_list(node, [i] * length)
                mcl_group = extend_list(mcl_group, [groups[i]] * length)
                central = extend_list(central, [centrality[i]] * length)
                bet_central = extend_list(bet_central, [betweenness_centrality[i]] * length)
                l_central = extend_list(l_central, [load_centrality[i]] * length)
                
                
                cluster_1 = extend_list(cluster_1, [line['cluster_1']] * length)
                cluster_2 = extend_list(cluster_2, [line['cluster_2']] * length)

    #         break
    #     break

    if check_if_lists_in_list_have_same_length([position, \
                                                f, \
                                                snps, \
                                                lineage, \
                                                based, \
                                                node, \
                                                mcl_group, \
                                                central, \
                                                bet_central, \
                                                l_central, \
                                                cluster_1, \
                                                cluster_2]):
        result_long = pd.DataFrame.from_records(zip(position, f, snps, lineage, based, cluster_1, cluster_2, node, mcl_group, central, bet_central, l_central), \
                                                columns = ['position', 'f', 'snps', 'lineage', 'based_on', 'cluster_1', 'cluster_2',\
                                                        'node', 'mcl_group', 'central', 'bet_central', 'l_central'])
        result_long.to_csv(data_folder / "nodes_long.csv", index=False)
    return 

def plot_collenction_dates_vs_nodes():
    x = []
    y = []
    lin = []


    cmap = matplotlib.cm.get_cmap('Set1', 9)
    colors = []
    for i in range(cmap.N):
        # rgb2hex accepts rgb or rgba
        colors.append(matplotlib.colors.rgb2hex(cmap(i)))
        
    cdict = {value : colors[i] for i,value in enumerate(lineages_of_interest)}

    for i in nodes.index:
        for j, val in nodes.loc[i, 'samples'].items():
            length = len(val)
            lin = extend_list(lin, [j] * length)
            x = extend_list(x, [i] * length)
            y = extend_list(y, df.loc[val, 'collection date'].tolist())

    y = [np.nan if d == 'not provided' else d for d in y]
    dates = pd.DataFrame.from_records(zip(y, x, lin), columns = ['dates', 'nodes', 'lineage'])
    dates.dropna(how = 'any')
    dates['dates'] =  pd.to_datetime(dates['dates'])


    fig, ax = plt.subplots()
    for lineage in dates['lineage'].unique():
        x = dates.loc[dates['lineage'] == lineage, 'nodes']
        y = dates.loc[dates['lineage'] == lineage, 'dates']
        
        ax.scatter(x, y, c=cdict[lineage], label=lineage)

    ax.set_ylabel('Collection date')
    ax.set_xlabel('Nodes')
    ax.legend()
    ax.grid(True)

    plt.savefig(graphs_folder / "dates_vs_nodes.pdf", format='pdf')
    plt.clf()
    
    return

if __name__ == '__main__':

    start_time = time.time()
    parser = init_argparse()
    args = parser.parse_args()

    # Init argument variables

    data_folder = args.datafolder
    dataset_name = args.datasetname
    graphs_folder = args.graphsfolder
    nodes_file_name = args.nodesfilename
    lineages_of_interest = args.lineages    # Lineages of interest
    export_arrow_edges = args.export_arrow_edges
    export_nodes_long_format = args.export_nodes_long_format
    plot_dates_vs_nodes = args.plot_dates_vs_nodes
    
    # Init programm constants
    based_on = [1, 2]

    # Create necessary directories
    if not os.path.exists(graphs_folder):
        os.makedirs(graphs_folder)

        if not os.path.exists(graphs_folder / "circular_images"):
            os.makedirs(graphs_folder / "circular_images")
    
    tic('Loading the dataset...')
    df = pd.read_csv(data_folder / dataset_name, index_col = 0, low_memory=False) # loading the dataset
    mut_cols = [col for col in df.columns.to_list() if col.isnumeric()] # retrive positional columns


    df = df.loc[df['lineage'].isin(lineages_of_interest)].copy()    # retrive the part of the dataset of the lineages of interest
    lineages = df.loc[:,'lineage'].copy()   # retrive the column with of the lineages
    toc('Datased loaded successfully.')


    tic('Clustering...')
    for bo in based_on:
        print('Based on {}'.format(bo))
        X = prepare_for_clustering(df, mut_cols, based_on = bo)

        X_transformed = TSNE(n_components=2,\
                        init='random', n_jobs=-1).fit_transform(X)
        X_transformed.shape

        # kmeans = KMeans(n_clusters=20, random_state=0).fit(X_transformed)
        # kmeans.labels_

        dbscan = DBSCAN(eps = 5, n_jobs = -1).fit(X_transformed)
        dbscan.labels_
        df['cluster_' + str(bo)] = dbscan.labels_
        n_clusters_ = len(set(dbscan.labels_)) - (1 if -1 in dbscan.labels_ else 0)
        n_noise_ = list(dbscan.labels_).count(-1)
        print('Estimated number of clusters', n_clusters_)
        print('Estimated number of noise', n_noise_)


        colors = [str(k) for k in dbscan.labels_.tolist()]
        centroids = np.zeros((len(set(dbscan.labels_)), 2))
        for i, label in enumerate(np.unique(dbscan.labels_)):
            msk = dbscan.labels_ == label
            points = [X_transformed[j,:] for j, value in enumerate(msk) if value]
            centroids[i] = np.mean(points, axis=0)




        palette = sns.color_palette("Paired") + sns.color_palette("Set3", 12)
        g = sns.scatterplot(x = X_transformed[:,0], y = X_transformed[:,1], hue = colors, style=lineages)
        g.scatter(x = centroids[:,0], y=centroids[:,1], color='black', marker='+', s=16, linewidths=0.7)
        g.set_title('TSNE->DBSCAN on the bin matrix ' + str(bo) + 's')
        plt.legend(bbox_to_anchor=(1,1), loc="upper left")
        plt.savefig(graphs_folder / str('TSNE->DBSCAN on the bin matrix ' + str(bo) + 's_final.pdf'), format='pdf', bbox_inches='tight')
        plt.clf()
    toc('Clustering completed successfully.')

    if not nodes_file_name:
        tic('Counting of snps... (this should take a while, be patient)')
        counter = count_snps(df, lineages_of_interest, based_on)
        counter.occurences_of_mutations()
        toc('Counting completed.')

    tic('Loading nodes...')

    if not nodes_file_name:
        # ------------- new method for all lineages -------------
        # Path is not provided and so calculate nodes

        counter.set_path(data_folder / "indices")
        # # ------------- new method for all lineages in parallel-------------
        def process(i):
            d_pos = {}
            d_counts = {}
            d_snps = {}
            for lineage in nodes.loc[i, 'samples'].keys():
                pos_col = "{lineage}_C{cluster}_pos_{bo}".format(lineage = lineage, cluster = nodes.loc[i, 'cluster_' + bo], bo = bo)
                count_col = "{lineage}_C{cluster}_counts_{bo}".format(lineage = lineage, cluster = nodes.loc[i, 'cluster_' + bo], bo = bo)


                d_pos[lineage] = remove_nan_from_list( \
                                        counter.counts.loc[:, pos_col].values.tolist())
                d_counts[lineage] = remove_nan_from_list( \
                                        counter.counts.loc[:, count_col].values.tolist())
                d_snps[lineage] = [counter.get_snps_from_lineages_clusters_position( \
                                            nodes.loc[i, 'samples'][lineage], \
                                            region = 'NC_045512.2', \
                                            position = int(k)) for k in d_pos[lineage]]
            return (d_pos, d_counts, d_snps)


        nodes = counter.clusters_intersection()
        # nodes.dropna(how = 'any', inplace = True)
        for bo in based_on:
            bo = str(bo)

            results = Parallel(n_jobs=-1, verbose = 11, backend = 'loky')(delayed(process)(i) for i in nodes.index)
            column_pos, column_counts, column_snps = zip(*results)

            nodes['based_on_' + bo + '_pos'] = column_pos
            nodes['based_on_' + bo + '_counts'] = column_counts
            nodes['snps_on_' + bo] = column_snps

        pickle.dump(nodes, open(data_folder / "nodes", "wb"))
    else:
        # Load the existing file
        nodes = pickle.load(open(data_folder / nodes_file_name, "rb"))
    toc('Loading nodes completed.')

    for i in nodes.index:
        counts = {}
        samples_dict = nodes.loc[i, 'samples']
        for key, value in samples_dict.items():
            increase_value_of_dict_by_value_new(counts, key, len(value))
        # Pie chart, where the slices will be ordered and plotted counter-clockwise:
        (labels, sizes) = zip(*counts.items())
        #     labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
        #     sizes = [15, 30, 45, 10]
        #     explode = (0, 0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')

        fig1, ax1 = plt.subplots()
        ax1.axis('tight')
        ax1.axis('off')
    #     fig1.set_size_inches(10, 10)
        ax1.pie(sizes, labels=labels, autopct='%1.1f%%',
                shadow=False, startangle=90, labeldistance = None, rotatelabels = True,
                wedgeprops={"edgecolor":"0",'linewidth': 0})
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

        plt.savefig(graphs_folder / 'circular_images/{}.png'.format(i), transparent=True, bbox_inches='tight', pad_inches = 0)
        plt.clf()
        plt.close()

    # If the corresponding argument is given, plot the the dates vs nodes graph
    if plot_dates_vs_nodes:
        plot_collenction_dates_vs_nodes()


    tic('Calculating fully connected graph edges')
    # ------------- new method for all lineages ---------------
    threshold = 1
    edges = []
    j_loop_indecies = nodes.index.tolist()
    # test = []
    for i in nodes.index:
        j_loop_indecies.remove(i)

        keys_i = nodes.loc[i, 'samples'].keys()
        
        for key_i in keys_i:
            pos_i_dict = list_to_index_val_dict(nodes.loc[i, 'based_on_1_pos'][key_i] + nodes.loc[i, 'based_on_2_pos'][key_i])
            counts_i = nodes.loc[i, 'based_on_1_counts'][key_i] + nodes.loc[i, 'based_on_2_counts'][key_i]
        #     print(pos_i_dict)
            for j in j_loop_indecies:
                for key_j in nodes.loc[j, 'samples']:
                    pos_j_dict = list_to_index_val_dict(nodes.loc[j, 'based_on_1_pos'][key_j] + nodes.loc[j, 'based_on_2_pos'][key_j])
                    counts_j = nodes.loc[j, 'based_on_1_counts'][key_j] + nodes.loc[j, 'based_on_2_counts'][key_j]
            #         print(pos_j_dict)
                    d = {}
                    intersection = list(set(pos_i_dict) & set(pos_j_dict))
                    for inter in intersection:
                        d[str(int(inter))] = min(get_list_items_from_idx_list(counts_i, pos_i_dict[inter]) + 
                                                get_list_items_from_idx_list(counts_j, pos_j_dict[inter]))

                    keys_to_be_deleted = []
                    for key, value in d.items():
                        if value < threshold:
                            keys_to_be_deleted.append(key)
                    for key in keys_to_be_deleted:
                        d.pop(key)
                    if d:
                        (pos, min_weight) = zip(*d.items())
                        length = len(pos)
                        _from = [key_i] * length
                        _to = [key_j] * length

                        if len(pos) == len(min_weight) == len(_from) == len(_to):
                            title = pd.DataFrame.from_records(zip(_from, _to, pos, min_weight), columns = ['from', 'to', 'position', 'weights']).to_html(index = False)
                            edges.append((i, j, {'weights' : d,
                                                'value' : len(d),
                                                'title' : title}))
                        else:
                            print('Something went wrong with the lengths')
            #         print(d)
            #         break
            #     break
    toc('Fully connected graph calculated successfully.')

    # Normalize and preprocess the graph's weights
    Q = 400
    cols = {
        'from' : [],
        'to' : [],
        'weight' : []
    }
    for edge in edges:
        cols['from'].append(edge[0])
        cols['to'].append(edge[1])
        cols['weight'].append(edge[2]['value'])
        
    x = cols['weight']
    _std = statistics.pstdev(cols['weight'])
    _mean = statistics.mean(cols['weight'])
    # _min = min(cols['weight'])
    # _max = max(cols['weight'])
    # cols['weight'] = [((w - _min) / (_max - _min)) * Q  for w in cols['weight']]
    cols['weight'] = [(w - _mean) / _std for w in cols['weight']]
    cols['weight'] = norm.cdf(cols['weight'], loc=0, scale=0.04) * Q
    pd.DataFrame.from_dict(cols).to_csv(data_folder / 'fully-connected-graph.csv', header = False, index = False, sep = " ")
    x = cols['weight']
    plt.scatter(range(len(x)), x)

    plt.title('Fully connected graph weights normalization.')
    plt.savefig(graphs_folder / "weights_normalization.pdf")
    plt.clf()

    # Do something to call mcl and perform calculations on
    # the fully connected graph
    normal = subprocess.run('mcl fully-connected-graph.csv --abc -I 30 -scheme 7 -o fully-connected-graph-clusters.out',
        cwd = data_folder,
        shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=True,
        text=True)

    print(normal.stdout)
    print(normal.stderr)



    tic('Calculating arrow edges...')
    # -----------------new arrow adges calculation for all lineages-----------------
    edges = []
    for bo in based_on:
        fo = str(based_on[~based_on.index(bo)])
        bo = str(bo)
        for i in nodes.index:


            pos_i = flatten_list_of_lists(list(nodes.loc[i, 'based_on_' + bo + '_pos'].values()))

            same_clusters = nodes.loc[(nodes['cluster_' + bo] == nodes.loc[i, 'cluster_' + bo]) & \
                                    (nodes['cluster_' + fo] != nodes.loc[i, 'cluster_' + fo]), :]
            # if bo == '2' and i == 10:
            #     print(same_clusters)
            if not same_clusters.empty:
        #             print(same_clusters[['cluster_1', 'cluster_2']])
                for j in same_clusters.index:
                    d = {}
                    i_snps = flatten_list_of_lists(list(nodes.loc[i, 'snps_on_'+ bo].values()))
                    for l, snp_i in enumerate(i_snps):
                        if snp_i[0] == 0:
                            pos_j = flatten_list_of_lists(list(same_clusters.loc[j, 'based_on_' + bo + '_pos'].values()))
                            j_snps = flatten_list_of_lists(list(same_clusters.loc[j, 'snps_on_'+ bo].values()))
                            
                            for k, snp_j in enumerate(j_snps):
                                if (snp_j[0] > 0) and (pos_j[k] == pos_i[l]):
                                    # if bo == '2' and i == 10 and j == 9:
                                    #     print('on {bo} empty pos {pos_i} of node {i} -> {pos_j} of node {j}'.format(bo = bo, pos_i = pos_i[l], pos_j = pos_j[k], i = i, j = j))
                                    d[pos_i[l]] = (snp_i, snp_j)
                if d:
                    title = pd.DataFrame.from_dict(d, orient = 'index', columns = ['from', 'to']) \
                                .reset_index() \
                                .rename(columns = {'index':'position'}) \
                                .to_html(index = False)
                    edges.append((i, j, {'color' : 'blue' if bo == "1" else 'red',
                                        'arrows' : 'to',
                                        'data' : d, 
                                        'title' : title
                                        }))

    if export_arrow_edges:
        #DO SOMETHING HERE TO EXPORT ARROW EDGES TO A PROPER FORMAT
        _from = []
        _to = []
        _type = []

        for edge in edges:
            _from.append(edge[0])
            _to.append(edge[1])
            _type.append(edge[2]['color'])
            
        pd.DataFrame(np.column_stack([_from, _to, _type]), 
                                       columns=['from', 'to', 'color']).to_csv(data_folder / "arrow_edges.csv", index=False)

        ## Different csv format --comment out to use it--
        # arrow_edges = {}
        # for edge in edges:
        #     if edge[0] not in arrow_edges:
        #         arrow_edges[edge[0]] = [(edge[1], edge[2]['color'])]
        #     else:
        #         if (edge[1], edge[2]['color']) not in arrow_edges[edge[0]]:
        #             arrow_edges[edge[0]].append((edge[1], edge[2]['color']))

        # pd.DataFrame.from_dict(arrow_edges, orient='index').to_csv(data_folder / "arrow_edges.csv")

    # Either way save the edges to a binary file
    pickle.dump(edges, open(data_folder / "arrow_edges", "wb"))
    toc('Arrow edges calculated successfylly.')

    tic('Graph rendering...')
    #-------------------new method for all lineages--------------
    d = pd.read_csv(data_folder / 'fully-connected-graph-clusters.out', sep='\t', header = None).T.to_dict(orient='list')
    d = remove_nan_from_dict_of_lists(d)
    d = int_dict_of_lists(d)
    groups = {}
    for key, value in d.items():
        for element in value:
            groups[element] = key

    nodes['id'] = nodes.index
            
    G = nx.DiGraph(name = lineages_of_interest)
    G.add_nodes_from(list(nodes.to_dict(orient='index').items()))
    G.add_edges_from(edges)

    centrality = nx.algorithms.centrality.degree_centrality(G)
    betweenness_centrality = nx.algorithms.centrality.betweenness_centrality(G)
    load_centrality =  nx.algorithms.centrality.load_centrality(G)

    if export_nodes_long_format:
        export_nodes_long()

    net = Network(height='750px', width='100%', bgcolor='#222222', font_color='white', directed = True, heading = " ".join(str(x) for x in G.name))
    net.set_options(json.dumps({
        "configure": {
            "enabled": False
        },
        "nodes":{
            "borderWidth": 10,
            "borderWidthSelected": 10,
            "size": 400,
            "chosen": True,
            "font": {
                "size": 11,
            }
        },
        "edges": {
            "scaling": {
                "min": 0.1,
                "max": 1,
            },
            "selectionWidth" : 0,
            "color": {
                "inherit": False
            },
            "smooth": {
                "enabled": False,
                "type": "continuous"
            }
        },
        "interaction": {
            "dragNodes": True,
            "hideEdgesOnDrag": False,
            "hideNodesOnDrag": False
        },
        "physics": {
            "enabled": False,
        },
        "layout": {
            "randomSeed": 100,
        }
    }))

    net.from_nx(G)
    for i in range(len(net.nodes)):
        pos = []
        snps = []
        index = []
        total = []
        lineage = []
        f = []
    #     total_len = 0
        for bo in based_on:
            for lin in net.nodes[i]['based_on_{}_pos'.format(bo)].keys():
                bo = str(bo)

                pos = extend_list(pos, [int(k) for k in net.nodes[i]['based_on_{}_pos'.format(bo)][lin]])
                f = extend_list(f, ["{:.2f}".format(k) for k in net.nodes[i]['based_on_{}_counts'.format(bo)][lin]])
                
                length = len(net.nodes[i]['based_on_{}_pos'.format(bo)][lin])
                index = extend_list(index, [bo] * length)
                lineage = extend_list(lineage, [lin] * length)


                total = extend_list(total, [k[0] for k in net.nodes[i]['snps_on_{}'.format(bo)][lin]])
                snps = extend_list(snps, [str(k[1]) for k in net.nodes[i]['snps_on_{}'.format(bo)][lin]])
                
        if len(index) == len(lineage) == len(pos) == len(f) == len(total) == len(snps):
            net.nodes[i]['title'] = pd.DataFrame.from_records(zip(index, lineage, pos, f, total, snps), columns = ['Based_on', 'Lineage', 'Position', 'f', 'Total', 'Snps']).to_html(index = False) 
        #     total_len = total_len + length 
            net.nodes[i]['label'] = " ".join(str(x) for x in list(set(lineage))) + '\n' + str(net.nodes[i]['id'])
            net.nodes[i]['group'] = groups[net.nodes[i]['id']]
            net.nodes[i]['shape'] = "circularImage"
            net.nodes[i]['image'] = str("circular_images/{}.png".format(net.nodes[i]['id']))
    #         net.nodes[i]['color'] = 'black'
        else:
            print('Something went wrong with the lengths')
    tmp = "-".join(lineages_of_interest)
    net.show(str(graphs_folder / "graph-{}.html".format(tmp)))
    toc('Graph rendering completed.')

