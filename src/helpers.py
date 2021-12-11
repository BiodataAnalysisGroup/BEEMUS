from IPython.core.display import display, HTML
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from wordcloud import WordCloud
import math


def force_show_all(df):
    with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', None):
        display(HTML(df.to_html()))

def negate_list(lst):
    return [not elem for elem in lst]

def bar_plot(data):
    # set width of bars
    barWidth = 0.25

    # separate and derive the columns to which contains the positions and the counts
    pos_columns = [k for k in data.columns.to_list() if k.startswith('mut')]
    count_columns = [k for k in data.columns.to_list() if k.startswith('count')]
    
    
    # data[count_columns] = data[count_columns].apply(lambda x: [0 if y < 0.9 else y for y in x])
    # data.dropna(how='all', subset = count_columns, inplace = True)

    
    
    fig = plt.figure(figsize=(20,10))
    for i, (pos_col, count_col) in enumerate(zip(pos_columns, count_columns)):
        tmp = data[[pos_col, count_col]].dropna(how='all')
        d = pd.Series(tmp[count_col].values, index = tmp[pos_col]).to_dict()
    
        wordcloud = WordCloud(width=1600, height=800, max_font_size=100, max_words=1000, background_color="white")
        wordcloud.generate_from_frequencies(frequencies=d)

        ax = fig.add_subplot(math.ceil(len(pos_columns) / 3),3,i+1)

        ax.imshow(wordcloud, interpolation="bilinear")
        ax.axis('off')
    plt.title('xa')
    plt.tight_layout()
    plt.savefig('cloud.pdf')

def count_occurences_of_mutations(data, lineages, based_on, samples_threshold = 5, n_clusters_to_proceed = 0, m_lineages_to_proceed = 0, normalize = False):
    # Leave n_clusters_to_proceed, m_lineages_to_proceed to their default values = 0 to proceedd them all
    # Identifiy which columns coresponds to mutations postions and not metadata (the column name must be numeric)
    mutations_columns = [col for col in data.columns.to_list() if col.isnumeric()]

    # Init some variables to store results throughout the loops
    counts = {}
    lengths = {
        'lineage' : [], \
        'cluster' : [], \
        'based_on' : [], \
        'total' : [],
        }
    for bo in based_on:
        print('Based on {}'.format(bo))
        # Retrive all clusters in the dataset
        clusters = data['cluster_' + str(bo)].unique()

        for m, lineage in enumerate(lineages):
            print("Lineage {} is now being processed...".format(lineage))
            for n, cluster in enumerate(clusters):
                # For each lineage and cluster combination split the data to the mutations/positions
                # and their metadata
                sub_frame = data.loc[(data['cluster_' + str(bo)] == cluster) & (data['lineage'] == lineage), mutations_columns]
                sub_frame_metadata = data.loc[(data['cluster_' + str(bo)] == cluster) & (data['lineage'] == lineage), ~data.columns.isin(mutations_columns)]

                # If the sub_frame (which coresponds to the samples of the given lineage and cluster) has more than
                # "samples_threshold" samples, continue, else skip this lineage, cluster combination 
                l = len(sub_frame)
                if l > samples_threshold:
                    # Keap track of the lengths in a dictionary 
                    lengths['lineage'].append(lineage)
                    lengths['cluster'].append(cluster)
                    lengths['based_on'].append(bo)
                    lengths['total'].append(l)
                    # lengths["{}_C{}_{}".format(lineage, cluster, bo)] = l

                    # Count occurences of the levels in the mutations "sub_frame" and sort to a descending order  
                    occurrences = sub_frame.apply(pd.value_counts, normalize=normalize).sort_values(by = bo, axis = 1, ascending = False)
                    # Take into account only the 1s
                    occurrences = occurrences.loc[bo, ~occurrences.loc[bo,:].isnull()]
                    # Store positions and counts to a dictionary
                    counts["{}_C{}_pos_{}".format(lineage, cluster, bo)] = [float(k) for k in occurrences.index]
                    counts["{}_C{}_counts_{}".format(lineage, cluster, bo)] = occurrences.values

                # for the inner loop
                if n == n_clusters_to_proceed - 1:
                    break
            # for the outer loop
            if m == m_lineages_to_proceed - 1:
                    break

    # Store results to dataFrames and return
    summary = pd.DataFrame.from_dict(lengths)
    result_counts = pd.DataFrame(dict([ (k, pd.Series(v)) for k, v in counts.items() ]))
    return result_counts, summary

def prepare_for_clustering(df, cols_of_interest, based_on):
    mutations = df.loc[:, cols_of_interest]
    samples = mutations.index
    positions = mutations.columns
    X = mutations.to_numpy()
    # X = np.transpose(mutations.to_numpy())
    print('Samples length =', len(samples))
    print('Features length =', len(positions))

    if based_on == 1:
        X[X == 2] = 0

    elif based_on == 2:
        X[X == 1] = 0
        X[X == 2] = 1

    return X

def base_0_to_1(x):
    return x + 1

def base_1_to_0(x):
    return x - 1

def increase_value_of_dict_by_value(d, key, value):
        if key in d:
            d[key] = d[key] + value
        else:
            d[key] = 1
        return d

def increase_value_of_dict_by_value_new(d, key, value):
        if key in d:
            d[key] = d[key] + value
        else:
            d[key] = value
        return d

def remove_nan_from_list_of_lists(list):
    try:
        return [[z_ij for z_ij in z_i if not pd.isna(z_ij)] for z_i in list]
    except TypeError:
        return list

def remove_nan_from_dict_of_lists(dict):
    try:
        for key, value in dict.items():
            dict[key] = [element for element in value if not pd.isna(element)]
        return dict
    except TypeError:
        return dict

def int_dict_of_lists(dict):
    try:
        for key, value in dict.items():
            dict[key] = [int(element) for element in value]
        return dict
    except TypeError:
        return dict

def remove_nan_from_list(list):
    try:
        return [z for z in list if not pd.isna(z)]
    except TypeError:
        return list

def extend_list(list1, data):
    try:
        list1.extend(data)
    except TypeError:
        list1.append(data)
    return list1
    
def flatten_list_of_lists(list):
    try:
        return [z for sublist in list for z in sublist]
    except TypeError:
        return list 

def list_to_index_val_dict(lst):
    msk = [1] * len(lst)
    d = {}
    for i, val in enumerate(lst):
        if msk[i]:
            indices = get_indices(val, lst)
            d[val] = indices
            for idx in indices:
                msk[idx] = 0
    return d

def get_indices(x, list):
    return [i for i, value in enumerate(list) if value == x]

def get_list_items_from_idx_list(list, indices):
    return [list[i] for i in indices]

def check_if_lists_in_list_have_same_length(list_of_lists):
    return all(len(i) == len(list_of_lists[0]) for i in list_of_lists)
