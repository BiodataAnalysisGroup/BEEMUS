from IPython.core.display import display, HTML
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from wordcloud import WordCloud
import math
from scipy.cluster.hierarchy import dendrogram
import itertools
import time


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
                    # Take into account only the bos
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

def search_elements_in_list(list, elements):
    list_of_lists_of_indices_of_elements = []
    for element in elements:
        list_of_lists_of_indices_of_elements.append([i for i, value in enumerate(list) if value == element])
    return list_of_lists_of_indices_of_elements

def get_list_items_from_idx_list(list, indices):
    return [list[i] for i in indices]

def check_if_lists_in_list_have_same_length(list_of_lists):
    return all(len(i) == len(list_of_lists[0]) for i in list_of_lists)

def add_prefix_to_list_of_strings(lst, prefix):
    return [prefix + x for x in lst if not str(x) == "nan"]

def search_in_edges(edges, from_to):
    res = []
    for touple in from_to:
        for edge in edges:
            if touple[0] == edge[0] and touple[1] == edge[1]:
                res.append(edge)
    return res

def search_in_list_of_dics_based_on_key(list, key, value):
    for i, dict in enumerate(list):
        if dict[key] == value:
            return i

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    plt.xticks(np.arange(data.shape[1]), labels=col_labels, fontsize=6)
    plt.yticks(np.arange(data.shape[0]), labels=row_labels, fontsize=6)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-50, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    # ax.spines[:].set_visible(False)

    # ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    # ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    # ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

def func(x, pos):
    return "{:.2f}".format(x).replace("0.00", "0").replace("1.00", "1")

def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)

def model_tree_to_nested_dics(model, labels):
    max_node_id = np.max(model.children_) + 1
    ii = itertools.count(model.n_leaves_)
    tree = [{'node_id': next(ii), 'left': x[0], 'right':x[1], 'length':model.distances_[i]} for i, x in enumerate(model.children_)]
    nested_dics = nested_tree(tree, max_node_id, model.n_leaves_, labels)
    return nested_dics

def nested_tree(tree, node_id, n, labels):
    j = search_in_list_of_dics_based_on_key(tree, 'node_id', node_id)
    left = tree[j]['left']
    right = tree[j]['right']
    length = tree[j]['length']
    # print('left = ', left)
    # print('right = ', right)
    tree_dict = {'name' : '', 'length' : length, 'branchset' : []}

    if left >= n:
        tree_dict['branchset'].append(nested_tree(tree, left, n, labels))
    else:
        tree_dict['branchset'].append({'name' : labels[left], 'length' : length, 'branchset' : []})
        
    if right >= n:
        tree_dict['branchset'].append(nested_tree(tree, right, n, labels))
    else:
        tree_dict['branchset'].append({'name' : labels[right], 'length' : length, 'branchset' : []})

    return tree_dict

def check_if_frequency_is_empty(f):
    if f != '':
        return "{:.2f}%".format(float(f))
    else:
        return ''


def tic(message):
    global start_time
    if message:
        print(message)
    start_time = time.time()
    return

def toc(message):
    global start_time
    if message:
        print(message)
    print("--- %s seconds ---" % (time.time() - start_time))

def coordinate_belongs_to_gene_color(coordinate, gene_coordinates, fade = False):
    found = False
    for i in gene_coordinates.index:
        start = gene_coordinates.loc[i, 'gene_start']
        end = gene_coordinates.loc[i, 'gene_end']
        
        if int(start) <= int(coordinate)+1 <= int(end):
            found = True
            break
    if found:
        # hsl = genes_coordinates.loc[i, 'gene_color' if not fade else 'gene_color_faded']
        # return "#" + "".join("%02X" % round(i*255) for i in colorsys.hls_to_rgb(hsl[0]/360, hsl[2]/100, hsl[1]/100))
        return gene_coordinates.loc[i, 'gene_color']
    else:
        return 'black' if not fade else '#BABABA'