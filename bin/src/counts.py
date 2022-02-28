from src.helpers import *
from src.vcfs_parser import parser

class count_snps:
    def __init__(self, df, lineages_of_interest, based_on = [1, 2]):
        self.df = df
        self.lineages_of_interest = lineages_of_interest
        self.based_on = based_on

        self.vcfs_indices_path = "data/clinical_variant_files-indices"

    def occurences_of_mutations(self):
        self.counts, self.summary = count_occurences_of_mutations(self.df, \
                                                        self.lineages_of_interest, \
                                                        based_on = self.based_on, \
                                                        samples_threshold = 0, \
                                                        normalize = True)

    def clusters_intersection(self):
        columns = ["cluster_{}".format(k) for k in self.based_on]
        levels = ["level_{}".format(i) for i, value in enumerate(columns)]
        
        s = pd.Series(self.df.groupby(columns).groups, name = 'samples') \
            .reset_index() \
            .rename(columns={i : k for i, k in zip(levels, columns)})
        s['samples'] = s['samples'].apply(list)


        for i in range(len(s)):
            samples = s.loc[i, 'samples']
            d = {key : [] for key in self.df.loc[samples, 'lineage'].unique()}
            for sample in samples:
                d[self.df.loc[sample, 'lineage']].append(sample)
            s.at[i, 'samples'] = d
        
        return s

    def get_snps_from_lineages_clusters_position(self, samples, region, position):
        pp = parser(self.vcfs_indices_path)
        return pp.get_snp_from_pos(region, base_0_to_1(position), samples)

    def set_path(self, path):
        self.vcfs_indices_path = path