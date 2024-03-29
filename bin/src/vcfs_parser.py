import numpy as np
import pandas as pd

import hgvs.parser
from hgvs.exceptions import HGVSError
hp = hgvsparser = hgvs.parser.Parser()

from pysam import VariantFile
from os import listdir
from os.path import isfile, join, split
import time

from src.helpers import *

class parser:

    def __init__(self, clinical_data_path, metadata = 0, lineages_data = 0):
        self.clinical_data_path = clinical_data_path
        self.metadata = metadata
        self.lineages_data = lineages_data
        self.ref_len = 29903

    def convert_to_bin(self):
        if (not self.metadata.empty) and (not self.lineages_data.empty):
            self.vcfs_paths, self.vcfs_prefixes = self.__get_files_paths_and_prefixes(self.clinical_data_path)
            # self.vcfs_paths = [f for f in listdir(self.clinical_data_path) if isfile(join(self.clinical_data_path, f))]
            # self.vcfs_prefixes = [k.split('_', 1)[0] for k in self.vcfs_paths]
            self.data = np.zeros((len(self.vcfs_paths), self.ref_len), dtype = 'int8') # init a bin table
            
            print(len(self.vcfs_paths), 'vcfs have been loaded.')
            print('Processing...')
            start = time.time()
            self.__run_through_files()
            print('Processing done in', time.time() - start, 'secs.')
        return

    def __get_files_paths_and_prefixes(self, dir):
        vcfs_paths = [f for f in listdir(dir) if isfile(join(dir, f))]
        vcfs_prefixes = [k.split('_', 1)[0] for k in vcfs_paths]
        return vcfs_paths, vcfs_prefixes
        
    def __is_there_metadata(self, idx):
        if idx in self.metadata.index and pd.notna(self.metadata.loc[idx, 'lineage']):
            return True
        else:
            return False
        
    def __annotate_old(self, data, i, pos, stop, flag):
        if pos == stop:
            pos = pos - 1        

        data[i, pos : stop] = flag
        return data
    
    def __annotate(self, data, i, rec, flag):
        pos = rec.pos
        stop = rec.stop
#         print('ref len =', len(rec.ref))
#         print('alt len =', len(rec.alts[0]))   #...more than one
        if len(rec.ref) == len(rec.alts[0]):
#             print('Substitution')
            pos = pos - 1
        elif len(rec.ref) > len(rec.alts[0]):
#             print('Deletion')
            pass
        elif len(rec.ref) < len(rec.alts[0]):
#             print('Insertion')
            return data
        
        data[i, pos : stop] = flag
        return data

    def __parse(self, record, i, lineage_mutations, is_there_metadata):
    #     retrive the aa information from the ANN field
        try:
            ann = record.info['ANN']
        except KeyError as e:
            print('Something went wrong!')
            print(e)
            return
        ann = [k.split('|') for k in ann]
        flag = 1
        if ann[0][10] != '' and is_there_metadata:
            try:
                v = hp.parse_hgvs_variant(record.chrom + ':' + ann[0][10])
            except HGVSError as e:
                print(e)
                return
            
            if v.posedit.edit.type == 'sub':
                pos_aa = v.posedit.pos.start.base
                ref_aa = v.posedit.pos.start.aa
                alt_aa = v.posedit.edit.alt
                if (not lineage_mutations.loc[(lineage_mutations['ref_aa'] == ref_aa) & \
                                          (lineage_mutations['alt_aa'] == alt_aa) & \
                                          (lineage_mutations['codon_num'] == pos_aa)].empty):
                    flag = 2
            elif v.posedit.edit.type == 'del':
                start_aa = v.posedit.pos.start.base
                end_aa = v.posedit.pos.end.base
                if (not lineage_mutations.loc[(lineage_mutations['codon_num'] == start_aa) & \
                                          (lineage_mutations['codon_end'] == end_aa)].empty):
                    flag = 2
            elif v.posedit.edit.type == 'ins':
                flag = 1
#     #             print('Insertion')
#     #             Do nothing as the insertion change the length
            elif v.posedit.edit.type == 'dup':
                flag = 1
# #                 print(record.chrom + ':' + ann[0][10])
# #                 print(record)
#     #             print('Duplication')
#     #             Do nothing as the duplication change the length the same way the insertion does
            elif v.posedit.edit.type == 'delins':
                flag = 1
#                 print(record.chrom + ':' + ann[0][10])
#                 print('start =', record.start)
#                 print('pos =', record.pos)
#                 print('stop =', record.stop)
#                 print(record)
            elif v.posedit.edit.type == 'fs':
                flag = 1
#                 print(record.chrom + ':' + ann[0][10])
#                 print('start =', record.start)
#                 print('pos =', record.pos)
#                 print('stop =', record.stop)
#                 print(record)
            else:
                flag = 1
#                 print(v.posedit.edit.type)
    #             print(record.chrom + ':' + ann[0][10])
        else:
            flag = 1
#             print('There is not metadata')

        self.data = self.__annotate(self.data, i, record, flag)
        return
    
    def __run_through_files(self):
        
        for i, file in enumerate(self.vcfs_paths):
            bcf_in = VariantFile(join(self.clinical_data_path, file))  # auto-detect input format
#             bcf_out = VariantFile('-', 'w', header=bcf_in.header)

            idx = self.vcfs_prefixes[i]
        #     print(idx)
            is_there_metadata = self.__is_there_metadata(idx)
            if is_there_metadata:
                ln = self.metadata.loc[idx, 'lineage']
                lin_mutations = self.lineages_data.loc[self.lineages_data['lineage'] == ln]

            for rec in bcf_in.fetch():
                self.__parse(rec, i, lin_mutations, is_there_metadata)
        #     break
        self.data = pd.DataFrame(data = self.data, index = self.vcfs_prefixes, dtype = 'int8')

    def get_snp_from_pos(self, chrom, position, samples):
        vcfs_paths, vcfs_prefixes = self.__get_files_paths_and_prefixes(self.clinical_data_path)
        
        try:
            idx = [vcfs_prefixes.index(k) for k in samples]
        except:
            print("Value is not in the list.")
            return
        
        vcfs_paths = [vcfs_paths[k] for k in idx]
        total = 0
        snps = {}
        for i, file in enumerate(vcfs_paths):
            bcf_in = VariantFile(join(self.clinical_data_path, file))  # auto-detect input format

            # print(file)
            region_string = "{chrom}:{start}-{stop}".format(chrom = chrom, \
                                                            start = position, \
                                                            stop = position)
            for rec in bcf_in.fetch(region = region_string):
                total = total + 1
                if len(rec.ref) == len(rec.alts[0]):
                    # print('Substitution')
                    ref = rec.ref
                    alt = rec.alts[0]
                    snps = increase_value_of_dict_by_value(snps, ref + '/' + alt, 1)
                elif len(rec.ref) > len(rec.alts[0]):
                    # print('Deletion')
                    ref = rec.ref[position - rec.pos]
                    alt = '-'
                    snps = increase_value_of_dict_by_value(snps, ref + '/' + alt, 1)
                elif len(rec.ref) < len(rec.alts[0]):
                    print('Insertion')                
        return total, snps
