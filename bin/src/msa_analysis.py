from Bio import SeqIO
import numpy as np
from scipy.stats import entropy
from sklearn import metrics
from tqdm.notebook import trange
from joblib import Parallel, delayed

class msa_parser:
    def __init__(self, fasta_path):
        self.msa_path = fasta_path

        self.names = []
        self.seqs = []
        self.lengths = []
        self.uniques = []
        if self.msa_path:
            for seq_record in SeqIO.parse(self.msa_path, "fasta"):
                self.names.append(seq_record.id)
                self.seqs.append(str(seq_record.seq))
                self.uniques.append(list(set(str(seq_record.seq))))
                self.lengths.append(len(seq_record))
            
            self.transposed_array_of_sequences = np.array([list(word) for word in self.seqs]).T
        else:
            print('Please provide a path to a fasta file.')

    def calculate_probabilities(self):
        probabilities = []
        for i in range(self.transposed_array_of_sequences.shape[0]):
            uniques, counts = np.unique(self.transposed_array_of_sequences[i], return_counts=True)
            d = {u : c for u, c in zip(uniques, counts)}
            d.pop("-", None)
            s = sum(d.values())
            p = []
            for key, value in d.items():
                p.append(value / s)
            probabilities.append(p)
        self.probabilities = probabilities
        return

    def calculate_entropies(self, n = 0):
        self.entropies = np.array([entropy(prob, base=2) for prob in self.probabilities])
        # msk = np.argwhere(entropies > 0)
        # msk = np.reshape(msk, (msk.shape[0],))
        if n != 0:
            idx = np.argpartition(self.entropies, -n)[-n:]
            msk = idx[np.argsort((-self.entropies)[idx])]

            self.entropies = np.reshape(self.entropies[msk], (self.entropies[msk].shape[0],))
            self.msk_entropies = msk
        return

    def mutualinfo_vect(self, array, length, i):
        res = np.zeros(shape=(length))
        for j in range(i, length):
            res[j] = metrics.normalized_mutual_info_score(array[i], array[j])
        return res

    def calculate_mutual_info(self, array, save_file_path = "", axis = 0):
        if axis == 1:
            array = array.T
        length = array.shape[0]

        mutual_info = Parallel(n_jobs=16)(delayed(self.mutualinfo_vect)(array, length, i) for i in trange(length))

        mutual_info = np.stack(mutual_info)
        
        if save_file_path:
            np.save(save_file_path, mutual_info)

        self.mutual_info = mutual_info
        return