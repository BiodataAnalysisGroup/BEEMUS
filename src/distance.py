import numpy as np

class hamming_distance:
    
    @staticmethod
    def compute_HammingDistance(X):
        return (2 * np.inner(X-0.5, 0.5-X) + X.shape[1] / 2)
    
    
class jukes_cantor_distance:

    @staticmethod
    def compute_jukes_cantor (seq1, seq2):
        # Jukes Cantor distance formula: (-3/4)ln[1-p*(4/3)]
        p = percent_difference_of_nucleotides(seq1, seq2)
        return -0.75 * math.log(1 - (p*4/3)) if p else 0


    def percent_difference_of_nucleotides (seq1, seq2, nucleobases=set('ACGT')):
        # percentage of nucleotide difference in two sequences

        diff_count = 0 # number of nucleotide differences
        valid_nucleotides_count = 0.0 # number of valid nucleotides (value is float for computing percentage)

        for a, b in zip(seq1, seq2):
            if a in nucleobases and b in nucleobases:
                valid_nucleotides_count += 1
                if a != b: diff_count += 1

        return diff_count / valid_nucleotides_count if valid_nucleotides_count else 0
    
    @staticmethod
    def hamming_to_jukes_cantor(X, length):
        d = np.ndarray(X.shape, dtype='float')
        d.fill(np.inf)
        #The equation does not exist for values of p greater than 3/4, neither for values of length == 0
        if(length == 0):
            return d
        
        p = X / length
        
        msk = (p < 3.0/4)
        d[msk] = -3.0/4 * np.log(1 - 4.0/3 * p[msk])
        
        msk = (d == -0.0)
        d[msk] = -d[msk]
        return d

