
import numpy as np
from seq_stat import align, cluster_seq
from sklearn.cluster import KMeans  # Clustering algorithm
from sklearn.decomposition import PCA  # Dimensionality reduction
from sklearn.feature_extraction.text import CountVectorizer  # Text to count matrix
from sklearn.metrics.pairwise import euclidean_distances  # Distance calculation
import matplotlib.pyplot as plt  # Plotting library
from collections import Counter
from scipy.spatial.distance import hamming


bases = ['A', 'T', 'C', 'G']


class NaiveSequencingModel:
    """
    Single strand sequencing model -
    Independent base IDS errors
    """

    def __init__(self, insertion_probability=0.05,
                 deletion_probability=0.04, subsitution_probability=0.02,
                 strand_length=20):
        self.insertion_probability = insertion_probability
        self.deletion_probability = deletion_probability
        self.subsitution_probability = subsitution_probability
        self.strand_length = strand_length


    def simulate_sequencing(self, strand):
        """
        Simulates the sequencing process of a DNA strand by introducing
        random insertions, deletions, and substitutions based on
        predefined probabilities.

        Args:
            strand (str): The original DNA strand to be sequenced.
        Returns:
            str: The sequenced DNA strand after introducing random errors.
       """

        sequenced_strand = ""
        for i in range(len(strand)):

            # If there is an insertion, add a random base
            if np.random.rand() < self.insertion_probability:
                sequenced_strand += np.random.choice(bases)

            # If there is a deletion, skip the base
            if np.random.rand() < self.deletion_probability:
                continue

            # If there is a subsitution, replace the base
            if np.random.rand() < self.subsitution_probability:
                sequenced_strand += np.random.choice(bases)
                continue

            sequenced_strand += strand[i]

        return sequenced_strand

    def simulate_full_sequencing(self, strands):

        sequenced_strands = []
        for strand in strands:
            sequenced_strands.append(self.simulate_sequencing(strand))

        return sequenced_strands

    def get_base_votes(self, sequenced_strands):
        """Get the votes per position of the strand"""

        A_votes = np.zeros(self.strand_length)
        T_votes = np.zeros(self.strand_length)
        C_votes = np.zeros(self.strand_length)
        G_votes = np.zeros(self.strand_length)

        for i in range(self.strand_length):
            for j in sequenced_strands:
                if i >= len(j):
                    continue
                if j[i] == 'A':
                    A_votes[i] += 1
                elif j[i] == 'T':
                    T_votes[i] += 1
                elif j[i] == 'C':
                    C_votes[i] += 1
                elif j[i] == 'G':
                    G_votes[i] += 1
                else:   # In case there is a blank post alignment
                    continue

        return A_votes, T_votes, C_votes, G_votes

    def consensus_decoding(self, sequenced_strands, original_strand,
                           alignment=False):
        """
        Sequenced copies of the same strand, using majority voting to decode.
        Introduced alignment option to align the strands to the original
        strand before voting, so only the aligned bases are considered.
        """

        if alignment:
            sequenced_strands = [
                align(original_strand, sequenced_strand)
                if len(sequenced_strand) > 0
                else "" * len(original_strand)
                for sequenced_strand in sequenced_strands]

        A_votes, T_votes, C_votes, G_votes = self.get_base_votes(
            sequenced_strands)

        consensus_strand = ""
        bases = ['A', 'T', 'C', 'G']
        for i in range(self.strand_length):
            votes = [A_votes[i], T_votes[i], C_votes[i], G_votes[i]]
            consensus_strand += bases[np.argmax(votes)]

        return sum([
        i==j for i, j in zip(original_strand, consensus_strand)
        ])/len(original_strand)

    def get_base_likelihoods(self, A_votes, T_votes, C_votes, G_votes):
        """Given the votes per position of the strand, get the symbol
        likelihood array"""

        assert len(A_votes) == len(T_votes) == len(C_votes) == len(G_votes)

        likelihoods = []
        for i in range(self.strand_length):
            votes = [A_votes[i], T_votes[i], C_votes[i], G_votes[i]]
            likelihoods.append(votes/np.sum(votes))

        return likelihoods

