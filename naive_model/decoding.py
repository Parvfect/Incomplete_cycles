

import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from collections import Counter
from itertools import product
from seq_stat import align, cluster_seq


def get_base_votes(sequenced_strands: list[str], strand_length:int):
        """Get the votes per position of the strand"""

        A_votes = np.zeros(strand_length)
        T_votes = np.zeros(strand_length)
        C_votes = np.zeros(strand_length)
        G_votes = np.zeros(strand_length)

        for i in range(strand_length):
            for j in sequenced_strands:
                if i >= len(j):
                    break
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


def consensus_decoding(sequenced_strands: list[str], strand_length: int):
        """
        Sequenced copies of the same strand, using majority voting to decode.
        Introduced alignment option to align the strands to the original
        strand before voting, so only the aligned bases are considered.
        """

        A_votes, T_votes, C_votes, G_votes = get_base_votes(
            sequenced_strands, strand_length)

        consensus_strand = ""
        bases = ['A', 'T', 'C', 'G']
        for i in range(strand_length):
            votes = [A_votes[i], T_votes[i], C_votes[i], G_votes[i]]
            consensus_strand += bases[np.argmax(votes)]

        return consensus_strand

def get_longest_strand(sequenced_strands, strand_length):

    longest_strand = ""
    for strand in sequenced_strands:
        if len(strand) > len(longest_strand) and len(strand) <= strand_length:
            longest_strand = strand
            
    return longest_strand


def aligned_consensus(sequenced_strands, strand_length):


    longest_strand = get_longest_strand(sequenced_strands, strand_length)

    aligned_strands = [
                    align(sequenced_strand, longest_strand)
                    if len(sequenced_strand) > 20
                    else "" * strand_length
                    for sequenced_strand in sequenced_strands]
    
    return consensus_decoding(aligned_strands, strand_length)

def weighted_aligned_cluster(sequenced_strands, strand_length, n_clusters=3, consensus=False):
    
    # Sample DNA sequences
    dna_sequences = sequenced_strands

    # Step 1: Define all possible k-mers (e.g., 2-mers here)
    kmers = [''.join(p) for p in product("ATCG", repeat=2)]

    # Step 2: Define a function to count k-mers in a sequence
    def kmer_count(sequence, k, kmers):
        counts = Counter([sequence[i:i+k] for i in range(len(sequence) - k + 1)])
        return np.array([counts[kmer] for kmer in kmers])

    # Step 3: Create a feature matrix of k-mer frequency vectors
    feature_matrix = np.array([kmer_count(seq, 2, kmers) for seq in dna_sequences])

    # Step 4: Apply KMeans clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    labels = kmeans.fit_predict(feature_matrix)
    print(labels)
    print("What is going on?")

    # Seperating the strands based on their label
    clusters = [[] for i in range(len(labels))]

    for label, strand in zip(labels, sequenced_strands):
        if len(strand) == 0:
            continue
        clusters[label].append(strand)

    print(clusters)

    consensus_strands = []

    for i in clusters:
        if any([len(j) == 0 for j in i]):
            continue
        consensus_strands.append(aligned_consensus(i, strand_length))
        

    if not consensus:
        return aligned_consensus(consensus_strands, strand_length)

    return consensus_decoding(consensus_strands, strand_length)
   
