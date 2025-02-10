

import math
import re
from collections import Counter

def kmer_entropy(seq, k=3):
    """
    Calculates Shannon entropy based on k-mers of size k.

    :param seq: DNA sequence string
    :param k: k-mer size (default = 3 for trinucleotides)
    :return: Shannon entropy value
    """
    if len(seq) < k:
        return 0  # If sequence is shorter than k, entropy is 0

    # Count k-mer frequencies
    kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
    kmer_counts = Counter(kmers)
    total_kmers = sum(kmer_counts.values())

    # Compute Shannon entropy
    entropy = -sum((count / total_kmers) * math.log2(count / total_kmers) for count in kmer_counts.values())

    return entropy

def sum_entropies(seq):
    entropy_sum = 0
    for i in range(6):
        entropy_sum += kmer_entropy(seq, i)

    return entropy_sum


def has_long_run(seq, run_length=5):
    """
    Checks if a DNA sequence contains a run of the same symbol of length >= run_length.
    
    :param seq: DNA sequence string
    :param run_length: Minimum length of consecutive identical nucleotides to flag
    :return: True if sequence contains a long run, otherwise False
    """
    return bool(re.search(r"(A{%d,}|T{%d,}|C{%d,}|G{%d,})" % (run_length, run_length, run_length, run_length), seq))


def shannon_entropy(seq):
    """Calculate the Shannon entropy of a DNA sequence."""
    freq = {nuc: seq.count(nuc) / len(seq) for nuc in "ATCG"}
    print(freq)
    return -sum(p * math.log2(p) for p in freq.values() if p > 0)

def filter_low_complexity(dna_sequences, threshold=1.5):
    """
    Filters out low-complexity DNA sequences based on Shannon entropy.
    
    :param dna_sequences: List of DNA sequences
    :param threshold: Entropy threshold below which sequences are removed
    :return: Filtered list of sequences
    """
    return [seq for seq in dna_sequences if shannon_entropy(seq) > threshold]