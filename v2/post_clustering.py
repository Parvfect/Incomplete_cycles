
from utils import get_recovery_percentage

def get_kmers(seq, k):
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]

def remove_adapter(subseq, adapter):
    """
    Given the strand, remove the largest direct match of the starting adapter.
    """
    
    kmer_length = len(subseq)

    while kmer_length >=3:
        kmers = get_kmers(subseq, kmer_length)

        if any([i for i in kmers if i in adapter]):
            for ind, i in enumerate(kmers):
                if i in adapter:    
                    print(f"Length {kmer_length} adapter found at position {ind}")
                    return ind+kmer_length

        kmer_length -= 1

    return -1


def undress_strand(seq, starting_adapter, ending_adapter, len_original, original_strand=None):
    """
    """

    overhang = len(seq) - len_original
    print(f"Overhang is {overhang}")
    starting_subseq = seq[:overhang]
    first_index = remove_adapter(starting_subseq, starting_adapter)
    seq_ = seq[first_index:]


    if original_strand:
        rec = get_recovery_percentage(seq_, original_strand)
        if rec == 1.0:
            print("Found and recovered")
            return seq_
        else:
            print(rec)
