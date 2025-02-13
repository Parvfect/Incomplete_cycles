
from utils import get_recovery_percentage, reverse_complement
from tqdm import tqdm
from heirarchal_clustering import make_prediction
from Levenshtein import ratio


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
                    return ind+kmer_length

        kmer_length -= 1

    return 0
start_adapter = "AATGTACTTCGTTCAGTTACGTATTGCT" 
end_adapter = "GCAATACGTAACTGAACGAAGT"

def remove_adapters_from_strands(strands, ids, original_strand_length,
                                starting_adapter="AATGTACTTCGTTCAGTTACGTATTGCT",
                                ending_adapter="GCAATACGTAACTGAACGAAGT"):
    """
    Removes the adapters and fixes the orientation of the strand by checking for the adapter found
    """

    cleaned_strands = []
    cleaned_ids = []

    print("Removing adapters")
    for strand, id in tqdm(zip(strands, ids), total=len(strands)):
        overhang = len(strand) - original_strand_length

        if overhang > len(starting_adapter) + len(ending_adapter):
            continue

        if overhang > 3:

            # remove the starting adapter
            start_index = remove_adapter(strand[:len(start_adapter)+ 2], start_adapter)
            strand = strand[start_index:]

        cleaned_strands.append(strand)
        cleaned_ids.append(id)
        
    return cleaned_strands, cleaned_ids

class FiniteStateMachine():
    """Not really, but I like the concept"""

    def __init__(self, original_strands):
        self.original_strands = original_strands
        self.n = len(original_strands)

    def operate(self, guesses, hard=False, reversals=False):
        """Non ordering checks for fsm"""

        assert len(guesses) == self.n
        
        found = 0
        for i in self.original_strands:
            if i in guesses:
                found += 1
                continue
            if reversals:
                if reverse_complement(i) in guesses:
                    found += 1
        
        if found == self.n:
            print("Found all")
            return True

        if hard:
            return False

        else:
            print(f"Found {found}")
            return False
        



def check_clusters(original_strands, original_strand_ids, clusters, sampled_reads, sampled_ids, make_guess=False):
    """
    Given the cluster indices, returns the traceback recoveries (checks for reversals)
    """

    recs = []
    guesses = []
    for cluster in tqdm(clusters):
        strand_id = sampled_ids[cluster[0]]
        
        if strand_id == 'junk_seq' or strand_id == 'random_seq':
            continue
        
        original_strand_index = original_strand_ids.index(strand_id)
        reference = original_strands[original_strand_index]
        cluster_head = sampled_reads[cluster[0]]

        if make_guess:
            guess_sequences = [sampled_reads[i] for i in cluster]
            guess = make_prediction(guess_sequences, min(5, len(guess_sequences)))
            cluster_head = guess

        rec_1 = ratio(reference, cluster_head)
        rec_2 = ratio(reference, reverse_complement(cluster_head))

        if rec_1 > rec_2:
            guesses.append(cluster_head)
        else:
            guesses.append(reverse_complement(cluster_head))
        
        recs.append(max([rec_1, rec_2]))

    return recs, guesses


    