
from Bio import SeqIO, Align
from seq_stat import align
import pandas as pd
import numpy as np
import random
import json
from tqdm import tqdm
from Levenshtein import distance
import regex as re
import matplotlib.pyplot as plt


def reverse_complement(dna: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna))

def parse_biopython(input_fastq):
    for record in SeqIO.parse(input_fastq, 'fastq'):
        yield record

def get_fastq_records(fastq_filepath):
    return [record for record in parse_biopython(fastq_filepath)]

def postprocess_badread_sequencing_data(fastq_filepath, synthesized_padded_dict=None, reverse_oriented=False, filter=False):
    """
    The record description contains the strand starting, ending and orientation
    Returns strands and records (unchanged)
    """
    records = []
    for i, record in enumerate(tqdm(parse_biopython(fastq_filepath))):
        
        """
        # Correcting orientation if it is wrong
        if reverse_oriented:
            try:
                orientation = record.description.split()[1].split(',')[1]
                if orientation == '-strand':
                    strand = reverse_complement(strand)
            except:
                continue
        """
        records.append(record)

    return records

def post_process_results(recoveries_strands, capping_flags, coupling_rates):

    columns = [
    'capping',
    'coupling_rate',
    'pool_recovery'
    ]
    
    return pd.DataFrame(np.array([capping_flags, coupling_rates, recoveries_strands]).T, columns=columns)
    

def read_synthesized_strands_from_file(file_path, ids=False):
    """Not returning ids currently"""

    sequences = []
    ids = []

    with open(file_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith('>'):
            ids.append(line[1:].strip())

        if line.startswith('A') or line.startswith('C') or line.startswith('G') or line.startswith('T'):
            sequences.append(line.strip())

    if ids:
        return sequences, ids
    
    return sequences

def read_fasta_data(fasta_filepath, ids=False):
    if ids:
        sequences, ids = read_synthesized_strands_from_file(fasta_filepath, ids=ids)
    sequences = read_synthesized_strands_from_file(fasta_filepath, ids=ids)
    return sequences

def get_original_strands(original_strand_filepath, plain=False):
    ids = []
    coupling_rates = []
    capping_flags = []
    strands = []

    with open(original_strand_filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                split_line = line.split()
                
                if len(split_line) > 1:
                    ids.append(split_line[0])
                    coupling_rates.append(split_line[1])
                    capping_flags.append(split_line[2])
                else:
                    strands.append(split_line[0])
    if plain:
        return strands
    return ids, coupling_rates, capping_flags, strands
    
def get_recovery_percentage(consensus_strand, original_strand):
    """Gets recovery percentage based on two strands. Chooses the length of the original strand to evaluate identity"""

    min_length = min(len(original_strand), len(consensus_strand))
    return sum([
                1 for i in range(min_length)
                if consensus_strand[i] == original_strand[i]]
                ) / len(original_strand)

def create_fasta_file(ids, strands, output_filepath='output.fasta'):
    with open(output_filepath, 'w') as f:
        for i, strand in enumerate(strands):
            f.write(f">{ids[i]}\n")
            f.write(strand + '\n\n')

def create_random_strand(strand_length):
    strand = ""
    base_choices = ['A', 'C', 'T', 'G']

    return "".join([random.choice(base_choices) for i in range(strand_length)])


def load_json_file(json_filepath):
    with open(json_filepath, 'r') as f:
        dict_ = json.load(f)
    return dict_

def get_badread_strand_id(record):
    return record.description.split()[1].split(',')[0]


def get_best_candidates_and_recoveries(original_strands, candidates):
    """
    For a given set of strands, finds the best candidates and returns a dictionary with
    the recoveries, the number of fully recovered strands, the best set of candidates and the
    original strands
    """

    fully_recovered_strands = 0
    recoveries = []
    partially_recovered_recoveries = []
    best_candidates = []

    for strand in tqdm(original_strands):
        if strand in candidates:
            fully_recovered_strands += 1
            recoveries.append(1.0)
            best_candidates.append(strand)
        else:
            best_recovery_within_candidates = 0.0
            best_candidate = ""
            for candidate in candidates:
                strand_recovery = get_recovery_percentage(
                    consensus_strand=candidate, original_strand=strand)

                if strand_recovery > best_recovery_within_candidates:
                    best_recovery_within_candidates = strand_recovery
                    best_candidate = candidate

            recoveries.append(best_recovery_within_candidates)
            partially_recovered_recoveries.append(best_recovery_within_candidates)
            best_candidates.append(best_candidate)

    return {
        "recoveries": recoveries,
        "fully_recovered_strands": fully_recovered_strands,
        "partially_recovered_recoveries": partially_recovered_recoveries,
        "best_candidates": best_candidates,
        "original_strands": original_strands
    }


def get_aligned_identity(seqA, seqB):
    aligned, identity, indices = align(seqA, seqB)
    return identity

def count_ids_errors(str1, str2):
    edit_operations = Levenshtein.editops(str1, str2)
    
    insertions = sum(1 for op in edit_operations if op[0] == 'insert')
    deletions = sum(1 for op in edit_operations if op[0] == 'delete')
    substitutions = sum(1 for op in edit_operations if op[0] == 'replace')

    return {'Insertions': insertions, 'Deletions': deletions, 'Substitutions': substitutions}


def align(seqA, seqB, identity=True):
    """
    Performs pairwise sequence alignment between two sequences.

    Args:
        seqA (str): The first sequence to be aligned.
        seqB (str): The second sequence to be aligned.

    Returns:
        tuple: A tuple containing the aligned sequences (target and query).
    """
    # Initialize the PairwiseAligner object
    aligner = Align.PairwiseAligner()
    
    # Set alignment mode to global
    aligner.mode = "local"
    
    # Set match, mismatch, and gap scores
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.gap_score = -5
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    
    # Perform sequence alignment and get the first (best) alignment
    alignment = aligner.align(seqA, seqB)[0]

    # Get alignment as lines
    alignment_lines = alignment.format().strip().split("\n")

    target_starting_index, target_ending_index = 0, 0
        
    # Regex pattern
    pattern = re.compile(r'[0-9\s]')

    # Extract the target sequence
    target_line = "".join([line.split("target")[1] for line in alignment_lines if line.strip().startswith("target")])
    target = re.sub(pattern, "", target_line)
    

    # Extract the query sequence
    query_line = "".join([line.split("query")[1] for line in alignment_lines if line.strip().startswith("query")])
    query = re.sub(pattern, "", query_line)
    identities = alignment.counts()[1]
    mismatches = alignment.counts()[2]
    length = alignment.length

    if identity:
        return identities/length

    #return (alignment.format(), target, query, identities, mismatches, length)
    return alignment

def len_histogram(arr: list[list]):
    plt.hist([len(i) for i in arr])
    plt.show()

def get_sort_by_sublists_length(main_list):
    # Get sorted indices based on the length of sublists in descending order
    return sorted(range(len(main_list)), key=lambda i: len(main_list[i]), reverse=True)


def get_sample_statistics(records, strand_ids, original_strands, original_strand_ids=None, distance_threshold=10, reference=False):
    """
    Given a sample of the data, I want to get the following
    1. Number of strands that don't match to shit
    2. Average number of unique matches per strand
    3. Std number of unique matches per strand
    3. Total number of good strands
    """

    strands_by_index = np.zeros(len(original_strands))
    reverse_strands = 0
    straight_strands = 0
    unmatched = 0
    
    for read, strand_id in tqdm(zip(records, strand_ids), total=len(records)):

        #seq = str(record.seq)
        seq = read
        revseq = str(reverse_complement(seq))
        found_flag = False

        if reference:
            try:
                #strand_id = get_badread_strand_id(record)
                #synthesized_id = strand_ids_synthesized[strand_id]
                #index = original_strand_ids.index(synthesized_id)
                index = original_strand_ids.index(strand_id)
                strand = original_strands[index]

                if distance(seq, strand) <= distance_threshold:
                    strands_by_index[index] += 1
                    straight_strands += 1
                    found_flag = True

                elif distance(revseq, strand) <= distance_threshold:
                    strands_by_index[index] += 1
                    reverse_strands += 1
                    found_flag = True

            except:
                continue

        
        else:
            for ind, strand in enumerate(original_strands):
                if distance(seq, strand) <= distance_threshold:
                    strands_by_index[ind] += 1
                    straight_strands += 1
                    found_flag = True
                elif distance(revseq, strand) <= distance_threshold:
                    strands_by_index[ind] += 1
                    reverse_strands += 1
                    found_flag = True
        
        if not found_flag:
            unmatched += 1

    return {
        'distance_threshold': distance_threshold,
        'strands_by_index': strands_by_index,
        'mean_strands_per_index': float(np.mean(strands_by_index)),
        'std_strands_per_index': float(np.std(strands_by_index)),
        'unique_matches': sum([1 for i in strands_by_index if i>0]),
        'n_straight': straight_strands,
        'n_reverse': reverse_strands,
        'unmatched': unmatched
    }



def sample_reads(reads, ids, sampling_rate=None, n_samples=1000):
    """
    Sample reads and ids """

    if sampling_rate:
        n_samples = int(len(reads) * sampling_rate)

    sample_indices = [random.randint(0, len(reads)) for i in range(n_samples)]
    sampled_reads = [reads[i] for i in sample_indices]
    sampled_ids = [ids[i] for i in sample_indices]

    return sampled_reads, sampled_ids

def sort_clusters(clusters, clustered_seqs, centroids):
    
    sort_indices = get_sort_by_sublists_length(clusters)

    sorted_clusters = [clusters[i] for i in sort_indices]
    sorted_centroids = [centroids[i] for i in sort_indices]
    sorted_clustered_seqs = [clustered_seqs[i] for i in sort_indices]

    centroids = sorted_centroids
    clusters = sorted_clusters
    clustered_seqs = sorted_clustered_seqs

    return clusters, clustered_seqs, centroids