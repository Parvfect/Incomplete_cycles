
from Bio import SeqIO, Align
from seq_stat import align
import pandas as pd
import numpy as np
import random
import json
from tqdm import tqdm
import Levenshtein
import regex as re
import matplotlib.pyplot as plt


def parse_biopython(input_fastq):
    for record in SeqIO.parse(input_fastq, 'fastq'):
        yield record

def parse_fastq_data(fastq_filepath):
    sequenced_strands = []
    for i, record in enumerate(parse_biopython(fastq_filepath)):
        sequenced_strands.append(record)

def get_fastq_records(fastq_filepath):
    return [record for record in parse_biopython(fastq_filepath)]

def postprocess_badread_sequencing_data(fastq_filepath, synthesized_padded_dict=None, reverse_oriented=False, filter=False):
    """
    The record description contains the strand starting, ending and orientation
    """
    sequenced_strands = []
    for i, record in enumerate(tqdm(parse_biopython(fastq_filepath))):
        strand = str(record.seq)
        # Correcting orientation if it is wrong
        sequenced_strands.append(strand)

    return sequenced_strands

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

def read_fasta_data(fasta_filepath):
    sequences = read_synthesized_strands_from_file(fasta_filepath, ids=False)
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
