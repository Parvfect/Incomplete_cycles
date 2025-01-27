
from Bio import SeqIO
from seq_stat import align
import pandas as pd
import numpy as np
import random


def parse_biopython(input_fastq):
    for record in SeqIO.parse(input_fastq, 'fastq'):
        yield record

def parse_fastq_data(fastq_filepath):
    sequenced_strands = []
    for i, record in enumerate(parse_biopython(fastq_filepath)):
        sequenced_strands.append(record)

def postprocess_badread_sequencing_data(fastq_filepath, synthesized_padded_dict=None, reverse_oriented=False, filter=False):
    """
    The record description contains the strand starting, ending and orientation
    """
    sequenced_strands = []
    for i, record in enumerate(parse_biopython(fastq_filepath)):
        strand = str(record.seq)
        # Correcting orientation if it is wrong
        if reverse_oriented:
            try:
                orientation = record.description.split()[1].split(',')[2]
                if orientation == '-strand':
                    strand = strand[::-1]
            except:
                continue

        # Aligning to the target strand if we are filtering        
        if filter:
            try:
                strand_id = record.description.split()[1].split(',')[0]
                if strand_id in synthesized_padded_dict.keys():
                        target_strand = synthesized_padded_dict[strand_id]
                else:
                    continue

                aligned, identity, indices = align(target_strand, strand)

                if identity > 0.7:
                    sequenced_strands.append(strand)
            except:
                continue

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