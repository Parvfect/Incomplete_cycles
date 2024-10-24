

#import pysam
import random
import numpy as np

def read_strands_from_file(file_path):
    # Open the FASTQ file using pysam

    sequences = []

    with pysam.FastxFile(file_path) as fastq_file:
        for entry in fastq_file:
            sequences.append(entry.sequence)

    sequences = [i.replace('\r', '') for i in sequences]

    return sequences


def read_strands_from_file_alt(file_path):

    sequences = []

    with open(file_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith('A') or line.startswith('C') or line.startswith('G') or line.startswith('T'):
            sequences.append(line.strip())

    return sequences


def get_reference_from_file(file_path):

    with open(file_path, 'r') as f:
        strand_id = f.read().split()[0]

    return strand_id[1:]


def get_recovery_percentage(consensus_strand, original_strand):

    return sum([
                1 for i in range(len(consensus_strand))
                if consensus_strand[i] == original_strand[i]]
                ) / len(consensus_strand)


def create_low_quality_model(filename):
    # Quality scores range
    max_quality_score = 100
    quality_scores = np.arange(max_quality_score + 1)
    
    # Create probabilities for low quality
    probabilities = []

    # Generate probabilities, favoring lower scores for low quality
    for score in quality_scores:
        if score == 0:
            probabilities.append(1.0)  # Very low probability for high quality
        elif score <= 20:
            probabilities.append(0.0)  # Higher probabilities for lower scores
        else:
            probabilities.append(0.0)  # Low probabilities for higher scores

    # Normalize probabilities to ensure they sum to 1
    total = sum(probabilities)
    probabilities = [p / total for p in probabilities]
    
    # Create the cumulative probabilities and additional columns
    cumulative_probabilities = np.cumsum(probabilities)
    output_lines = []
    
    for i in range(max_quality_score + 1):
        line = [i] + [probabilities[i]] + [cumulative_probabilities[i]] + [0] * (len(quality_scores) - 3)
        output_lines.append('\t'.join(map(str, line)))
    
    # Write to the file
    with open(filename, 'w') as f:
        f.write('\n'.join(output_lines))
