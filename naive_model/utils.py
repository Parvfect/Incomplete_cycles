

import pysam

def read_strands_from_file(file_path):
    # Open the FASTQ file using pysam

    ids = []
    sequences = []

    with pysam.FastxFile(file_path) as fastq_file:
        for entry in fastq_file:

            ids.append(entry.name)
            sequences.append(entry.sequence)

    return ids, sequences

ids, sequences = read_strands_from_file('synthesized_strands.txt')

