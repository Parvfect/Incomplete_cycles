
from Bio import SeqIO
from tqdm import tqdm


def reverse_complement(dna: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna))

def parse_biopython(input_fastq):
    for record in SeqIO.parse(input_fastq, 'fastq'):
        yield record

def get_fastq_records(fastq_filepath):
    records = []
    for i, record in enumerate(tqdm(parse_biopython(fastq_filepath))):
        records.append(record)
    return records