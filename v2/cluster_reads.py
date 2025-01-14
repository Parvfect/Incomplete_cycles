
from aligned_clustering import conduct_align_clustering
from utils import get_original_strands, postprocess_badread_sequencing_data, post_process_results
import argparse
from datetime import datetime
import os
import json

# Parse input filename for reads and original strands. Can filter later. Return recoveries and candidates

parser = argparse.ArgumentParser(
                    prog='Clustering reads for dna storage',
                    epilog='Use parser to set some parameters for clustering')

parser.add_argument('--reads_filepath', type=str, help="Path to reads.fasta")
parser.add_argument('--info_filepath', type=str, help="Path to original_strands.txt. Also the output filepath.")


def get_run_information_from_files(info_filepath):
    """Gets the original strands from generated file"""

    original_strand_ids, coupling_rates, capping_flags, original_strands = get_original_strands(
        original_strand_filepath=os.path.join(info_filepath, 'original_strands.txt'))
    
    # TODO: Get experiment description files
    experiment_description = None

    return original_strand_ids, coupling_rates, capping_flags, original_strands

def extract_reads_from_fastq(reads_filepath, reverse_oriented=True):

    # TODO: Add filter option
    # TODO: Use synthesized padded dict to generate some run statistics
    sequenced_strands = postprocess_badread_sequencing_data(fastq_filepath=reads_filepath,
                                                             reverse_oriented=reverse_oriented)
    return sequenced_strands

def cluster_reads(sequenced_strands, original_strands=None):
    """Conducts aligned clustering for a given read bunch and compared with the original strands"""

    if original_strands:
        recoveries = conduct_align_clustering(
            original_strand=original_strands,
            trimmed_seqs=sequenced_strands,
            display=False,
            multiple=True
        )
    else:
        recoveries = conduct_align_clustering(
            original_strand=original_strands,
            trimmed_seqs=sequenced_strands,
            display=False,
            multiple=True
        )
    return recoveries

def create_clustering_report_file(
        recoveries, output_filepath, original_strands=None, coupling_rates=None, capping_flags=None):

    # Creating results summary if original strands are provided
    if original_strands is not None:
        df = post_process_results(recoveries_strands=list(
            recoveries['recoveries'].values()),
            capping_flags=capping_flags, coupling_rates=coupling_rates)
        df.to_csv(os.path.join(output_filepath, "results_summary.csv"))

    # Writing the recoveries object
    out_file = open(os.path.join(output_filepath, "recoveries.json"), "w")
    json.dump(recoveries, out_file, indent = 6)
    out_file.close()
    

args = parser.parse_args()

if __name__ == '__main__':
    reads_filepath = args.reads_filepath
    info_filepath = args.info_filepath
    original_strands, experiment_description = None, None

    # If no original strands, just output the recoveries
    if info_filepath:
        original_strand_ids, coupling_rates, capping_flags, original_strands = get_run_information_from_files(
            info_filepath=info_filepath)

    sequenced_strands = extract_reads_from_fastq(reads_filepath=reads_filepath)

    recoveries = cluster_reads(sequenced_strands=sequenced_strands, original_strands=original_strands)

    create_clustering_report_file(recoveries=recoveries, output_filepath=info_filepath,
                                  original_strands=original_strands, coupling_rates=coupling_rates,
                                  capping_flags=capping_flags)
