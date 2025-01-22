
from aligned_clustering import conduct_align_clustering
from utils import get_original_strands, postprocess_badread_sequencing_data, post_process_results, read_fasta_data
import argparse
from datetime import datetime
import os
import json
import random
import uuid
import time

# Parse input filename for reads and original strands. Can filter later. Return recoveries and candidates

parser = argparse.ArgumentParser(
                    prog='Clustering reads for dna storage',
                    epilog='Use parser to set some parameters for clustering')

parser.add_argument('--reads_filepath', type=str, help="Path to reads.fasta")
parser.add_argument('--info_filepath', type=str, help="Path to original_strands.txt. Also the default output filepath.")
parser.add_argument('--output_filepath', type=str, help="Output filepath.")
parser.add_argument('--sampling_rate', type=float, help="Sampling rate - defaults to 1.0")
parser.add_argument('--badread_data', action='store_true', help="Badread data flag")
parser.add_argument('--hpc', action='store_true', help="Running on HPC flag")

parser.set_defaults(reads_filepath=None, info_filepath=None, output_filepath=None, sampling_rate=1.0, hpc=False)

def get_run_information_from_files(info_filepath):
    """Gets the original strands from generated file"""

    original_strand_ids, coupling_rates, capping_flags, original_strands = get_original_strands(
        original_strand_filepath=os.path.join(info_filepath, 'original_strands.txt'))
    
    # TODO: Get experiment description files
    experiment_description = None

    return original_strand_ids, coupling_rates, capping_flags, original_strands

def extract_reads(reads_filepath, badread_data_flag=False, reverse_oriented=True, sampling_rate=1.0):

    if not badread_data_flag:
        sequenced_strands = read_fasta_data(reads_filepath)
    else:
        sequenced_strands = postprocess_badread_sequencing_data(fastq_filepath=reads_filepath,
                                                                reverse_oriented=reverse_oriented)
        
    return random.sample(sequenced_strands, int(len(sequenced_strands) * sampling_rate))

def cluster_reads(sequenced_strands, original_strands=None, sampling_rate=1.0, hpc=False):
    """Conducts aligned clustering for a given read bunch and compared with the original strands"""

    if original_strands:
        recoveries = conduct_align_clustering(
            original_strand=original_strands,
            trimmed_seqs=sequenced_strands,
            multiple=True,
            running_on_hpc=hpc
        )
    else:
        recoveries = conduct_align_clustering(
            original_strand=original_strands,
            trimmed_seqs=sequenced_strands,
            multiple=True,
            running_on_hpc=hpc
        )
    return recoveries

def create_clustering_report_file(
        recoveries, output_filepath, original_strands=None, coupling_rates=None, capping_flags=None):

    uid = str(uuid.uuid4())
    timestamp = str(datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H_%M_%S'))
    save_path = os.path.join(output_filepath, timestamp)

    os.mkdir(save_path) # Timestamp make results directory

    # Creating results summary if original strands are provided
    if original_strands is not None:
        df = post_process_results(recoveries_strands=list(
            recoveries['recoveries'].values()),
            capping_flags=capping_flags, coupling_rates=coupling_rates)
        df.to_csv(os.path.join(save_path, f"results_summary_{uid}.csv"))

    # Writing the recoveries object
    out_file = open(os.path.join(save_path, f"recoveries_{uid}.json"), "w")
    json.dump(recoveries, out_file, indent = 6)
    out_file.close()

    # Info file of the run
    

args = parser.parse_args()

if __name__ == '__main__':
    reads_filepath = args.reads_filepath
    info_filepath = args.info_filepath
    sampling_rate = float(args.sampling_rate)
    badread_data_flag = args.badread_data
    output_filepath = args.output_filepath
    hpc = args.hpc

    #print(reads_filepath)
    #print(info_filepath)
    #print()

    if output_filepath is None:
        output_filepath = info_filepath

    original_strands, experiment_description = None, None

    if args.reads_filepath is None:
        print("No reads filepath provided!")
        exit()

    # If no original strands, just output the recoveries
    if info_filepath:
        original_strand_ids, coupling_rates, capping_flags, original_strands = get_run_information_from_files(
            info_filepath=info_filepath)
        
    print(f"Original strands loaded \n {original_strands}\n\n")

    sequenced_strands = extract_reads(reads_filepath=reads_filepath, badread_data_flag=badread_data_flag, sampling_rate=sampling_rate)
    
    print(f"Number of strands in the pool = {len(sequenced_strands)}\n\n")

    recoveries = cluster_reads(sequenced_strands=sequenced_strands, original_strands=original_strands, sampling_rate=sampling_rate, hpc=hpc)

    print(recoveries['recoveries'])

    create_clustering_report_file(recoveries=recoveries, output_filepath=output_filepath,
                                  original_strands=original_strands, coupling_rates=coupling_rates,
                                  capping_flags=capping_flags)
