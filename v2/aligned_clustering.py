
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
import subprocess
from random import shuffle
from clustering import create_clusters
from cluster_merging import merge_clusters
import os
# from decoding import consensus_decoding

def multiple_alignment_muscle(cluster, out = False, running_on_hpc = False):
    
    # write cluster to file
    file = open("clm.fasta","w") 
    
    for i,c in enumerate(cluster):
        file.write(">S%d\n" % i)
        file.write(c)
        file.write("\n")
    file.close()
    
    if running_on_hpc:
        muscle_exe = os.path.join(os.environ['HOME'], 'muscle.exe')
    else:
        muscle_exe = r"C:\Users\Parv\Doc\RA\Projects\incomplete_cycles\muscle-windows-v5.2.exe"

    print(running_on_hpc)
    print(muscle_exe)
    output_alignment = "clmout.fasta"

    #!.\muscle-windows-v5.2.exe -align clm.fasta -output clmout.fasta
    output_message = subprocess.run(
        args=[
            f"{muscle_exe}", "-align", "clm.fasta", "-output", "clmout.fasta"
        ]
    )

    msa = AlignIO.read(output_alignment, "fasta")
    if out:
        print(msa)
    alignedcluster = []
    for i in msa:
        alignedcluster += [i.seq]
    return alignedcluster


def align_clusters(trimmed_seqs, clusters, masize = 15, running_on_hpc = False):

    fresults = []
    ### align clusters, generate candidates
    for i, clusterinds in enumerate(clusters):
        cluster = [trimmed_seqs[i] for i in clusterinds]
        if len(cluster) < 3:
            continue
        if len(cluster) > masize:
            for j in range(5):
                shuffle(cluster)
                ma = multiple_alignment_muscle(cluster[:masize])
                fresults.append(ma)
        else:
            ma = multiple_alignment_muscle(cluster[:masize], running_on_hpc)
            fresults.append(ma)

    return fresults


def get_recovery_percentage(consensus_strand: str, original_strand: str):

    min_length = min(len(original_strand), len(consensus_strand))
    return sum([
                1 for i in range(min_length)
                if consensus_strand[i] == original_strand[i]]
                ) / len(original_strand)

def filter_sequences(trimmed_seqs, length_filtering, original_strand_length):
    """Filter seqeunces that go into the aligned clustering based on the length"""

    return [i for i in trimmed_seqs if len(i) > length_filtering * original_strand_length]


def conduct_align_clustering(
        original_strand, trimmed_seqs, trivial_clustering = True,
        display=True, multiple=False, best_recovery=False, length_filtering=0, running_on_hpc = False):
    
    if length_filtering:
        trimmed_seqs = filter_sequences(trimmed_seqs, length_filtering, len(original_strand))

    clusters = create_clusters(
        trimmed_seqs=trimmed_seqs, TRIVIAL=trivial_clustering)

    fresults = align_clusters(
        trimmed_seqs=trimmed_seqs,
        clusters=clusters,
        running_on_hpc=running_on_hpc
    )

    candidates = merge_clusters(
        fresults=fresults
    )

    if not multiple:
        recoveries = [
                get_recovery_percentage(candidate, original_strand)
                for candidate in candidates
                ]
    else:
        recoveries = {strand: [get_recovery_percentage(candidate, strand) for candidate in candidates] for strand in original_strand}
        

    if display:
        print("Evaluating recovery percentage")
        print(f"Best recovery percentage in candidates = {max(recoveries)}")

    if best_recovery:
        return max(recoveries)

    return {
        "candidates":candidates,
        "fresults": fresults,
        "recoveries": recoveries
    }
    
"""
def k_best_candidates_merging(strand: str, candidates: list[str], cutoff_percentage=0.3):
    recoveries = {i: get_recovery_percentage(i, strand) for i in candidates}

    recoveries = sorted(recoveries.items(), key=lambda x: x[1], reverse=True)

    consensus_strands = [i[0] for i in recoveries if i[1] > cutoff_percentage]
    final_consensus = consensus_decoding(consensus_strands, len(strand))

    return get_recovery_percentage(final_consensus, strand)
"""