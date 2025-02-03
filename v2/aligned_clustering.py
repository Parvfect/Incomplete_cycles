
from Bio import AlignIO
import subprocess
from random import shuffle
from clustering import create_clusters
from cluster_merging import merge_clusters
import os
from tqdm import tqdm
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
        # Will need to install muscle on new device and change the path
        muscle_exe = os.path.join(os.environ['HOME'], 'muscle.exe')
    else:
        muscle_exe = r"C:\Users\Parv\Doc\RA\Projects\incomplete_cycles\muscle-windows-v5.2.exe"

    output_alignment = "clmout.fasta"

    #!.\muscle-windows-v5.2.exe -align clm.fasta -output clmout.fasta

    try:
        output_message = subprocess.run(
            args=[
                f"{muscle_exe}", "-align", "clm.fasta", "-output", "clmout.fasta"
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            timeout=120
        )
    except Exception as e:
        print(e)
        return None

    msa = AlignIO.read(output_alignment, "fasta")
    if out:
        print(msa)
    alignedcluster = []
    for i in msa:
        alignedcluster += [i.seq]
    return alignedcluster


def align_clusters(trimmed_seqs, clusters, masize = 15, running_on_hpc = False, min_cluster_length = 70, max_cluster_length=160, cluster_ids=None):

    fresults = []
    candidate_cluster_ids = []

    # Collecting the clusters
    clusters_arr = [[trimmed_seqs[i] for i in clusterinds] for clusterinds in clusters]

    # Filtering the clusters based on the average strand length of the cluster
    #clusters_arr = [
    #    i for i in clusters_arr if sum([len(j) for j in i]) / len(i) > min_cluster_length and sum([len(j) for j in i]) / len(i) < max_cluster_length]
    
    print(f"Number of clusters = {len(clusters_arr)}")

    for i in tqdm(range((len(clusters_arr)))):

        cluster = clusters_arr[i]
        cluster_id = cluster_ids[i]

        if len(cluster) < 3:
            continue
        
        if len(cluster) > masize:
            for j in range(5):
                shuffle(cluster)
                ma = multiple_alignment_muscle(cluster[:masize], running_on_hpc=running_on_hpc)

                if ma is not None:
                    fresults.append(ma)
                    candidate_cluster_ids.append(cluster_id)
        else:
            ma = multiple_alignment_muscle(cluster[:masize], running_on_hpc=running_on_hpc)

            if ma is not None:
                fresults.append(ma)
                candidate_cluster_ids.append(cluster_id)

    return fresults, candidate_cluster_ids


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
        trimmed_seqs, original_strand = None, trivial_clustering = False,
        multiple = False, best_recovery = False, length_filtering = 0, running_on_hpc = False,
        min_cluster_length=50, max_cluster_length=200, nbeg=14, multiple_align=True, target_clusters=None, m=10):
    
    if length_filtering:
        trimmed_seqs = filter_sequences(trimmed_seqs, length_filtering, len(original_strand))

    clusters, cluster_ids = create_clusters(
        trimmed_seqs=trimmed_seqs, TRIVIAL=trivial_clustering, nbeg=nbeg, target_clusters=target_clusters, m=10)
    
    #return clusters

    fresults, candidate_cluster_ids = align_clusters(
        trimmed_seqs=trimmed_seqs,
        clusters=clusters,
        running_on_hpc=running_on_hpc,
        min_cluster_length=min_cluster_length,
        max_cluster_length=max_cluster_length,
        cluster_ids=cluster_ids
    )

    candidates = merge_clusters(
        fresults=fresults
    )

    if original_strand is None:
        {
            "candidates": candidates,
            "candidate_cluster_ids": candidate_cluster_ids,
            "cluster_ids": cluster_ids,
            "clusters": clusters
        }

    if not multiple:
        recoveries = [
                get_recovery_percentage(candidate, original_strand)
                for candidate in candidates
                ]
    else:
        # Iterates through all the original strands and gives the best recovery from all the candidates for all of them
        recoveries = {}
        for strand in original_strand:
            recoveries_strand = [get_recovery_percentage(candidate, strand) for candidate in candidates]
            if len(recoveries_strand) > 0:
                recoveries[strand] = max(recoveries_strand)
            else:
                recoveries[strand] = 0

    if best_recovery and not multiple:
        return max(recoveries)

    return {
        "recoveries": recoveries,
        "candidates": candidates
    }
    
"""
def k_best_candidates_merging(strand: str, candidates: list[str], cutoff_percentage=0.3):
    recoveries = {i: get_recovery_percentage(i, strand) for i in candidates}

    recoveries = sorted(recoveries.items(), key=lambda x: x[1], reverse=True)

    consensus_strands = [i[0] for i in recoveries if i[1] > cutoff_percentage]
    final_consensus = consensus_decoding(consensus_strands, len(strand))

    return get_recovery_percentage(final_consensus, strand)
"""