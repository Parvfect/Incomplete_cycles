
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
import subprocess
from random import shuffle
from clustering import create_clusters
from cluster_merging import merge_clusters

def multiple_alignment_muscle(cluster, out=False):
    
    # write cluster to file
    file = open("clm.fasta","w") 
    
    for i,c in enumerate(cluster):
        file.write(">S%d\n" % i)
        file.write(c)
        file.write("\n")
    file.close()

    muscle_exe = r"muscle-windows-v5.2.exe" 
    output_alignment = "clmout.fasta"

    #!.\muscle-windows-v5.2.exe -align clm.fasta -output clmout.fasta
    output_message = subprocess.run(
        args=[
            ".\muscle-windows-v5.2.exe", "-align", "clm.fasta", "-output", "clmout.fasta"
        ],
        capture_output=True
    )

    msa = AlignIO.read(output_alignment, "fasta")
    if out:
        print(msa)
    alignedcluster = []
    for i in msa:
        alignedcluster += [i.seq]
    return alignedcluster


def align_clusters(trimmed_seqs, clusters, masize = 15):

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
            ma = multiple_alignment_muscle(cluster[:masize])
            fresults.append(ma)

    return fresults


def get_recovery_percentage(consensus_strand, original_strand):

    min_length = min(len(original_strand), len(consensus_strand))
    return sum([
                1 for i in range(min_length)
                if consensus_strand[i] == original_strand[i]]
                ) / len(original_strand)


def conduct_align_clustering(
        original_strand, trimmed_seqs, trivial_clustering=True,
        display=True):
    
    clusters = create_clusters(
        trimmed_seqs=trimmed_seqs, TRIVIAL=trivial_clustering)

    fresults = align_clusters(
        trimmed_seqs=trimmed_seqs,
        clusters=clusters
    )

    candidates = merge_clusters(
        fresults=fresults
    )

    recoveries = [
            get_recovery_percentage(candidate, original_strand)
            for candidate in candidates
            ]    

    if display:
        print("Evaluating recovery percentage")
        print(f"Best recovery percentage in candidates = {max(recoveries)}")

    return {
        "candidates":candidates,
        "fresults": fresults,
        "recoveries": recoveries
    }
    