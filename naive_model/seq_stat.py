# Required Libraries
import os  # OS interaction
import re  # Regular expressions
import glob  # Filename pattern matching
import subprocess  # Process management
import argparse  # Command-line argument parsing
from Bio import SeqIO  # Sequence I/O
from Bio import Align  # Sequence alignment
import numpy as np  # Numerical operations
from sklearn.cluster import KMeans  # Clustering algorithm
from sklearn.decomposition import PCA  # Dimensionality reduction
from sklearn.feature_extraction.text import CountVectorizer  # Text to count matrix
from sklearn.metrics.pairwise import euclidean_distances  # Distance calculation
import matplotlib.pyplot as plt  # Plotting library


# Function to extract k-mers from sequences
def get_kmers(sequence, k=5):
    """
    Generate k-mers from a sequence.

    Args:
        sequence (str): The input sequence.
        k (int): The length of each k-mer.

    Returns:
        list: A list of k-mers.
    """
    return [sequence[i:i+k] for i in range(len(sequence)-k+1)]

def cluster_seq(sequences, fname):
    """
    Cluster sequences using k-mer frequency vectors and plot the clusters.

    Args:
        sequences (list): List of sequences to be clustered.
        synthesized_sequence (str): The synthesized sequence to be highlighted.
        fname (str): The name of the output file (without extension).
        output_dir (str): The directory where the output file will be saved.
    """
    # Create a list of k-mer strings for all sequences
    kmer_list = [' '.join(get_kmers(seq)) for seq in sequences]

    # Use CountVectorizer to convert k-mer strings to frequency vectors
    vectorizer = CountVectorizer()
    features = vectorizer.fit_transform(kmer_list).toarray()
    # Apply K-means clustering
    num_clusters = 3  # Adjust based on your data
    kmeans = KMeans(n_clusters=num_clusters, random_state=42)
    labels = kmeans.fit_predict(features)
    # Dimensionality reduction using PCA
    pca = PCA(n_components=2)
    X = pca.fit_transform(features)
    # Plot the clusters
    plt.figure(figsize=(10, 6))
    # Plot each cluster with a different color
    for cluster in range(num_clusters):
        cluster_points = X[labels == cluster]
        plt.scatter(cluster_points[:, 0], cluster_points[:, 1], label=f'Cluster {cluster}', alpha=0.8)
    plt.title(f'Clusters in {fname} Data')    
    # Highlight the synthesized sequence
    synthesized_features = vectorizer.transform([synthesized_kmer_string]).toarray()
    # Find the closest point to the synthesized sequence
    distances = euclidean_distances(features, synthesized_features)
    synthesized_index = np.argmin(distances)
    plt.scatter(X[synthesized_index, 0], X[synthesized_index, 1], color='red', marker='*', s=300)
    plt.legend()
    plt.savefig(f'{output_dir}/{fname}.png', dpi=300, bbox_inches='tight')


def aligned_strand(indices, seqA, seqB):
    """
    Given the indices of alignment, returns the best reconstructed sequence
    """

    aligned_strand = ["-" for i in range(len(seqA))]
    for i,j in zip(indices[0], indices[1]):
        if i != -1 and j!=-1:
            aligned_strand[i] = seqB[j]
    return "".join(aligned_strand)
    

def align(seqA, seqB):
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


    #return (alignment.format(), target, query, identities, mismatches, length)
    return aligned_strand(alignment.indices, seqA, seqB)
    #return target, query

if __name__ == "__main__":
    # Command-line argument parsing
    parser = argparse.ArgumentParser(description="Process Fastq Files.")
    parser.add_argument('-i', '--fastq_dir', type=str, required=True, help="Directory containing fastq sequences")
    parser.add_argument('-l', '--leading', type=str, default="15", help='Trimmomatic parameters: leading')
    parser.add_argument('-t', '--trailing', type=str, default="15", help='Trimmomatic parameters: trailing')
    parser.add_argument('-m', '--minlength', type=str, default="75", help='Trimmomatic parameters: MINLEN')
    parser.add_argument('-s', '--skip', type=str, default="True", help='Skip Trimmomatic or Not')
    args = parser.parse_args()

    # Read Fastq Filepaths
    fastq_paths = glob.glob(f"{args.fastq_dir}/*.fastq*")

    # Create output directory
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Read reference sequence
    reference = str(SeqIO.read("reference.fa", "fasta").seq)
    synthesized_kmer_string = ' '.join(get_kmers(reference))

    # Process and write output
    with open("Final_Output.txt", "w") as out_file:
        for fastq_path in fastq_paths:
            # Extract base filenames
            fname = os.path.basename(fastq_path).split(".fastq")[0]
            out_file.write(f"@FILENAME: {fname}\n")
            
            p1 = subprocess.run(["fastqc", "--memory", "1024", fastq_path, "-o", output_dir])

            if args.skip.lower() == "false":
                numthreads = "4"
                TRUSEQ = "TruSeq3-SE.fa"
                trim_path = f"{output_dir}/{fname}_trimmed.fastq"
                p2 = subprocess.run(["TrimmomaticSE", "-threads", numthreads, "-phred33",
                                fastq_path, trim_path,
                                f"ILLUMINACLIP:{TRUSEQ}:2:30:10",
                                f"LEADING:{args.leading}", f"TRAILING:{args.trailing}",
                                "SLIDINGWINDOW:4:15", f"MINLEN:{args.minlength}", "CROP:75"])
                
                p3 = subprocess.run(["fastqc", trim_path, "-o", output_dir])
                fastq_path = trim_path
        
            sequences = [str(record.seq) for record in SeqIO.parse(fastq_path, "fastq")]
            try:
                # Generate cluster
                cluster_seq(sequences, fname)
                # Iterate over each sequence, align with reference and calculate stats
                for i, sequence in enumerate(sequences, start=1):
                    alignment, target, query, identities, mismatches, length = align(reference, sequence)
                    insertions = target.count("-")
                    deletions = query.count("-")
            
                    identity_percent = round((identities / length) * 100, 3)
                    sub_percent = round((mismatches / length) * 100, 3)
                    ins_percent = round((insertions / length) * 100, 3)
                    del_percent = round((deletions / length) * 100, 3)
            
                    out_file.write(f"#SEQ {i}. {sequence}\n")
                    out_file.write(f"Alignment:\n{alignment}\n")
                    out_file.write(f"\tIdentities (%): {identity_percent}%\n")
                    out_file.write(f"\tSubstitutions (%): {sub_percent}%\n")
                    out_file.write(f"\tInsertions (%): {ins_percent}%\n")
                    out_file.write(f"\tDeletions (%): {del_percent}%\n")

            except ValueError:
                continue
