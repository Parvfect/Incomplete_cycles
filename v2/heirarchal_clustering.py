
from Levenshtein import ratio
from tqdm import tqdm
from aligned_clustering import multiple_alignment_muscle
from cluster_merging import majority_merge
import random


def filter_junk_reads(records, ids=None, similarity_threshold=0.85):
    """
    Removes all the sequences that are not similar to any others.
    Prints out percentage of sequences removed
    """
    filtered_records = []
    filtered_seqs = set()

    for record in tqdm(records):
        for record_ in records:
            
            # Checking if its the same record
            if record is record_:
                continue
            
            # Checking if the record is already in the filtered pool
            if not record_.seq in filtered_seqs:
                if ratio(record.seq, record_.seq) > similarity_threshold:
                    filtered_records.append(record_)
                    filtered_seqs.add(record_.seq)

    print(f"{100 - (len(filtered_records) * 100 / len(records))} percent sequences filtered out")

    return filtered_records


def cluster_trivial(strand_pool, similarity_threshold=0.8):
    """
    Can be improved by doing centroids of the cluster rather than originating string
    """
    clusters = {}
    within_clusters = False

    for ind, i in tqdm(enumerate(strand_pool)):
        for j in clusters.keys():
            if ratio(j, i) > similarity_threshold:
                clusters[j].append(ind)
                within_clusters = True
                break
        if not within_clusters:
            clusters[i] = [ind]

        within_clusters = False
    
    return clusters
                    
def make_prediction(cluster, sample_size=5):
    cluster = random.sample(cluster, sample_size)
    return majority_merge(multiple_alignment_muscle(cluster))
