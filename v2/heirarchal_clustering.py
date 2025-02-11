
from Levenshtein import ratio, distance
from tqdm import tqdm
from aligned_clustering import multiple_alignment_muscle
from cluster_merging import majority_merge
import random
import numpy as np
from utils import get_recovery_percentage, reverse_complement

def get_edit_distance_matrix(strands):
    """
    Returns the edit distance matrix for the strands
    O(n^2)
    """
    n_strands = len(strands)
    edit_distance_matrix = np.zeros([n_strands, n_strands])
    for i in range(n_strands - 1):
        for j in range(i + 1, n_strands):
            edit_distance_matrix[i,j] = edit_distance_matrix[j, i] = distance(strands[i], strands[j])

    return edit_distance_matrix

def calculate_centroid(strands: list[str]):
    edit_distance_matrix = get_edit_distance_matrix(strands)

    distances = [sum(edit_distance_matrix[i, :]) for i in range(len(edit_distance_matrix))]
    return strands[distances.index(min(distances))]

def calculate_centroid(strands: list[str], edit_distance_matrix: np.ndarray):
    distances = [sum(edit_distance_matrix[i, :]) for i in range(len(strands))]
    return strands[distances.index(min(distances))]

def update_distance_matrix(
        added_strand: str, cluster_strands: list[str], distance_matrix: np.ndarray):
    """Adds the distances of the added strand to the cluster distance matrix"""

    # List of lists. keep going till next. if sum < main centroids, start going to next hmm

    new_strand_index = len(cluster_strands)
    for ind, cluster_strand in enumerate(cluster_strands):
        distance_matrix[ind, new_strand_index] = distance_matrix[new_strand_index, ind] = distance(added_strand, cluster_strand)

    return distance_matrix


def cluster_trivial(strand_pool, distance_threshold=40, similarity_threshold=0.85, use_centroids=False, analysis=False,
                    sort_order=True):
    """
    Can be improved by doing centroids of the cluster rather than originating string
    """
    clusters_by_strand = []
    clusters_by_index = []
    centroids = []
    within_clusters = False

    if use_centroids:
        distance_matrices = []

    print(f"Total strands {len(strand_pool)}")

    for ind, strand in enumerate(tqdm(strand_pool)):

        within_clusters = False
        rev_strand = reverse_complement(strand)
        for cluster_ind, centroid in enumerate(centroids):
            
            # If it is greater than or equal to the similarity threshold
            if distance(centroid, strand) <= distance_threshold:
                within_clusters = True
            elif distance(centroid, rev_strand) <= distance_threshold:
                within_clusters = True
                strand = rev_strand

            #if ratio(centroid, strand) >= similarity_threshold:# or ratio(centroid, rev_strand) >= similarity_threshold:
            #if distance.Levenshtein.distance(centroid, strand) < similarity_threshold:
            #if get_recovery_percentage(centroid[:14], strand[:14]) >= similarity_threshold:
            if within_clusters:            
                # Add to the cluster
                clusters_by_index[cluster_ind].append(ind)
                clusters_by_strand[cluster_ind].append(strand)

                if use_centroids:
                    # Okay so I need to pairwise check against each element until the edit distance sum surpasses the previous centroid's. If so, then I store it as much as the previous centroid + an arbitary value more / I have some flag that needs to compute the edit distances of the remaining - could be a matrice thing cause I go iteratively. If its None I suppose

                    # Update cluster distance matrix
                    distance_matrices[cluster_ind] = update_distance_matrix(added_strand=strand, cluster_strands=clusters_by_strand[cluster_ind], distance_matrix=distance_matrices[cluster_ind])
                    
                    # Update centroid
                    centroids[cluster_ind] = calculate_centroid(strands=clusters_by_strand[cluster_ind], edit_distance_matrix=distance_matrices[cluster_ind])
                break
            
        if not within_clusters:
            
            clusters_by_index.append([ind])
            centroids.append(strand)
            clusters_by_strand.append([strand])

            if use_centroids:
                # Create new edit distance matrice to compute the centroid
                distance_matrices.append(np.zeros([1000, 1000]))

        if sort_order:
            if ind % 100 == 0:
                # Sorting clusters by size
                # Sort based on the length of sublists in list1
                clusters_by_strand, clusters_by_index, centroids = zip(*sorted(zip(clusters_by_strand, clusters_by_index, centroids), key=lambda x: len(x[0]), reverse=True))

                # Convert back to lists
                clusters_by_strand = list(clusters_by_strand)
                clusters_by_index = list(clusters_by_index)
                centroids = list(centroids)
                
    if analysis:
        return clusters_by_index, centroids, distance_matrices
    
    
    print(f"Number of clusters = {len(centroids)}")

    return clusters_by_index, clusters_by_strand, centroids
                    

def make_prediction(cluster, sample_size=5):
    cluster = random.sample(cluster, sample_size)
    return majority_merge(multiple_alignment_muscle(cluster))
