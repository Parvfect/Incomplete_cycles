
from Levenshtein import ratio, distance
from tqdm import tqdm
from aligned_clustering import multiple_alignment_muscle
from cluster_merging import majority_merge
import random
import numpy as np

def get_edit_distance_matrix(strands):
    """
    Returns the edit distance matrix for the strands
    O(n^2)
    """
    n_strands = len(strands)
    edit_distance_matrix = np.zeros([n_strands, n_strands])
    for i in range(n_strands - 1):
        for j in range(i + 1, n_strands):
            edit_distance_matrix[i,j] = edit_distance_matrix[j, i] = ratio(strands[i], strands[j])

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

    new_strand_index = len(cluster_strands)
    for ind, cluster_strand in enumerate(cluster_strands):
        distance_matrix[ind, new_strand_index] = distance_matrix[new_strand_index, ind] = distance(added_strand, cluster_strand)

    return distance_matrix


def cluster_trivial(strand_pool, similarity_threshold=0.8, use_centroids=True, analysis=False):
    """
    Can be improved by doing centroids of the cluster rather than originating string
    """
    clusters_by_strand = []
    clusters_by_index = []
    centroids = []
    within_clusters = False

    if use_centroids:
        distance_matrices = []

    for ind, strand in tqdm(enumerate(strand_pool)):

        for cluster_ind, centroid in enumerate(centroids):
            
            # If it is greater than or equal to the similarity threshold
            if ratio(centroid, strand) >= similarity_threshold:
                
                # Add to the cluster
                clusters_by_strand[cluster_ind].append(strand)
                clusters_by_index[cluster_ind].append(ind)
                
                if use_centroids:
                    # Update cluster distance matrix
                    distance_matrices[cluster_ind] = update_distance_matrix(added_strand=strand, cluster_strands=clusters_by_strand[cluster_ind], distance_matrix=distance_matrices[cluster_ind])
                    
                    # Update centroid
                    centroids[cluster_ind] = calculate_centroid(strands=clusters_by_strand[cluster_ind], edit_distance_matrix=distance_matrices[cluster_ind])

                within_clusters = True
                break
            
        if not within_clusters:
            clusters_by_strand.append([strand])
            clusters_by_index.append([ind])
            centroids.append(strand)

            if use_centroids:
                # Create new edit distance matrice to compute the centroid
                distance_matrices.append(np.zeros([1000, 1000]))

        within_clusters = False

    if analysis:
        return clusters_by_index, centroids, distance_matrices
    
    return clusters_by_index, centroids
                    

def make_prediction(cluster, sample_size=5):
    cluster = random.sample(cluster, sample_size)
    return majority_merge(multiple_alignment_muscle(cluster))
