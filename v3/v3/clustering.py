
from typing import List, Tuple, Dict
from pool_preprocessing import remove_adapters_from_strands
from heirarchal_clustering import cluster_strands
from utils import reverse_complement
from strand_reconstruction import get_clustered_seqs, get_candidates, get_candidate_orientation
from evaluation import evaluate_candidates
import numpy as np


class Clustering:

    def __init__(
            self, strand_pool: List[str], reference_length: int, original_strands: List[str] = None,
            strand_pool_ids: List[str] = None, original_strand_ids: List[str] = None):

        self.strand_pool = strand_pool
        self.n_strands_pool = len(strand_pool)
        self.reference_length = reference_length
        self.original_strands = original_strands
        self.n_reference_strands = len(self.original_strands)
        self.strand_pool_ids = strand_pool_ids
        self.original_strands_ids = original_strand_ids

    def filter_by_length(
            self, max_length:int = 250, min_length: int = 190, ids: bool = False) -> List[str]:

        filtered_indices = [ind for ind in range(self.n_strands_pool) if len(
            self.strand_pool[ind]) < max_length and len(self.strand_pool[ind]) > min_length]
                
        strand_pool = [self.strand_pool[ind] for ind in filtered_indices]
        self.strand_pool = strand_pool
        self.n_strands_pool = len(self.strand_pool)
        
        if ids:
            strand_pool_ids = [self.strand_pool_ids[ind] for ind in filtered_indices]
            self.strand_pool_ids = strand_pool_ids

            return self.strand_pool, self.ids
        
        return self.strand_pool
    
    def remove_adapters(self, overwrite: bool = True) -> Tuple[List[str], List[str]]:

        strand_pool, strand_pool_ids = remove_adapters_from_strands(
            strands=self.strand_pool, original_strand_length=self.reference_length,
            ids=self.strand_pool_ids)
        
        if overwrite:
            self.strand_pool, self.strand_pool_ids = strand_pool, strand_pool_ids
            self.n_strands_pool = len(self.strand_pool)
            
        return strand_pool, strand_pool_ids

    def cluster_strand_pool(self, distance_threshold:int = 40, strand_pool: List[str] = None):

        if strand_pool:
            cluster_dict = cluster_strands(
                strand_pool=strand_pool, distance_threshold=distance_threshold)
        else:
            cluster_dict = cluster_strands(
                strand_pool=self.strand_pool, distance_threshold=distance_threshold)            

        self.clusters = self.cluster_dict['clusters']
        self.reversed_markers = self.cluster_dict['reversed_markers']
        self.cluster_heads = self.cluster_dict['cluster_heads']

        assert [len(self.clusters[i]) > len(self.clusters[i + 1]) for i in range(len(self.clusters) - 1)]
        print("Clusters are sorted")
        
        self.clustered_seqs = get_clustered_seqs(
            clusters=self.clusters, reversed_markers=self.reversed_markers, strand_pool=self.strand_pool)
        print("Orientation fixed in the strand pool")

        return self.clustered_seqs
    
    def generate_candidates(
            self, n_candidates: int, n_samples: int = 15, clustered_seqs: List[List[str]] = None,
            fix_orientation=True) -> List[str]:

        clustered_seqs = clustered_seqs[:n_candidates]

        if clustered_seqs:
            self.candidates = get_candidates(clustered_seqs=clustered_seqs, n_samples=n_samples)
        else:
            self.candidates = get_candidates(clustered_seqs=self.clustered_seqs, n_samples=n_samples)

        if fix_orientation:
            print("Fixing candidate orientations")
            reversed_markers = get_candidate_orientation(original_strands=self.original_strands, candidates=candidates)
            n_reversed = sum(reversed_markers)
            print(f"{n_reversed/len(n_candidates)} candidates are reversed")
            candidates = [
                reverse_complement(candidates[ind]) if reversed_markers[ind] else candidates[ind] for ind in range(len(self.candidates))]
            self.candidates = candidates
            
        return self.candidates
    
    def evaluate_candidates(self, candidates: List[str] = None) -> Dict[str, np.ndarray]:
        
        if candidates:
            self.evaluation_dict = evaluate_candidates(
                original_strands=self.original_strands,
                candidates=candidates
            )
        else:
            self.evaluation_dict = evaluate_candidates(
                original_strands=self.original_strands,
                candidates=self.candidates
            )

        self.reference_recoveries = self.evaluation_dict['reference_recoveries']
        self.reference_strand_indices = self.evaluation_dict['reference_strand_indices']
        self.recovery_rates = self.evaluation_dict['recovery_rates']

        return self.evaluation_dict
    
    def fsm(self, candidates: List[str] = None, hard = False) -> bool:
        assert len(candidates) == self.n_reference_strands
        
        found = 0
        for i in self.original_strands:
            if i in candidates:
                found += 1
                continue
        
        if found == self.n:
            print("Found all")
            return True

        if hard:
            return False

        else:
            print(f"Found {found}")
            return False