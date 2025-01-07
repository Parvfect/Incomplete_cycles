
import numpy as np
from aligned_clustering import conduct_align_clustering

def conduct_analysis(
        strand_id: str, coupling_rate: float, 
        capping: bool, synthesized_strands: list[str], original_strand: str,
        deletions_per_strand: list[int], clustering: bool = False, length_filtering:float = 0) -> dict:
    """
    Takes in the synthesized strands and the original strand. Gives out the mean length, std, max length, deletions (got to change the synthesis function for this) for the synthesized strands. Clustering to see the best candidate recovery percentage after simulating synthesis.
    """

    synthesized_strand_lengths = [len(i) for i in synthesized_strands]
    
    analysis_dict = {
        "type": "capping" if capping else "no_capping",
        "coupling_rate": coupling_rate,
        "best_recovery_clustering": 0.0,
        "mean_length": np.mean(synthesized_strand_lengths),
        "max_length": max(synthesized_strand_lengths),
        "std_length": np.std(synthesized_strand_lengths),
        "mean_deletions": np.mean(deletions_per_strand),
        "std_deletions": np.std(deletions_per_strand),
        "strand_id": strand_id,
        "strand": original_strand
    } # Remember, recovery could change to Levenshtien distance. But let's do that later.
    
    if clustering:
        analysis_dict['best_recovery_clustering'] = conduct_align_clustering(
            original_strand=original_strand,
            trimmed_seqs=synthesized_strands,
            display=False, best_recovery=True,
            length_filtering=length_filtering)

    return analysis_dict