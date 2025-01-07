
import numpy as np
from aligned_clustering import conduct_align_clustering

def conduct_analysis(
        strand_id, coupling_rate, 
        capping, synthesized_strands, original_strand,
        deletions_per_strand, clustering = False, length_filtering = 0,
        running_on_hpc = False) -> dict:
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
            length_filtering=length_filtering,
            running_on_hpc=running_on_hpc)

    return analysis_dict