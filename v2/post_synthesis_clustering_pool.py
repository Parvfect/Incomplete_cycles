
from aligned_clustering import conduct_align_clustering

# Read original strands from the file
original_strands = []

# Read synthesised strands from file - 360,000 of these
synthesized_strands = []

# sample 1/100th and see how much time it takes and increase the number in the pool in that fashion. Ideally, this is where I would like to start profiling the clustering code
# If I can make this fast, I can do a lot of stuff in conjunction remember, this is a key step - need to get it right.

recoveries = conduct_align_clustering(
    original_strand=original_strands,
    trimmed_seqs=synthesized_strands,
    display=False,
    multiple=True
)