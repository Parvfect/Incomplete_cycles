

# reduce to kmers of 3
# get hash signatures
# normalise and clean
# kmeans it to the required number of clusters
# since there are no adapters, should be easier right

from clustering import kmerDNA, minhashsig
from tqdm import tqdm
from sklearn.cluster import KMeans
import numpy as np
from sklearn.preprocessing import normalize, StandardScaler


def normalise_frequency_counts(frequency_counts):
    scaler = StandardScaler()
    return scaler.fit_transform(frequency_counts)

def get_kmer_frequency_counts(seqs, k=3):
    freq_counts = np.zeros([len(seqs), 4**k])

    for i, seq in tqdm(enumerate(seqs)):
        kmers = kmerDNA(seq)
        for kmer in kmers:
            freq_counts[i, kmer] += 1
    
    return freq_counts
    
def get_hash_signatures(seqs, m=8, k=3, nbeg=40):

    hash_tables = minhashsig(m=m, k=k)

    hash_signatures = []

    for i in tqdm(seqs):
        hash_signatures.append(
            hash_tables.generate_signature(kmerDNA(i[:nbeg])))
        

    return hash_signatures

def cluster(hash_signatures, n_clusters):

    kmeans = KMeans(n_clusters=n_clusters).fit(hash_signatures)
    return kmeans
