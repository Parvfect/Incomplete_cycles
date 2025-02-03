


import time
import numpy as np
import itertools
from skbio.alignment import local_pairwise_align_ssw
from skbio import DNA
from scipy.cluster.hierarchy import linkage, fcluster
from collections import defaultdict
from tqdm import tqdm
import uuid


#===== assign numbers to shingles of each sequence=====#
def kmerDNA(seq,k=3):
    kmer = []
    for ell in range(len(seq)-k+1):
        nstr = seq[ell:ell+k]
        index = 0
        for j,c in enumerate(nstr):
            if c == 'A':
                i = 0
            elif c == 'C':
                i = 1
            elif c == 'G':
                i = 2
            elif c == 'T':
                i = 3
            else:
                index = -1
                break
            index += i*(4**j)
        kmer += [index]
    return kmer

#=====min-hash object=====#
class minhashsig():
    # min-hash of k-mers
    def __init__(self,m,k):
        # m is the number of signatures
        self.tables = [np.random.permutation(4**k) for i in range(m)]
        self.k = k
    def generate_signature(self,seq):
        kmer = kmerDNA(seq,self.k)
        sig = [ min([table[i] for i in kmer]) for table in self.tables]
        return sig
    
#=====pair detection=====#
def extract_similar_pairs(sigs,m,k_lsh,ell_lsh,maxsig):
    # sigs: minhash signatures
    # ell_lsh: number of LSH signatures
    # k_lsh: number of MH signatures to be concatenated
    # we use generatrs to yield a number of pairs at a time for the sake of memory efficiency
    
    pairs = set([])
    
    # generate ell_lsh random indices
    for ell in range(ell_lsh):
        pair_count = 0
        s = time.time()
        lshinds = np.random.permutation(m)[:k_lsh]
        # generate LSh signatures
        lshsigs = []
        for sig in sigs:
            lshsig = 0
            for i,lshind in enumerate(lshinds):
                lshsig += sig[lshind]*(maxsig**i)
            lshsigs += [lshsig]
        d = {}
        for ind,sig in enumerate(lshsigs):
            if sig in d:
                d[sig] += [ind]
            else:
                d[sig] = [ind]
        for candidates in d.values():
            cent = set([])
            if len(candidates) > 1:
                for pair in itertools.combinations(candidates,2):
                    cent.add(pair[0])
                    if len(cent)==1:
                        pairs.add(pair)
                    else:
                        break
                        
        yield pairs,ell
        pair_count += len(pairs)
        pairs = set([])

#=====form clusters based on pairs=====#
def center_cluster(pairs):
    clusters = {}
    hold = 0
    t_counter = 0
    ell_copy = 0
    pairsize = 0
    while not hold:
        
        try:
            out = next(pairs)
            pairs_sort = list(out[0])
            ell = out[1]
            pairsize += len(pairs_sort)
            pairs_sort.sort()
            s = time.time()
            for (u,v) in pairs_sort:
                if u in clusters:
                    clusters[u] += [v]
        
                if v in clusters:
                    clusters[v] += [u]
        
                if v not in clusters and u not in clusters:
                    clusters[u] = [v]
    
        except StopIteration:
            hold = 1
            #print("clustering completed","---",pairsize,"pairs clustered")
        if ell==ell_copy:
            t_counter += time.time()-s
        else:
            #print("Clustering time for LSH",ell_copy,":",t_counter,'\n')
            t_counter = time.time()-s
            ell_copy = ell

    return clusters  # Returns a list of limited clusters
 


#=====LSH clustering (main function)=====#
def lsh_cluster(seqs,m,k,k_lsh=2,ell_lsh=4):
    # This is the main function
    maxsig = 4**k
    minhash = minhashsig(m,k)
    sigs = [minhash.generate_signature(seq[:190]) for seq in seqs]
    
    pairs = extract_similar_pairs(sigs ,m, k_lsh, ell_lsh, maxsig)
    #clusters = center_cluster(pairs)
    clusters = center_cluster(pairs=pairs)

    return clusters


def filter_nonunique(seqs):
    # filter out sequences that appear many times
    d = {}
    ctr = 0
    for i,seq in enumerate(seqs):
        if seq in d:
            d[seq] += [i]
        else:
            d[seq] = [i]
    import operator
    sorted_d = sorted(d.items(), key=operator.itemgetter(1))
    sorted_d.reverse()
    return d,ctr


def filter_lsh_clusters(clusters, reads):
    clusts = [ [c] + list(set(clusters[c])) for c in clusters if len(clusters[c]) > 3 ]

    #=====max matching=====# 
    def max_match(seq1, seq2):
        # This function checks whether seq1 and seq2 are similar or not
        # Checking all pairs within a cluster dramatically increases the time complexity, 
        # so by default, in the next cell, we call this function to only check the pairs
        # that one of their members is the cluster center
        
        alignment, score, start_end_positions \
            = local_pairwise_align_ssw(DNA(seq1) , DNA(seq2) , match_score=2,mismatch_score=-3)
        a = str(alignment[0])
        b = str(alignment[1])
        ctr = 0
        for i,j in zip(a,b):
            if i==j:
                ctr += 1
        return ctr
    
    th = 10 # filtering threshold

    k = len(clusts)
    s = time.time()
    fclusts = []
    for i,c in enumerate(clusts):
        cent = reads[c[0]]
        cc = [c[0]]
        for e in c[1:]:
            score = max_match(cent,reads[e])
            if score >= th:
                cc += [e]
        fclusts += [cc]
        if i%1000 == 0:
            print("%",round(i*100/len(clusts),2),"of the clusters are filtered.")
    print("filtering time for",k,"clusters:",round(time.time()-s,2),"s")

    return fclusts

def create_clusters(trimmed_seqs, TRIVIAL=False, nbeg=14, target_clusters=None, m=10):

    if TRIVIAL:
        start = time.time()
        d,ctr = filter_nonunique([seq[:nbeg] for seq in trimmed_seqs])
        clusters = [d[a] for a in d if len(d[a]) > 3]
        end = time.time()
        fclusts = clusters.copy()
    else:
        # set up the parameters and call the lsh_cluster function
        k_lsh = 4
        sim = 0.5
        ell_lsh = int(1 / (sim ** k_lsh))
        k = 5
        start = time.time()
        clusters = lsh_cluster(trimmed_seqs, m , k, k_lsh, ell_lsh)
        print(f"Initial LSH clusters: {len(clusters)}")

        end = time.time()

        print("Runtime:",round(end-start,1),"s")
        print(len(clusters),"number of clusters created")

        fclusts = filter_lsh_clusters(clusters=clusters, reads=trimmed_seqs)
        #print(fclusts)

    
    cluster_ids = [str(uuid.uuid4()) for i in range(len(fclusts))]
    return fclusts, cluster_ids

"""
Trivial clustering
1. Looks at first 14 bases of the strand and puts unique one into a new cluster
This is really really not good. But still does remarkably well. We need to try LSH Hashing.
"""