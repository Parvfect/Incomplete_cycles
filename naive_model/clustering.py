


import time
import numpy as np
import itertools

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
            print("clustering completed","---",pairsize,"pairs clustered")
        if ell==ell_copy:
            t_counter += time.time()-s
        else:
            print("Clustering time for LSH",ell_copy,":",t_counter,'\n')
            t_counter = time.time()-s
            ell_copy = ell
 
    return clusters

#=====LSH clustering (main function)=====#
def lsh_cluster(seqs,m,k,k_lsh=2,ell_lsh=4):
    # This is the main function
    maxsig = 4**k
    minhash = minhashsig(m,k)
    sigs = [minhash.generate_signature(seq[:40]) for seq in seqs]
    pairs = extract_similar_pairs(sigs,m,k_lsh,ell_lsh,maxsig)
    clusters = center_cluster(pairs)
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


def create_clusters(trimmed_seqs, TRIVIAL=True):
    TRIVIAL = True # when real synthesized data is selected, set this flag to "True" (this dramatically reduces the runtime)

    if TRIVIAL:
        start = time.time()
        nbeg = 14
        d,ctr = filter_nonunique([seq[:nbeg] for seq in trimmed_seqs])
        clusters = [d[a] for a in d if len(d[a]) > 3]
        end = time.time()
        print("Runtime:",round(end-start,1),"s")
        print(len(clusters),"number of clusters created.")
        fclusts = clusters.copy()
    else:
        # set up the parameters and call the lsh_cluster function
        k_lsh = 4
        sim = 0.5
        ell_lsh = int(1/(sim**k_lsh))
        m,k=50,5
        start = time.time()
        clusters = lsh_cluster(trimmed_seqs,m,k,k_lsh,ell_lsh)
        end = time.time()

        print("Runtime:",round(end-start,1),"s")
        print(len(clusters),"number of clusters created")


    return clusters