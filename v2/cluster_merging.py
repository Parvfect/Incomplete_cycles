
import operator

# This function returns the fraction of origignal squences recovered given a number of candidates
def fraction_recovered(candidates,orig_seqs):
    d = {}
    for seq in orig_seqs:
        d[seq] = 0
    for cand in candidates:
        if cand in d:
            d[cand] += 1
    av = sum([ d[seq]>0 for seq in d]) / len(d)
    print("Fraction of recovered sequences: ", av )
    if av>0:
        print("Fraction of recovered sequences: ", sum([ d[seq] for seq in d]) / len(d) / av )

def majority_merge(reads,weight = 0.4):
    # assume reads have the same length
    res = ""
    for i in range(len(reads[0])):
        counts = {'A':0,'C':0,'G':0,'T':0,'-':0,'N':0}
        for j in range(len(reads)):
            if i >= len(reads[j]):
                continue
            counts[reads[j][i]] +=1
        counts['-'] *= weight
        mv = max(counts.items(), key=operator.itemgetter(1))[0]
        if mv != '-':
            res += mv
    return res

def merge_clusters(fresults):
    candidates = []
    for ma in fresults:
        candidates.append(majority_merge(ma,weight=0.5))
    
    return candidates