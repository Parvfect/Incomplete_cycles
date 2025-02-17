
import Levenshtein
from tqdm import tqdm
from utils import reverse_complement
import numpy as np


def get_recovery_percentage(consensus_strand, original_strand):
    """Gets recovery percentage based on two strands. Chooses the length of the original strand to evaluate identity"""

    min_length = min(len(original_strand), len(consensus_strand))
    return sum([
                1 for i in range(min_length)
                if consensus_strand[i] == original_strand[i]]
                ) / len(original_strand)


def count_ids_errors(str1, str2):
    edit_operations = Levenshtein.editops(str1, str2)
    
    insertions = sum(1 for op in edit_operations if op[0] == 'insert')
    deletions = sum(1 for op in edit_operations if op[0] == 'delete')
    substitutions = sum(1 for op in edit_operations if op[0] == 'replace')

    return {'Insertions': insertions, 'Deletions': deletions, 'Substitutions': substitutions}

def evaluate_candidates(original_strands, candidates):

    oriented_candidates = []
    for candidate in candidates:

        rev = reverse_complement(candidate)

        if candidate in original_strands:
            oriented_candidates.append(candidate)
            continue

        if rev in original_strands:
            oriented_candidates.append(rev)
            continue

        flag = False
        best_ratio = 0.0
        for strand in original_strands:

            rec_straight = ratio(candidate, strand)
            rec_rev = ratio(rev, strand)

            if rec_straight > best_ratio:
                flag = False

            if rec_rev > best_ratio:
                flag = True

        if flag:
            oriented_candidates.append(rev)
        else:
            oriented_candidates.append(candidate)

    return oriented_candidates        