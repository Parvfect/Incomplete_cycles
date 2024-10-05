
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from synthesis import NaiveSynthesisModel
from sequencing import NaiveSequencingModel


class NaiveDNAStorageModel:

    def __init__(self, coupling_rate, strand_length, repeats,
                  insertion_probatility=0.015, deletion_probability=0.055, subsitution_probability=0.075,
                  capping=True):
        self.synthesis_model = NaiveSynthesisModel(coupling_rate, strand_length, repeats, capping)
        self.strand = self.synthesis_model.strand
        self.sequencing_model = NaiveSequencingModel(insertion_probatility, deletion_probability, subsitution_probability, strand_length)
        self.strand_length = strand_length

    def simulate_storage(self):
        synthesized_strands = self.synthesis_model.simulate_synthesis()
        sequenced_strands = self.sequencing_model.simulate_full_sequencing(synthesized_strands)
        consensus_strand = self.sequencing_model.consensus_decoding(sequenced_strands)

        return consensus_strand
    

    



