

from synthesis import NaiveSynthesisModel
from sequencing import NaiveSequencingModel
from utils import read_strands_from_file


class NaiveDNAStorageModel:

    def __init__(
            self, coupling_rate, strand_length, repeats,
            insertion_probatility=0.04,
            deletion_probability=0.1,
            subsitution_probability=0.15, capping=True, write_file=False):
        self.synthesis_model = NaiveSynthesisModel(
            coupling_rate, strand_length, repeats, capping,
            write_file=write_file)       
        self.strand = self.synthesis_model.strand
        self.repeats = repeats
        self.sequencing_model = NaiveSequencingModel(
            insertion_probatility, deletion_probability,
            subsitution_probability, strand_length)
        self.strand_length = strand_length
        self.write_file = write_file

    def simulate_synthesis(self):
        self.synthesized_strands = self.synthesis_model.simulate_synthesis()
        return self.synthesized_strands

    def simulate_storage(self):
        self.synthesized_strands = self.synthesis_model.simulate_synthesis()
        self.sequenced_strands = self.sequencing_model.simulate_full_sequencing(
            self.synthesized_strands)
        self.consensus_strand = self.sequencing_model.consensus_decoding(
            self.sequenced_strands, self.strand, alignment=True)

        return self.consensus_strand
    
    def read_sequenced_strands_from_file(self, file_path):
        ids, sequences = read_strands_from_file(file_path)
        self.sequenced_strands = sequences
        return sequences
