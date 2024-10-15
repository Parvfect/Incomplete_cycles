
import numpy as np
from uuid import uuid4

bases = ['A', 'T', 'C', 'G']


class NaiveSynthesisModel:
    """
    Single strand synthesis model -
    n length strand, k repeats, with coupling rate alpha per base
    Option for capping or not, to end synthesis if a base is not added
    """

    def __init__(self, coupling_rate, strand_length, repeats,
                 capping=True, write_file=False):
        self.coupling_rate = coupling_rate
        self.strand_length = strand_length
        self.repeats = repeats
        self.strand = "".join(np.random.choice(bases, strand_length))
        self.capping = capping
        self.write_file = write_file

    def simulate_synthesis(self):

        self.synthesized_strands = []
        for i in range(self.repeats):
            synthesizing_strand = ""
            for j in range(self.strand_length):
                if np.random.rand() < self.coupling_rate:
                    synthesizing_strand += self.strand[j]
                else:
                    if self.capping:
                        break
                    continue
            self.synthesized_strands.append(synthesizing_strand)

        if self.write_file:
            with open('synthesized_strands.txt', 'w') as f:
                for strand in self.synthesized_strands:
                    f.write(f">strand a id {uuid4().hex}\n")
                    f.write(strand + '\n\n')

        return self.synthesized_strands
