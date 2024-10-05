
import numpy as np
from dna_storage import NaiveDNAStorageModel

coupling_rate = 0.90
strand_length = 200
repeats = 1000
simulations = 20
capping=False

model = NaiveDNAStorageModel(coupling_rate, strand_length, repeats, capping=capping)
synthesized_strands = model.synthesis_model.simulate_synthesis()
sequenced_strands = model.sequencing_model.simulate_full_sequencing(synthesized_strands)

sequencing_model = model.sequencing_model

A_votes, T_votes, C_votes, G_votes = sequencing_model.get_base_votes(sequenced_strands)
likelihoods = sequencing_model.get_base_likelihoods(A_votes, T_votes, C_votes, G_votes)

print(likelihoods[:10])