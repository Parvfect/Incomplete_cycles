
from synthesis import NaiveSynthesisModel
from datetime import datetime
from synthesis_analysis import conduct_analysis
import os
from tqdm import tqdm
import cProfile
import re

synthesis_models = []
running_on_hpc = True

# Parameters (argparser eventually, okay for now)
coupling_rates = [0.99]
sim_repeats_per_coupling_rate = 1
strand_repeats = 1000
strand_length = 200

# Initiating sim run data path
timestamp = str(datetime.now()).replace(':', '.')
preceeding_path = os.path.join('runs', timestamp)

if running_on_hpc:
    home_dir = os.environ['HOME']
    preceeding_path = os.path.join(home_dir, preceeding_path)

# Assuming that the runs folder already exists in the directory - otherwise it breaks
os.mkdir(preceeding_path)

synthesized_strands_write_path = os.path.join(preceeding_path, 'synthesized.fasta')
original_strand_write_path = os.path.join(preceeding_path, 'original_strands.txt')
parameters_path = os.path.join(preceeding_path, 'run_info_file.txt')

# Starting a new file
with open(original_strand_write_path, 'w') as f:
    f.write("")

with open(synthesized_strands_write_path, 'w') as f:
    f.write("")

with open(parameters_path, 'w') as f:
    f.write(f"\nRun on {timestamp}\n")
    f.write(f"Coupling Rates = {coupling_rates}\n")
    f.write(f"Simulation Repeats = {sim_repeats_per_coupling_rate}\n")
    f.write(f"Strand repeats = {strand_repeats}\n")
    f.write(f"Strand length = {strand_length}\n")
    

# Creating all the synthesis models
for coupling_rate in coupling_rates:
    for _ in range(sim_repeats_per_coupling_rate):
        synthesis_models.append(NaiveSynthesisModel(
            coupling_rate, strand_length=strand_length, repeats=strand_repeats, capping=True, write_file=False))

        synthesis_models.append(NaiveSynthesisModel(
            coupling_rate, strand_length=strand_length, repeats=strand_repeats, capping=False, write_file=False))
        
"""
# Get all the original strands and write them to the file
for model in synthesis_models:
    with open(original_strand_write_path, 'a') as f:
        f.write(
            f'{model.strand_id} {model.coupling_rate} {model.capping}\n{model.strand}\n\n')
"""
                 
# Synthesise strands and write them - add analysis?
synthesized_strands_arr = []

for model in tqdm(synthesis_models):

    synthesized_strands, strand_deletions = model.simulate_synthesis(return_deletions=True)
    strand_id = str(model.strand_id)

    strand_analysis_dict = (conduct_analysis(
        strand_id=strand_id, coupling_rate=model.coupling_rate,
        capping=model.capping,
        synthesized_strands=synthesized_strands,
        deletions_per_strand=strand_deletions,
        original_strand=model.strand,
        clustering=True,
        length_filtering=0.2,
        running_on_hpc=running_on_hpc
    ))

    print(strand_analysis_dict)

    with open(parameters_path, 'a') as f:
        f.write(str(strand_analysis_dict))
        f.write('\n')
    
### Write to file
"""
# So one file for each seperate model - about 20ish files
split_strands = [synthesized_strands[i:i + 9000] for i in range(0, len(synthesized_strands) - 9001, 9000)]

for i, strands_ in enumerate(split_strands):  
    write_path = os.path.join(preceeding_path, f'synthesized_{strand_id}_{i}.fasta')
    with open(write_path, 'w') as f:
        for strand in strands_:

            if len(strand) < 100: # PBSim does not accept strands that are less than 100 bases long
                continue
            f.write(f">{strand_id}\n")
            f.write(strand + '\n\n')
"""