
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
coupling_rates = [0.7, 0.75, 0.8, 0.85, 0.9, 0.925, 0.95, 0.975, 0.99]
sim_repeats_per_coupling_rate = 3
strand_repeats = 10000
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
       
# Get all the original strands and write them to the file
for model in synthesis_models:
    with open(original_strand_write_path, 'a') as f:
        f.write(
            f'{model.strand_id} {model.coupling_rate} {model.capping}\n{model.strand}\n\n')

                 
# Synthesise strands and write them - add analysis?
synthesized_strands_arr = []

for model in tqdm(synthesis_models):

    synthesized_strands, strand_deletions = model.simulate_synthesis(return_deletions=True)
    strand_id = str(model.strand_id)

    try:
        strand_analysis_dict = (conduct_analysis(
            strand_id=strand_id, coupling_rate=model.coupling_rate,
            capping=model.capping,
            synthesized_strands=synthesized_strands,
            deletions_per_strand=strand_deletions,
            original_strand=model.strand,
            clustering=True,
            length_filtering=0,
            running_on_hpc=running_on_hpc
        ))

        with open(parameters_path, 'a') as f:
            f.write(str(strand_analysis_dict))
            f.write('\n')

        with open(synthesized_strands_write_path, 'a') as f:
            for strand in synthesized_strands:
                if strand and len(strand) > 1:
                    f.write(f">{strand_id}\n")
                    f.write(strand + '\n\n')
    
    except Exception as e:
        continue
