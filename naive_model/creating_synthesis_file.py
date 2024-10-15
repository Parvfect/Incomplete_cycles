
from synthesis import NaiveSynthesisModel

coupling_rates = [8, 8.5, 9, 9.5, 9.99]

synthesis_models = []

repeats = 1


# Starting a new file
with open('original_strands.txt', 'w') as f:
    f.write("\n")

with open('synthesized_strands.txt', 'w') as f:
    f.write("\n")

filepath = "synthesized_strands.txt"

# Creating all the synthesis models
for coupling_rate in coupling_rates:
    for _ in range(repeats):
        synthesis_models.append(NaiveSynthesisModel(
            coupling_rate, strand_length=250, repeats=1500, capping=True, write_file=True))

        synthesis_models.append(NaiveSynthesisModel(
            coupling_rate, strand_length=250, repeats=1500, capping=False, write_file=True))
 

# Get all the original strands and write them to a file
for model in synthesis_models:

    with open('original_strands.txt', 'a') as f:
        f.write(
            f'{model.strand_id} {model.coupling_rate} {model.capping}\n{model.strand}\n\n')

# Synthesise strands and write them
for model in synthesis_models:
    model.simulate_synthesis(f"{filepath}")