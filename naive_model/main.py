
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from utils import get_recovery_percentage

from dna_storage import NaiveDNAStorageModel


def recovery_percentage_coupling_rate(
        coupling_rates, strand_length, repeats, simulations, capping=True):
    recovery_percentages = []
    for coupling_rate in tqdm(coupling_rates):
        local_recovery_percentages = []
        for i in tqdm(range(simulations)):
            model = NaiveDNAStorageModel(coupling_rate, strand_length, repeats, capping=capping)
            consensus_strand = model.simulate_storage()
            recovery_percentage = get_recovery_percentage(consensus_strand, model.strand)
            local_recovery_percentages.append(recovery_percentage)

        recovery_percentages.append(np.mean(local_recovery_percentages))

    return recovery_percentages


def compare_capping_no_capping():
    
    strand_length = 300
    repeats = 1000
    simulations = 10

    # Assuming recovery percentages for different coupling rates are already calculated
    coupling_rates = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99]

    recovery_percentages_no_capping = recovery_percentage_coupling_rate(
        coupling_rates=coupling_rates,
        strand_length=strand_length, repeats=repeats,
        simulations=simulations, capping=False)

    recovery_percentages_with_capping = recovery_percentage_coupling_rate(
        coupling_rates=coupling_rates,
        strand_length=strand_length, repeats=repeats,
        simulations=simulations, capping=True)

    plt.figure(figsize=(10, 6))
    plt.plot(coupling_rates, recovery_percentages_no_capping,
             marker='o', linestyle='-', label='No Capping')
    plt.plot(
        coupling_rates, recovery_percentages_with_capping,
        marker='o', linestyle='-', label='With Capping')
    
    plt.xlabel('Coupling Rate', fontsize=12)
    plt.ylabel('Recovery Percentage', fontsize=12)
    plt.title(
    f'Recovery Percentage vs Coupling Rate with {repeats} repeats, for strand length {strand_length}', fontsize=14)

    plt.legend(fontsize=10)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.show()


if __name__ == '__main__':

    coupling_rate = 0.90
    strand_length = 300
    repeats = 1000
    capping = True
    model = NaiveDNAStorageModel(coupling_rate, strand_length, repeats, capping=capping, write_file=True)
    
    # Creating the synthesis file
    model.simulate_synthesis()

    # Run the Badread and read from file
    sequenced_strands = model.read_sequenced_strands_from_file('synthesized_strands.txt')

    # Run the consensus decoding
    consensus_strand = model.sequencing_model.consensus_decoding(sequenced_strands, model.strand, alignment=True)

    print(f"Percentage of correct bases: {sum([1 for i in range(len(consensus_strand)) if consensus_strand[i] == model.strand[i]]) / len(consensus_strand)}")

