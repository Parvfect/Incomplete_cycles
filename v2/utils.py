

def read_synthesized_strands_from_file(file_path):

    sequences = []
    ids = []

    with open(file_path, 'r') as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith('>'):
            ids.append(line[1:].strip())
        if line.startswith('A') or line.startswith('C') or line.startswith('G') or line.startswith('T'):
            sequences.append(line.strip())

    return sequences, ids


def get_original_strands(original_strand_filepath):
    ids = []
    coupling_rates = []
    capping_flags = []
    strands = []

    with open(original_strand_filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                split_line = line.split()
                
                if len(split_line) > 1:
                    ids.append(split_line[0])
                    coupling_rates.append(split_line[1])
                    capping_flags.append(split_line[2])
                else:
                    strands.append(split_line[0])

    return ids, coupling_rates, capping_flags, strands