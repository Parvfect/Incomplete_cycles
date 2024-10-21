



def get_original_strands():
    ids = []
    coupling_rates = []
    capping_flags = []
    strands = []

    with open('original_strands.txt', 'r') as f:
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

