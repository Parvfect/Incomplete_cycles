
from Levenshtein import ratio


def filter_junk_reads(records, ids=None, similarity_threshold=0.85):
    """
    Removes all the sequences that are not similar to any others.
    Prints out percentage of sequences removed
    """
    filtered_records = []
    filtered_seqs = set()

    for record in records:
        for record_ in records:
            
            # Checking if its the same record
            if record is record_:
                continue
            
            # Checking if the record is already in the filtered pool
            if not record_.seq in filtered_seqs:
                if ratio(record.seq, record_.seq) > similarity_threshold:
                    filtered_records.append(record_)
                    filtered_seqs.add(record_.seq)

    print(f"{100 - (len(filtered_records) * 100 / len(records))} percent sequences filtered out")

    return filtered_records
                    

