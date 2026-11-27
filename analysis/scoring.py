
def scan_sequences(sequences, target, max_mismatches=3, step=1):

    results = []

    for record in sequences:

        # Get sequence ID and DNA sequence
        if isinstance(record, dict):
            sequence_id = record.get("id", "unknown")
            dna_sequence = record.get("sequence", "")
        else:
            sequence_id = "unknown"
            dna_sequence = str(record)

        dna_sequence = dna_sequence.upper()

        # Scan DNA sequence
        for i in range(0, len(dna_sequence) - len(target) + 1, step):

            chunk = dna_sequence[i:i + len(target)]

            # Count mismatches
            mismatches = sum(
                1 for a, b in zip(chunk, target)
                if a != b
            )

            if mismatches <= max_mismatches:
                results.append({
                    "id": sequence_id,
                    "position": i,
                    "sequence": chunk,
                    "target": target,
                    "mismatches": mismatches
                })

    return results
