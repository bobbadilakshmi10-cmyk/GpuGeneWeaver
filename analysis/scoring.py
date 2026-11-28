# Function to scan DNA sequences and identify matching regions
# based on the allowed number of mismatches.
def scan_sequences(sequences, target, max_mismatches=3, step=1):

    # Create an empty list to store the matching DNA regions.
    results = []

    # Iterate through each DNA sequence record.
    for record in sequences:

        # Check whether the record is stored as a dictionary.
        if isinstance(record, dict):

            # Extract the sequence ID; use "unknown" if the ID is missing.
            sequence_id = record.get("id", "unknown")

            # Extract the DNA sequence; use an empty string if missing.
            dna_sequence = record.get("sequence", "")

        else:
            # If the record is not a dictionary, convert it into a string.
            sequence_id = "unknown"
            dna_sequence = str(record)

        # Convert the DNA sequence to uppercase for consistent comparison.
        dna_sequence = dna_sequence.upper()

        # Slide a window across the DNA sequence to extract target-sized chunks.
        # The step value determines how many positions the window moves each time.
        for i in range(0, len(dna_sequence) - len(target) + 1, step):

            # Extract a DNA chunk having the same length as the target sequence.
            chunk = dna_sequence[i:i + len(target)]

            # Count the number of mismatched bases between the chunk and target.
            mismatches = sum(
                1 for a, b in zip(chunk, target)
                if a != b
            )

            # Check whether the mismatch count is within the allowed limit.
            if mismatches <= max_mismatches:

                # Store the matching region and its details in the results list.
                results.append({
                    "id": sequence_id,       # Original DNA sequence ID
                    "position": i,            # Starting position of the match
                    "sequence": chunk,        # Extracted DNA chunk
                    "target": target,         # Target DNA sequence
                    "mismatches": mismatches  # Number of mismatched bases
                })

    # Return all DNA regions that satisfy the mismatch condition.
    return results