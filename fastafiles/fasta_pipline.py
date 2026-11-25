
# Import Path to validate and manage FASTA file paths.
from pathlib import Path

# Import StringIO to read FASTA text directly from memory.
from io import StringIO

# Import NumPy for efficient numerical representation of DNA bases.
import numpy as np

# Import SeqIO from BioPython to parse FASTA format files.
from Bio import SeqIO


# Dictionary that maps each DNA nucleotide to a numerical value.
# A, C, G, T represent the four standard DNA bases.
# N represents an unknown or undetermined DNA base.
DNA_ENCODING = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "N": -1,
}


# --------------------------------------------------
# FUNCTION 1: CLEAN AND VALIDATE DNA SEQUENCE
# --------------------------------------------------
def clean_dna(sequence):

    # Remove unnecessary spaces, tabs and newline characters.
    # Convert lowercase DNA characters into uppercase.
    sequence = "".join(sequence.split()).upper()

    # Identify DNA characters that are not present in our encoding map.
    invalid = set(sequence) - set(DNA_ENCODING)

    # Raise ValueError if the sequence contains invalid DNA characters.
    if invalid:
        raise ValueError(f"Invalid DNA bases: {sorted(invalid)}")

    # Return the validated and cleaned DNA sequence.
    return sequence


# --------------------------------------------------
# FUNCTION 2: ENCODE DNA SEQUENCE INTO NUMBERS
# --------------------------------------------------
def encode_dna(sequence):

    # First clean and validate the input DNA sequence.
    sequence = clean_dna(sequence)

    # Convert each DNA nucleotide into its corresponding numerical value.
    # NumPy int8 stores the encoded values efficiently using 8-bit integers.
    return np.array(
        [DNA_ENCODING[base] for base in sequence],
        dtype=np.int8,
    )


# --------------------------------------------------
# FUNCTION 3: READ FASTA FILE AND PREPARE DNA CHUNKS
# --------------------------------------------------
def prepare_data(file_path, chunk_size=5, step=1):

    # Validate chunk_size and step to ensure both are positive.
    if chunk_size <= 0 or step <= 0:
        raise ValueError("chunk_size and step must be positive")

    # Convert the given file path into a Path object.
    path = Path(file_path)

    # Verify that the specified FASTA file exists.
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    # Open the FASTA file using UTF-8 encoding.
    # utf-8-sig handles files containing a UTF-8 byte-order mark.
    with path.open("r", encoding="utf-8-sig") as file:

        # Read non-empty lines and remove unnecessary whitespace.
        lines = [
            line.strip()
            for line in file
            if line.strip()
        ]

    # Locate the first FASTA header.
    # FASTA headers always begin with the '>' symbol.
    # Return None if no header is available.
    first_header = next(
        (i for i, line in enumerate(lines)
         if line.startswith(">")),
        None,
    )

    # Validate the FASTA structure before parsing.
    if first_header is None:
        raise ValueError(
            "No FASTA header found. Each sequence must start with >."
        )

    # Remove any unwanted text appearing before the first FASTA header.
    lines = lines[first_header:]

    # Combine the FASTA lines into a single text string.
    fasta_text = "\n".join(lines) + "\n"

    # Parse the FASTA content using BioPython's SeqIO module.
    # StringIO makes the in-memory string behave like a file.
    # list() collects all parsed FASTA records.
    fasta_records = list(
        SeqIO.parse(StringIO(fasta_text), "fasta")
    )

    # Check whether at least one valid FASTA record was found.
    if not fasta_records:
        raise ValueError("No valid FASTA sequences found.")

    # Initialize an empty list to store the processed DNA chunks.
    records = []

    # Iterate through every FASTA record in the input file.
    for fasta_record in fasta_records:

        # Extract the DNA sequence from the current FASTA record.
        # Convert it into a string and perform cleaning and validation.
        dna = clean_dna(str(fasta_record.seq))

        # Reject empty DNA sequences to avoid invalid processing.
        if not dna:
            raise ValueError(f"Empty DNA sequence: {fasta_record.id}")

        # Convert the complete DNA sequence into its numerical representation.
        encoded = encode_dna(dna)

        # Generate starting positions for DNA chunk creation.
        # If the sequence is shorter than chunk_size, create one short chunk.
        # Otherwise, generate overlapping chunks based on the given step.
        starts = (
            [0]
            if len(dna) < chunk_size
            else range(0, len(dna) - chunk_size + 1, step)
        )

        # Iterate through each starting position to extract DNA chunks.
        for start in starts:

            # Calculate the ending position without exceeding sequence length.
            end = min(start + chunk_size, len(dna))

            # Store the processed chunk information as a dictionary.
            records.append({

                # Store the unique identifier of the FASTA sequence.
                "id": fasta_record.id,

                # Store the starting position using 1-based indexing.
                "position": start + 1,

                # Extract the original DNA nucleotide sequence for this chunk.
                "sequence": dna[start:end],

                # Extract the corresponding numerical DNA representation.
                # copy() creates an independent array for each chunk.
                "encoded": encoded[start:end].copy(),

                # Store the actual number of DNA bases in this chunk.
                "length": end - start,
            })

    # Return all processed DNA chunks for further analysis.
    return records
