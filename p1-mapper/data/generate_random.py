#!/usr/bin/env python3


import random
import os

def generate_reference(filename="reference.fasta", length=1_000_000):
    """
    Generates a random DNA sequence of the specified length (default 1 Mb) 
    and saves it to a FASTA file.
    """
    print(f"[1/2] Generating reference genome of {length:,} bases...")
    bases = ['A', 'T', 'C', 'G']
    # list comprehension and join for fast string concatenation
    sequence = "".join(random.choice(bases) for _ in range(length))
    
    # Note: now includes the 'data/' directory from the caller.
    with open(filename, "w") as f:
        f.write(">Ref_Chr1_Simulated_1Mbp\n")
        # Wrap sequence at 80 characters for standard FASTA format
        for i in range(0, length, 80):
            f.write(sequence[i:i + 80] + "\n")
    
    print(f"    - Success: Reference genome saved to {filename}")
    return sequence

def generate_reads(reference_sequence, ref_filename="reference.fasta", reads_filename="reads.fasta", num_reads=10000, read_length=150, error_rate=0.01):
    """
    Generates simulated short reads (Illumina-style) from a reference sequence.
    Reads include simulated sequencing errors (substitutions).
    
    Args:
        reference_sequence (str): The sequence from which to draw reads.
        reads_filename (str): The output file path for the reads.
        num_reads (int): The number of reads to generate.
        read_length (int): The length of each read (e.g., 150 bp).
        error_rate (float): The substitution error rate (e.g., 0.01 = 1%).
    """
    print(f"[2/2] Generating {num_reads:,} reads of length {read_length} with {error_rate*100:.1f}% error...")
    
    ref_len = len(reference_sequence)
    bases = ['A', 'T', 'C', 'G']
    
    if ref_len < read_length:
        print("    - Error: Reference genome is shorter than the requested read length. Aborting reads generation.")
        return

    def introduce_error(base):
        """Simulates a substitution error."""
        if random.random() < error_rate:
            # random substitution error
            return random.choice([b for b in bases if b != base])
        return base

    # Note: The file path now includes the 'data/' directory from the caller.
    with open(reads_filename, "w") as f:
        for i in range(num_reads):
            # 1. Choose a random start position
            start_pos = random.randint(0, ref_len - read_length)
            
            # 2. Extract the perfect sequence
            perfect_read = reference_sequence[start_pos:start_pos + read_length]
            
            # 3. Introduce errors (substitution)
            error_read = "".join(introduce_error(base) for base in perfect_read)
            
            # 4. Write the read to the file
            # Header format includes the original start position for validation
            f.write(f">Read_{i}_StartPos_{start_pos}\n")
            f.write(error_read + "\n")
            
    print(f"    - Success: Reads file saved to {reads_filename}")

def main():
    """Main function to orchestrate the data generation."""
    
    # --- Configuration ---
    DATA_DIR = "data"
    
    # Use os.path.join to correctly combine the directory and file name
    REF_FILENAME = os.path.join(DATA_DIR, "reference.fasta")
    READS_FILENAME = os.path.join(DATA_DIR, "reads.fasta")
    
    REFERENCE_LENGTH = 1_000_000  # 1 Mb
    NUM_READS = 10_000
    READ_LENGTH = 150
    ERROR_RATE = 0.01 # 1% error rate

    print("--- Starting Simulated Genome Data Generation ---")
    print("This script creates test data for a sequence mapper.")
    
    # --- Directory Creation ---
    # Create the output directory if it doesn't exist. exist_ok=True prevents errors 
    # if the directory is already there.
    os.makedirs(DATA_DIR, exist_ok=True)
    
    # 1. Generate the reference genome
    reference_sequence = generate_reference(REF_FILENAME, REFERENCE_LENGTH)
    
    # 2. Generate the reads from the sequence created in step 1
    if reference_sequence:
        generate_reads(reference_sequence, REF_FILENAME, READS_FILENAME, NUM_READS, READ_LENGTH, ERROR_RATE)
    
    print("\n--- Generation Complete ---")
    print(f"You can now run your mapper using '{REF_FILENAME}' as the reference and '{READS_FILENAME}' as the query reads. ")

if __name__ == "__main__":
    main()
