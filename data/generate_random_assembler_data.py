#!/usr/bin/env python3
import random
import os

def generate_reference(length=50_000, filename="data/reference.fasta"):
    """Generates a random DNA genome and writes to a FASTA file."""
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    bases = ['A', 'T', 'C', 'G']
    sequence = "".join(random.choice(bases) for _ in range(length))
    with open(filename, "w") as f:
        f.write(">RandomGenome\n")
        for i in range(0, length, 80):
            f.write(sequence[i:i+80] + "\n")
    return sequence

def generate_reads(reference, reads_filename="data/input_reads.fasta", num_reads=1000, read_length=150, error_rate=0.01):
    """Generates reads from a reference genome with optional substitution errors."""
    os.makedirs(os.path.dirname(reads_filename), exist_ok=True)
    bases = ['A', 'T', 'C', 'G']

    def mutate(base):
        return random.choice([b for b in bases if b != base]) if random.random() < error_rate else base

    with open(reads_filename, "w") as f:
        for i in range(num_reads):
            start = random.randint(0, len(reference) - read_length)
            read_seq = "".join(mutate(b) for b in reference[start:start+read_length])
            f.write(f">Read_{i}_pos{start}\n{read_seq}\n")

    print(f"Generated {num_reads} reads from a {len(reference)}-bp genome in {reads_filename}")

if __name__ == "__main__":
    REF_FILE = "data/reference.fasta"
    READS_FILE = "data/input_reads.fasta"
    REFERENCE_LENGTH = 50_000   # smaller genome for quick testing
    NUM_READS = 2000
    READ_LENGTH = 150
    ERROR_RATE = 0.01

    genome = generate_reference(REFERENCE_LENGTH, REF_FILE)
    generate_reads(genome, READS_FILE, NUM_READS, READ_LENGTH, ERROR_RATE)
