# Assembler and Mapper

**Algorithms for Sequencing Data - an implementation of a basic mapping and assembly algorithm for high-throughput sequencing data**

Younginn Park

## Mapper

Implements a basic read mapping algorithm that aligns high-throughput sequencing reads to a given reference genome, generating alignment information for each read.

Usage
>`python3 mapper.py reference.fasta reads.fasta output.txt`

## Assembler

Performs de novo assembly of high-throughput sequencing reads into longer contiguous sequences (contigs), reconstructing the genome without requiring a reference.

Usage
>`./assembly input_reads.fasta output_contigs.fasta`
