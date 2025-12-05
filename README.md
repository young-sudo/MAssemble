# htskit-mini

**Algorithms for High-Throughput Sequencing Data - a lightweight implementation of a basic mapping and assembly algorithms from scratch in Python**

Younginn Park

![Python](https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white)
![Biopython](https://img.shields.io/badge/Biopython-0096D8?style=for-the-badge&logo=python&logoColor=white)
![Snakemake](https://img.shields.io/badge/Snakemake-333333?style=for-the-badge&logo=python&logoColor=white)
![Nextflow](https://img.shields.io/badge/Nextflow-23CC85?style=for-the-badge&logo=nextflow&logoColor=white)
![Docker](https://img.shields.io/badge/Docker-2496ED?style=for-the-badge&logo=docker&logoColor=white)
![Apptainer](https://img.shields.io/badge/Apptainer-2E6CE6?style=for-the-badge&logo=linuxcontainers&logoColor=white)

## Usage

Clone the repository:

```bash
git clone https://github.com/young-sudo/htskit-mini.git
cd htskit-mini
```

### Run with Snakemake

To run the workflow **via Snakemake**, install it first (Conda recommended):

```bash
conda install -c bioconda -c conda-forge snakemake
```

Alternatively, install via pip:

```bash
pip install snakemake
```

**Basic run:**

```bash
snakemake --config mode=mapper
```

**Available modes:**
The `mode` parameter accepts two options: `mapper` or `assembler`
- `mode=mapper` to run the read-mapping workflow
- `mode=assembler` to run the *de-novo* assembly workflow

Use parameter `--use-conda` to automatically create and use the Conda environment defined in the Snakefile or `env.yml` when running the workflow.

#### (Optional) Define parameters

The default parameters are:

```bash
snakemake \
  --use-conda \
  --config mode=mapper \
  input_reference=data/reference.fasta \
  input_reads=data/reads.fasta \
  output_mapping=results/output.txt
```

### Run with Nextflow

```bash

```

## Mapper

Implements a basic read mapping algorithm that aligns high-throughput sequencing reads to a given reference genome, generating alignment information for each read.

### Method

The NGS mapper uses an **FM-index-based approach** for efficient sequence alignment. The key steps and components are as follows:

#### FM-Index Construction
- **Downsampling:** By default, every 50 nucleotides.
- **BWT:** Burrows-Wheeler Transform of the sequence.
- **SA:** Suffix array of the sequence.
- **First column (`first`):** Counts of alphabet characters in the first column of the BWT.
- **T:** Original sequence text.

#### Implementation Details
- **SA Construction:** Karkkainen-Sanders algorithm.
- **BWT from SA:** `bwtFromSa` function.
- **FM-Index Class:** `FmIndex`.

#### Local Alignment
- **Smith-Waterman:** Local sequence alignment using dynamic programming (`smithWaterman` and `traceback_with_indices` functions).
- **Edit Cost:** Defined as +1 for a match and -1 for a mismatch.

This approach allows fast and memory-efficient mapping of sequencing reads to the reference genome using the FM-index for indexing and Smith-Waterman for accurate local alignment.

You can run the mapper directly:

```bash
python3 mapper.py reference.fasta reads.fasta output.txt
```

## Assembler

Performs de novo assembly of high-throughput sequencing reads into longer contiguous sequences (contigs), reconstructing the genome without requiring a reference.

### Method

The assembler uses a **hashmap-based k-mer approach** for de novo sequence assembly. The main ideas and steps are:

#### Core Concept
- Uses a **Python dictionary (hashmap)** to quickly check k-mer membership.
- Stores k-mer counts for all reads.
- Consumes more memory than a Bloom filter but preserves exact counts and avoids false positives.

#### Pre-processing
- Count all k-mers (k=17) across reads and store them in the dictionary.
- Filter out rare k-mers (counts ≤ 2).
- Select a set of **seeds**: the three most frequent k-mers ("solid k-mers").

#### Greedy Assembly Algorithm
- Extend each seed **in both directions** (forward and backward).
- Generate the next/previous k-mer based on the alphabet (A, C, T, G).
- Check if this k-mer exists in the counts.
- Extend the seed by adding the nucleotide from the **most frequent k-mer** among valid options.
- Reduce the k-mer count (e.g., divide by 2) to avoid cycles.
- After extension, check that the contig:
  - Has not been discovered previously.
  - Has length > 300 bp.
- Stop extending a seed when no further valid k-mers are found in the counts.

This greedy, count-based approach allows building contigs efficiently while avoiding repetitive loops and ensuring that only significant k-mers contribute to assembly.


Usage
```bash
./assembly input_reads.fasta output_contigs.fasta
```
