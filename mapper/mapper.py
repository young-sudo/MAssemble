#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mapper
---

2023 - Younginn Park
"""

import argparse
import timeit

from Bio import SeqIO
from fmindex import FmIndexApprox, FmIndexWithT
from fmhelper import editCost, traceback_with_indices

def argument_parser():
    """ Parsing command line arguments """
    parser = argparse.ArgumentParser(
        description="A script for mapping reads to a reference.",
        epilog="Enjoy :)"
    )

    # Positional arguments
    parser.add_argument("reference", help="Path to the reference.fasta file")
    parser.add_argument("reads", help="Path to the reads.fasta file")
    parser.add_argument("output", help="Path to the output.txt file")

    args = parser.parse_args()

    reference_file = args.reference
    reads_file = args.reads
    output_file = args.output

    return reference_file, reads_file, output_file


def naive_mapper(reference_file, reads_file, output_file):
    """ Mapper that naively assumes that mapped length = read length """
    # It's faster, maps 100 reads <1s even for references 20M
    # But it can be inaccurate (max. error ~17 nt)

    reference_sequence = str(SeqIO.read(reference_file, "fasta").seq)
    reads = SeqIO.parse(reads_file, "fasta")

    # FM index construction
    print("Started FM index construction")
    start_time = timeit.default_timer()
    fm_index = FmIndexApprox(reference_sequence)
    end_time = timeit.default_timer()

    execution_time = end_time - start_time
    print(f"Executed FM construction in: {execution_time:.3f} seconds")

    # Mapping
    print("Started mapping")
    start_time = timeit.default_timer()
    
    mapping = {}

    not_mapped = 0 # how many reads were not mapped
    num_reads = 0

    for read in reads:
        num_reads += 1

        # klen - kmer length, kint - kmer interval
        # found to work best for ~1kbp reads with ~10% error rate (23, 17)
        # shorter kmer - more probable to find exact match, but longer execution
        # longer kmer - faster, less exact occurrences, but less likely to find match
        # klen and kmer reflect prior belief of error occurrences and distribution
        klen = 23 # this could be a parameter
        kint = 17
        rlen = len(read.seq)

        kmers, inds = fm_index.generate_spaced_kmers(read, k=klen, ki=kint)

        kmer_occs = {}
        for it in range(len(kmers)):
            kmer = kmers[it]
            i = inds[it]
            kmer_occs[i] = fm_index.occurrences(kmer.seq)

        found = sum(map(len, kmer_occs.values()))
        if not found:
            # If no exact match for any of the k-mers is found
            not_mapped += 1
            continue

        # "mapping"
        proposals = {}
        for kmer, match_list in kmer_occs.items():
            for match_position in match_list:

                pstart = match_position - kmer
                pstop = match_position + (rlen - kmer)

                if pstart < 0 or pstop > fm_index.slen:
                    continue

                proposal = (pstart, pstop)
                if (pstart, pstop) in proposals.keys():
                    proposals[proposal] += 1
                else:
                    proposals[proposal] = 1
                    
        mapping[read.id] = max(proposals, key=proposals.get)

    end_time = timeit.default_timer()
    execution_time = end_time - start_time
    print(f"Executed mapping in: {execution_time:.5f} seconds")

    # Save results
    with open(output_file, "w") as f:
        for read_id, positions in mapping.items():
            f.write(f"{read_id}\t{positions[0]}\t{positions[1]}\n")

    mapped = num_reads - not_mapped
    print(f"Successfully mapped {mapped} out of {num_reads} reads: {100*(mapped)/num_reads}%")

def accurate_mapper(reference_file, reads_file, output_file):
    """A much accurate mapper """

    reference_sequence = str(SeqIO.read(reference_file, "fasta").seq)
    reads = SeqIO.parse(reads_file, "fasta")

    # FM index construction
    print("Started FM index construction")
    start_time = timeit.default_timer()
    fm_index = FmIndexWithT(reference_sequence)
    end_time = timeit.default_timer()

    execution_time = end_time - start_time
    print(f"Executed FM construction in: {execution_time:.3f} seconds")

    # Mapping
    start_time = timeit.default_timer()
    s = editCost
    mapping = {}

    not_mapped = 0 # how many reads were not mapped
    num_reads = 0

    for read in reads:
        num_reads += 1

        # klen - kmer length, kint - kmer interval
        # found to work best for ~1kbp reads with ~10% error rate (23, 17)
        # shorter kmer - more probable to find exact match, but longer execution
        # longer kmer - faster, less exact occurrences, but less likely to find match
        # klen and kmer reflect prior belief of error occurrences and distribution
        klen = 23
        kint = 17
        rlen = len(read.seq)

        kmers, inds = fm_index.generate_spaced_kmers(read, k=klen, ki=kint)

        kmer_occs = {}
        for it in range(len(kmers)):
            kmer = kmers[it]
            i = inds[it]
            kmer_occs[i] = fm_index.occurrences(kmer.seq)

        found = sum(map(len, kmer_occs.values()))
        if not found:
            # If no exact match for any of the k-mers is found
            not_mapped += 1
            continue

        proposals = {}
        proposals_starts = {} # starting positions of kmers involved in a proposal
        for kmer, match_list in kmer_occs.items():
            for match_position in match_list:

                pstart = match_position - kmer
                pstop = match_position + (rlen - kmer)

                if pstart < 0 or pstop > fm_index.slen:
                    continue

                proposal = (pstart, pstop)
                if (pstart, pstop) in proposals.keys():
                    proposals[proposal] += 1
                    proposals_starts[proposal].append(match_position)
                else:
                    proposals[proposal] = 1
                    proposals_starts[proposal] = [match_position]
        
        # Edge refinement
        # best proposal - the one with the most confirmation
        kmer_starts = proposals_starts[max(proposals, key=proposals.get)]
        naive_proposal = max(proposals, key=proposals.get)

        # indices of T and read that will be aligned
        safety = 5
        tl_start = max(naive_proposal[0] - safety, 0)
        tl_end = min(kmer_starts)
        tr_start = max(kmer_starts) + klen
        tr_end = min(naive_proposal[1] + safety, fm_index.slen)

        pl_start = 0
        pl_end = tl_end - tl_start + safety
        pr_start = tr_start - tl_start - safety
        pr_end = rlen

        tl = fm_index.t[tl_start : tl_end]
        pl = read.seq[pl_start : pl_end]
        tr = fm_index.t[tr_start : tr_end]
        pr = read.seq[pr_start : pr_end]

        tl_inds, _ = traceback_with_indices(tl, pl, s)
        tr_inds, _ = traceback_with_indices(tr, pr, s)

        new_start = tl_inds[0] + tl_start
        new_end = min(tr_inds[1] + tr_start, fm_index.slen)

        new_proposal = (new_start, new_end)

        mapping[read.id] = new_proposal

    end_time = timeit.default_timer()
    execution_time = end_time - start_time
    print(f"Executed mapping in: {execution_time:.5f} seconds")

    # Save results
    with open(output_file, "w") as f:
        for read_id, positions in mapping.items():
            f.write(f"{read_id}\t{positions[0]}\t{positions[1]}\n")

    mapped = num_reads - not_mapped
    print(f"Successfully mapped {mapped} out of {num_reads} reads: {100*(mapped)/num_reads}%")


def main():
    reference_file, reads_file, output_file = argument_parser()

    # naive_mapper(reference_file, reads_file, output_file) # faster, less accurate

    accurate_mapper(reference_file, reads_file, output_file) # slower, more accurate


if __name__ == "__main__":
    main()
