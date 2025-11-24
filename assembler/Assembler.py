#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from copy import copy

class Assembler():

    def __init__(self, k, d, alpha='ACTG'):
        
        self.kmer_counts = {}
        self.k = k # length of k-mers
        self.d = d # min k-mer count threshold
        self.alpha = alpha # alphabet



    def count_kmers_from_read(self, read: str):
        """"""
        for i in range(len(read) - self.k + 1):
            kmer = read[i:i+self.k]
            if kmer in self.kmer_counts:
                self.kmer_counts[kmer] += 1
            else:
                self.kmer_counts[kmer] = 1

    def filter_dict_by_threshold(self):
        # keep if kmer count >= thresh
        self.kmer_counts = {key: value for key, value in self.kmer_counts.items() if value >= self.d}

    def hist2distr(self):
        vals = list(self.kmer_counts.values())
        binwidth = 1
        plt.hist(vals,
                 color="lightsteelblue",
                 alpha=0.7,
                 bins=range(min(vals) - 1,
                            max(vals) + binwidth + 1, binwidth))
        plt.xlabel("K-mer Length")
        plt.ylabel("Frequency")
        plt.show()
    
    
    def extend_seed(self, seed, greedy=False):

        if not greedy:
            tmp_kmer_counts = copy(self.kmer_counts)
        else:
            tmp_kmer_counts = self.kmer_counts
            
        og_seed = seed # original seed

        # forward pass
        keep_extending = True
        kmer = og_seed
        while keep_extending:

            counts = {}
            for c in self.alpha:
                option = ''.join((kmer[1:], c))
                if option in tmp_kmer_counts:
                    counts[option] = tmp_kmer_counts[option]
                
            if counts:
                kmer = max(counts, key=counts.get)
                seed = ''.join([seed, kmer[-1]])
                if not greedy:
                    tmp_kmer_counts[kmer] /= 2
            else:
                keep_extending = False
            
            
        # backward pass
        keep_extending = True
        kmer = og_seed

        while keep_extending:
            counts = {}
            for c in self.alpha:
                option = ''.join((c, kmer[:-1]))
                if option in tmp_kmer_counts:
                    counts[option] = tmp_kmer_counts[option]
                
            if counts:
                kmer = max(counts, key=counts.get)
                seed = ''.join([kmer[0], seed])
                if not greedy:
                    tmp_kmer_counts[kmer] /= 2
            else:
                keep_extending = False
            

        return seed