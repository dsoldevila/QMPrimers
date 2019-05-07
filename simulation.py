#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 12:39:40 2019

@author: david
"""

import pandas as pd
from common import *
import numpy as np

class Simulation():
    def __init__(self, template, sample_size):
        """
        @param template
        @param sample_size: Number of genomes to take from the template.
        """
        self.template = template
        self.sample_size = sample_size
        
        self.full_sample = self.template["fastaid"].unique()
        self.primer_pairs = self.template["primerPair"].unique()
        
        return

    def simulate(self, k, B, N=100):
        """
        @param k: Paramater used to calculte the concentration of each sample. Range(0,1).
        @param B: Amplification eficiency = B^(-missmatches)
        @returns For each primer, the mean and median of the Pearson Coef. between the original and the final concentration
        """
        self.k = k
        self.B = B
        self.raw_stats = pd.DataFrame(index=self.primer_pairs, columns=range(N))
        
        for i in range(N):
            self.step = i #needed by update_stats
            sample = self.get_random_sample(self.full_sample, self.sample_size, k)
            for pp in self.primer_pairs:
                amplified_sample = self.amplify(self.template, pp, sample, B)
                self.update_stats(amplified_sample, pp)

        self.cooked_stats = self.cook_stats()
        print(self.raw_stats)
        print(self.cooked_stats)
        return

    def get_random_sample(self, full_sample, sample_size, k):
        sample = pd.DataFrame(columns=["fastaid", "oprop"])
        sample["fastaid"] = np.random.choice(full_sample, size=sample_size)
        #k = self.k
        Ck = 1/(1-(1-k)**sample_size)
        for i in range(sample.shape[0]):
            prop = Ck*k*(1-k)**i
            sample.loc[i, "oprop"] = prop
        return sample

    def amplify(self, template, pp, sample, B):
        reduced_template = template.loc[template["fastaid"].isin(sample["fastaid"])]
        
        a = 0.0
        sample["fprop"] = 0.0
        pp_matches = reduced_template.loc[reduced_template["primerPair"] == pp]
        #TODO Could be the case that there's been a multiple alignment, use the best case
        
        if(not pp_matches.empty):
            for i in sample.index.values:
                row = pp_matches.loc[pp_matches["fastaid"]==sample.loc[i, "fastaid"]]
                try:
                    m = row["mismFT"].values[0] + row["mismRT"].values[0] #TODO make it more clean(?)
                    amp_eff = 1/(B**m)
                except:
                    amp_eff = 0
                    print("No matches")
                a += amp_eff*sample.loc[i, "oprop"]
                sample.loc[i, "fprop"] = amp_eff
            
            #compute final proportion, note that the previous fprop is amplification efficiency
            sample["fprop"] = sample["fprop"]*(sample["oprop"]/a)
            
        return sample

    def update_stats(self, sample, primer_pair): 
        tmp = sample[["oprop", "fprop"]]
        tmp = tmp.astype('float64') #TODO patch, create a better solution of slow
        tmp = tmp.corr()
        tmp = tmp.loc["oprop", "fprop"]
        self.raw_stats.loc[primer_pair, self.step] = tmp
        return
    
    def cook_stats(self):
        cooked_stats = pd.DataFrame(index=self.primer_pairs, columns=["mean"])
        cooked_stats["mean"] = self.raw_stats.sum(axis=1)/self.raw_stats.shape[1]
        return cooked_stats



