#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 12:39:40 2019

@author: david
"""

import pandas as pd
from common import *
import numpy as np
import math

class Simulation():
    def __init__(self, template, sample_size):
        """
        @param template
        @param sample_size: Number of genomes to take from the template.
        """
        self.template = template
        self.sample_size = sample_size
        
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
        columns = list(range(N))
        columns.append("ncombinations")
        self.raw_stats = pd.DataFrame(columns=range(N))
        
        for pp in self.primer_pairs:
            template = self.template.loc[self.template["primerPair"] == pp]
            full_sample = template["fastaid"].unique()
            
            if(full_sample.shape[0]<N):
                print("Warning: sample too small with primer pair ", str(pp), ". Skipping...")
                
            else:
                self.raw_stats.loc[pp] = 0.0
                n_combinations = math.factorial(full_sample.shape[0])/(math.factorial(self.sample_size)
                *math.factorial(full_sample.shape[0]-self.sample_size))
                self.raw_stats.loc[pp, "ncombinations"] = n_combinations
                if(N > 0.5*n_combinations):
                    print("Warning: Only ", str(n_combinations), " possible combinations with primer pair ", str(pp))
                    
                for i in range(N):
                    self.step = i #needed by update_stats
                    sample = self.get_random_sample(full_sample, self.sample_size, k)
                    amplified_sample = self.amplify(template, pp, sample, B)
                    self.update_stats(amplified_sample, pp)

        self.cooked_stats = self.cook_stats(0.95)
        print(self.raw_stats)
        print(self.cooked_stats)
        return

    def get_random_sample(self, full_sample, sample_size, k):
        sample = pd.DataFrame(columns=["fastaid", "oprop", "fprop"])
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

        pp_matches = reduced_template.loc[reduced_template["primerPair"] == pp]
        #TODO Could be the case that there's been a multiple alignment, use the best case
        
        for i in sample.index.values:
            row = pp_matches.loc[pp_matches["fastaid"]==sample.loc[i, "fastaid"]]
            m = row["mismFT"].values[0] + row["mismRT"].values[0] #TODO make it more clean(?)
            amp_eff = 1/(B**m)
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
    
    def cook_stats(self, ci):
        cooked_stats = pd.DataFrame(index=self.raw_stats.index, columns=["mean", "median", "ncombinations", "CI"])
        raw_stats = self.raw_stats[self.raw_stats.columns[:-1]]
        
        raw_stats = np.sort(raw_stats)
    
        low_interval = int(((1-ci)/2)*raw_stats.shape[1])
        high_interval = int((1-((1-ci)/2))*raw_stats.shape[1])
    
        raw_stats = raw_stats[:,low_interval:high_interval]
        
        cooked_stats["CI"] = ci
        cooked_stats["ncombinations"] = self.raw_stats["ncombinations"]
        cooked_stats["mean"] = raw_stats.mean(axis=1)
        cooked_stats["median"] = np.median(raw_stats, axis=1)
    
        return cooked_stats



