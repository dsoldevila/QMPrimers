#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 12:39:40 2019

@author: david
"""

import pandas as pd
from common import *
import numpy as np



def simulate(template, sample_size, K, B, N=100):
    """
    @param template
    @param sample_size: Number of genomes to take from the template.
    @param k: Paramater used to calculte the concentration of each sample. Range(0,1).
    @param B: Amplification eficiency = B^(-missmatches)
    @returns For each primer, the mean and median of the Pearson Coef. between the original and the final concentration
    """
    full_sample = template["fastaid"].unique()
    primer_pairs = template["primerPair"].unique()
    
    for i in N:
        sample = get_random_sample(template, sample_size, K)
        amplified_sample = amplify(template, sample, B)
        update_stats(sample, amplified_sample)
    return

def get_random_sample(full_sample, sample_size, k):
    sample = pd.DataFrame(columns=["fastaid", "oprop"])
    sample["fastaid"] = np.random.choice(full_sample, size=sample_size)
    
    Ck = 1/(1-(1-k)**sample_size)
    for i in range(sample.shape[0]):
        prop = Ck*k*(1-k)**i
        sample.loc[i, "oprop"] = prop
    print("Sum O =", str(sample["oprop"].sum()))
    return sample

def amplify(template, primer_pairs, sample, B):
    reduced_template = template.loc[template["fastaid"].isin(sample["fastaid"])]
    
    for pp in primer_pairs:
        a = 0.0
        sample["fprop"] = 0.0
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
        #update_stats(pp, sample)
    return

def update_stats(sample, primer_pair):
    sample = sample.astype('float64') 
    tmp = sample[["oprop", "fprop"]]
    tmp = tmp.astype('float64') #TODO patch, create a better solution of slow
    tmp = tmp.corr()
    #append the result to the global statitstics
    return sample.corr()



