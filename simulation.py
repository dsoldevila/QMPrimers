#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 12:39:40 2019

@author: david
"""

import pandas as pd
from common import *



def simulate(template, sample_size, K, B, N=100):
    """
    @param template
    @param sample_size: Number of genomes to take from the template.
    @param k: Paramater used to calculte the concentration of each sample. Range(0,1).
    @param B: Amplification eficiency = B^(-missmatches)
    @returns For each primer, the mean and median of the Pearson Coef. between the original and the final concentration
    """
    for i in N:
        sample = get_random_sample(template, sample_size, K)
        amplified_sample = amplify(template, sample, B)
        get_stats(sample, amplified_sample)
    return

def get_random_sample(template, sample_size, K):
    #Get random genomes
    #Calculate concentration with K
    #return a table containing for each genome, its concentration
    return

def amplify(template, sample, B):
    #calculate amplification efficiency
    #calculate final concentration
    return

def get_stats():
    #calculate Perason coef for each primer
    #append the result to the global statitstics
    return



