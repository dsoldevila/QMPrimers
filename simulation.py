#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  5 12:39:40 2019

@author: david
"""

import pandas as pd
import numpy as np
import math
import re
from common import*

class Simulation():
    def __init__(self, template, sample_size):
        """
        @param template
        @param sample_size: Number of genomes to take from the template.
        """
        self.template = template
        self.sample_size = sample_size
        
        try:
            self.primer_pairs = self.template["primerPair"].unique()
            
            f_reg = re.compile("^mismF[T|N\d+]$")
            self.mismF = list(filter(f_reg.match, template.columns.values))[0]
            
            r_reg = re.compile("^mismR[T|N\d+]$")
            self.mismR = list(filter(r_reg.match, template.columns.values))[0]
        except:
            raise ValueError("Template not valid")
        
        return

    def simulate(self, k, B, N):
        """
        @param k: Paramater used to calculte the concentration of each sample. Range(0,1).
        @param B: Amplification eficiency = B^(-missmatches)
        @returns For each primer, the mean and median of the Pearson Coef. between the original and the final concentration
        """
        self.k = k
        self.B = B
        index = list(range(N))
        index.append("ncombinations")
        self.raw_stats = pd.DataFrame(index=index)
        
        pplen = self.primer_pairs.shape[0]
        count = 0
        for pp in self.primer_pairs:
            template = self.template.loc[self.template["primerPair"] == pp]
            full_sample = template["fastaid"].unique()
            
            if(full_sample.shape[0]<N or full_sample.shape[0] < self.sample_size):
                logging.warning("Sample too small with primer pair "+str(pp)+". Skipping...")
                
            else:
                self.raw_stats[pp] = 0.0
                n_combinations = math.factorial(full_sample.shape[0])/(math.factorial(self.sample_size)
                *math.factorial(full_sample.shape[0]-self.sample_size))
                self.raw_stats.loc["ncombinations", pp] = n_combinations
                if(N > 0.5*n_combinations):
                    logging.warning("Only "+str(n_combinations)+" possible combinations with primer pair "+str(pp))
                    
                for i in range(N):
                    self.step = i #needed by update_stats
                    sample = self.get_random_sample(full_sample, self.sample_size, k)
                    amplified_sample = self.amplify(template, pp, sample, B)
                    self.update_stats(amplified_sample, pp)
            count+=1
            print("{0:.2f}".format((count/pplen)*100)+"%")
                

        self.cooked_stats = self.cook_stats(0.95)
        print("Simulation done!")
        return self.raw_stats, self.cooked_stats

    def get_random_sample(self, full_sample, sample_size, k):
        sample = pd.DataFrame(columns=["fastaid", "oprop", "fprop"])
        sample["fastaid"] = np.random.choice(full_sample, size=sample_size)
        #k = self.k
        Ck = 1/(1-(1-k)**sample_size)
        for i in range(sample.shape[0]):
            prop = Ck*k*(1-k)**i
            sample.loc[i, "oprop"] = prop
        return sample

    def amplify(self, template, pp, sample, B, mode=1):
        reduced_template = template.loc[template["fastaid"].isin(sample["fastaid"])]
        
        a = 0.0

        pp_matches = reduced_template.loc[reduced_template["primerPair"] == pp]
        #TODO Could be the case that there's been a multiple alignment, use the best case
        
        for i in sample.index.values:
            row = pp_matches.loc[pp_matches["fastaid"]==sample.loc[i, "fastaid"]]
            #TODO use N end missmatches if mode=2
            m = row[self.mismF].values[0] + row[self.mismR].values[0] #TODO make it more clean(?)
            amp_eff = 1/(B**m)
            a += amp_eff*sample.loc[i, "oprop"]
            sample.loc[i, "fprop"] = amp_eff
            
        #compute final proportion, note that the previous fprop is amplification efficiency
        sample["fprop"] = sample["fprop"]*(sample["oprop"]/a)
            
        return sample

    def update_stats(self, sample, primer_pair): 
        tmp = sample[["oprop", "fprop"]]
        tmp = tmp.astype('float64') #TODO patch, create a better solution if slow
        tmp = tmp.corr()
        tmp = tmp.loc["oprop", "fprop"]
        self.raw_stats.loc[self.step, primer_pair] = tmp
        return
    
    def cook_stats(self, ci):
        cooked_stats = pd.DataFrame(index=self.raw_stats.columns, columns=["min", "max", "mean", "median", "ncombinations", "CI"])
        raw_stats = self.raw_stats.loc[self.raw_stats.index[:-1]]
        
        raw_stats = np.sort(raw_stats, axis=0)
    
        low_interval = int(((1-ci)/2)*raw_stats.shape[0])
        high_interval = int((1-((1-ci)/2))*raw_stats.shape[0])
    
        raw_stats = raw_stats[low_interval:high_interval, :]
        
        cooked_stats["min"] = raw_stats[0,:]
        cooked_stats["max"] = raw_stats[-1, :]
        cooked_stats["CI"] = ci
        cooked_stats["ncombinations"] = self.raw_stats.loc["ncombinations"]
        cooked_stats["mean"] = raw_stats.mean(axis=0)
        cooked_stats["median"] = np.median(raw_stats, axis=0)
    
        return cooked_stats
    
    @staticmethod
    def store_raw_data(output_file, raw_stats, cooked_stats, sample_size, k, B, N):
        raw_stats.to_csv(output_file+".csv", index_label="Step")
        
        parameters_used = "SampleSize = "+str(sample_size)+"\t\tk = "+str(k)+"\t\tBeta = "+str(B)+"\t\tN = "+str(N)+"\n"
        
        with open(output_file+".txt",'w') as outfile:
            outfile.write(parameters_used)
            cooked_stats.to_string(outfile)
            
        print("Saved!")
        return
    
    def plot_results(self, cooked_stats):
        #boxplot = cooked_stats.boxplot()
        return
    
    def store_plot(self):
        return



