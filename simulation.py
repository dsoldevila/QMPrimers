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
    def __init__(self):
        return
    
    def get_mism(self, template):
        self.template = template
        try:

            self.mismF = get_missmatch_column_name(template.columns.values, primer="f")
            self.mismR = get_missmatch_column_name(template.columns.values, primer="r")
            self.primer_pairs = self.template["primerPair"].unique()
            
            #np.random.seed(0) #TODO For debugging purposes
        except:
            raise Exception("Template not valid")
        
    def check_data(self, sample_size, k, B, N, ci):
        assert(sample_size>0), "Sample Size must be greater than 0"
        assert(k>0 and k<1), "K must be in the range (0,1)"
        assert(B>0), "Beta must be greater than 0"
        assert(N>0), "N must be greater than 0"
        assert(ci>0 and ci<1), "Confidence interval must be in the range (0,1)"
        
    def simulate(self, template, sample_size, k, B, N, ci):
        """
        @param k: Paramater used to calculte the concentration of each sample. Range(0,1).
        @param B: Amplification eficiency = B^(-missmatches)
        @returns For each primer, the mean and median of the Pearson Coef. between the original and the final concentration
        """
        #check input data:
        try:
            self.get_mism(template)
            self.check_data(sample_size, k, B, N, ci)
        except(Exception) as e:
            logging.error(e)
            return pd.DataFrame(), pd.DataFrame()

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
            
            try:
                n_combinations = math.factorial(full_sample.shape[0])/(math.factorial(sample_size)
                *math.factorial(full_sample.shape[0]-sample_size))
            except: #exception if full_sample < sample_size
                n_combinations = 0
                pass
                
            if(n_combinations<sample_size):
                logging.error("Sample size too small with primer pair "+str(pp)+". Skipping...")
                
            else:
                if(n_combinations<N):
                    logging.warning("With primer pair "+str(pp)+" only "+str(n_combinations)+" possible combinations")
                self.raw_stats[pp] = 0.0
                self.raw_stats.loc["ncombinations", pp] = n_combinations
                    
                for i in range(N):
                    self.step = i #needed by update_stats
                    sample = self.get_random_sample(full_sample, sample_size, k)
                    amplified_sample = self.amplify(template, pp, sample, B)
                    self.update_stats(amplified_sample, pp)
            count+=1
            print("{0:.2f}".format((count/pplen)*100)+"%")
                
        self.cooked_stats = self.cook_stats(ci)
        print("Simulation done!")
        self.raw_stats = self.raw_stats.drop("ncombinations")
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
        if(ci>1.0): ci=1.0
        cooked_stats = pd.DataFrame(index=self.raw_stats.columns, columns=["min", "max", "mean", "median", "ncombinations", "CI", "min_ci", "max_ci"])
        raw_stats = self.raw_stats.loc[self.raw_stats.index[:-1]]
        #TODO better way to get data within CI?
        raw_stats = np.sort(raw_stats, axis=0)
        cooked_stats["min"] = raw_stats[0,:]
        cooked_stats["max"] = raw_stats[-1, :]
        
        low_interval = int(((1-ci)/2)*raw_stats.shape[0])
        high_interval = int((1-((1-ci)/2))*raw_stats.shape[0])+1
        raw_stats = raw_stats[low_interval:high_interval, :]

        cooked_stats["min_ci"] = raw_stats[0,:]
        cooked_stats["max_ci"] = raw_stats[-1, :]
        cooked_stats["CI"] = ci
        cooked_stats["ncombinations"] = self.raw_stats.loc["ncombinations"]
        cooked_stats["mean"] = raw_stats.mean(axis=0)
        cooked_stats["median"] = np.median(raw_stats, axis=0)
    
        return cooked_stats
    
    @staticmethod
    def store_data(output_file, raw_stats, cooked_stats, infile_name, sample_size, k, B, N):
        
        try:
            parameters_used = "Template = "+infile_name+"     S = "+str(sample_size)+"     k = "+str(k)+"     Beta = "+str(B)+"     N = "+str(N)+"\n"
        except:
            logging.error("Could not save files, bad parameters")
            return
        
        try:
            with open(output_file+".txt",'w') as outfile:
                outfile.write(parameters_used)
                outfile.write(str(datetime.datetime.now())+"\n")
                cooked_stats.to_string(outfile)
                print("Statistics saved")
        except:
            logging.error("Could not save statistics")
                
        try:
            with open(output_file+".csv",'w') as outfile:  
                outfile.write(parameters_used)
                outfile.write(str(datetime.datetime.now())+"\n")
                raw_stats.to_csv(outfile, mode='a', index_label="Step")
            print("Raw data saved!")
        except:
            logging.error("Could not save raw data")
            
        return
    
    def plot_results(self, cooked_stats):
        #boxplot = cooked_stats.boxplot()
        return
    
    def store_plot(self):
        return



