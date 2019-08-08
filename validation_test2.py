#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 12:58:15 2019

@author: david
"""

import QMPrimers as qm
import os
import interface as ie
import pandas as pd
import numpy as np


rootd = os.getcwd()
ind = os.path.join(rootd, "test_input")
outd = os.path.join(rootd, "test_output")
vald = os.path.join(rootd, "test_validated")


def check_template(template, val_template):
    template = ie.load_template_only(template)
    val_template = ie.load_template_only(val_template)
    
    result = pd.merge(template, val_template, on=list(template.columns.values), how='left', indicator='Correct')
    del result["amplicon"] #amplicon too long
    result['Correct'] = np.where(result.Correct == 'both', True, False)
    
    wrong_results = result[result['Correct']==False]
    
    return wrong_results.empty, wrong_results
    
        
    

def matching_test():
    argv = ["--match", "-gf", "GENOME", "-pf", "PRIMER", "-o", "OUTPUT", "-fm", 5, "-rm", 5]
    match_input = [["sbog_test.fasta", "PP.csv", "valtest_test"]]
    match_output =[]
    
    for m in match_input:
        fasta = os.path.join(ind, m[0])
        primer = os.path.join(ind, m[1])
        output = os.path.join(outd, m[2])
        argv[2] = fasta
        argv[4] = primer
        argv[6] = output
        qm.main(argv)
        
        is_correct, wrong_results = check_template(output+"_corrupt.csv", output+"_positive.csv")
        match_output.append([m[0], m[1], is_correct, wrong_results])
    
    for m in match_output:
        print(m[:-1])
        print(m[-1].to_string())
    

def check_sim():
    #load sim file
    #compare
    pass

def sim_test():
    argv = ["--sim", "-i", "TEMPLATE", "-s", 10, "-b", 4, "-k", 0.5, "-n", 10 , "-o", "OUTPUT", "-ci", 0.95]
    match_input = [["insects_positive.csv", "valtest_test_sim"]]
    
    for m in match_input:
        template = os.path.join(outd, m[0])
        output = os.path.join(outd, m[1])
        argv[2] = template
        argv[12] = output
        qm.main(argv)
    
    

if(__name__=="__main__"):
    matching_test()