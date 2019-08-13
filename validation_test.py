#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 12:58:15 2019

@author: david
@brief Passing this test should be enough to publish a new version of QMPrimers. Currently there are no validated results
to compare to.
"""

import QMPrimers as qm
import os
import interface as ie
import pandas as pd
import numpy as np
import csv
import sys


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

    
def get_tests_input(infile):
    tests = []
    with open(infile, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        for c in csvreader:
            tests.append(c)
    return tests
    
def sim_test(test):
    pass

def matching_test(test):
    val_positive = os.path.basename(test[6])
    val_positive = os.path.join(vald, val_positive)
    correct, wrong_results = check_template(test[6]+"_positive.csv", val_positive+"_positive.csv")
    if not correct:
        with open(os.path.join(outd, "test.log"),'a') as outfile:
            outfile.write(str(test)+"\n")
            outfile.write(wrong_results.to_string()+"\n")

    return correct

def perform_test(test):
    if True:
        matching_test(test)
    else:
        sim_test(test)
    return
    

if(__name__=="__main__"):
    if len(sys.argv) > 1:
        infile = sys.argv[1]
    else:
        infile = "test_input/testflow.csv"
        
    tests = get_tests_input(infile)
    
    for i in range(len(tests)):
        print("TEST ", i)
        print(tests[i])
        qm.main(tests[i])
        """
        if(perform_test(tests[i])):
            print("TEST ", i," SUCCESS")
        else:
            print("TEST ", i," FAIL")
            print("TEST FAILED :(")
            break
        """