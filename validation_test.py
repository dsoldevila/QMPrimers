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


def check_match_result(template_file, val_template_file):
    template = ie.load_template_only(template_file+"_positive.csv")
    negative = pd.read_csv(template_file+"_negative.csv", index_col=0, skiprows=3)
    val_template = ie.load_template_only(val_template_file+"_positive.csv")
    val_negative = ie.load_template_only(val_template_file+"_negative.csv")
    template_check = pd.merge(template, val_template, on=list(template.columns.values), how='left', indicator='Correct')
    del template_check["amplicon"] #amplicon too long
    template_check['Correct'] = np.where(template_check.Correct == 'both', True, False)
    
    negative_check = pd.merge(negative, val_negative, on=list(negative.columns.values), how='left', indicator='Correct')
    negative_check['Correct'] = np.where(negative_check.Correct == 'both', True, False)
    
    bad_template_entries = template_check[template_check['Correct']==False]
    bad_negative_entries = negative_check[negative_check['Correct']==False]
    return bad_template_entries.empty, bad_template_entries, bad_negative_entries.empty, bad_negative_entries

    
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
    tcorrect, template, ncorrect, negative = check_match_result(test[6], val_positive)
    if not tcorrect or not ncorrect:
        with open(os.path.join(outd, "test.log"),'a') as outfile:
            outfile.write(str(test)+"\n")
            if not tcorrect:
                outfile.write(template.to_string()+"\n")
            if not ncorrect:
                outfile.write(negative.to_string()+"\n")
    
    

    return tcorrect and ncorrect

def perform_test(test):
    if True:
        return matching_test(test)
    else:
        return sim_test(test)
    

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

        if(perform_test(tests[i])):
            print("TEST ", i," SUCCESS")
        else:
            print("TEST ", i," FAIL")
            print("TEST FAILED :(")
            break
        