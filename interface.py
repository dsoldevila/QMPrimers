#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 13:31:13 2019

@author: david
"""

import load_data as ld
import matching as m
from common import *

def compute(parameters):
    """
    @Brief calls the matching algorithms for both GUI and command line modes
    """
    template = None
    gen_record = load_gen_record(parameters.loc["gen", "value"], parameters.loc["check_integrity", "value"], 
                                 parameters.loc["check_uppercase", "value"], parameters.loc["hanging primers", "value"])
    primer_pairs = load_primer_pairs(parameters.loc["primer_pairs", "value"])
    if(gen_record!=None and primer_pairs!=None):
        template = m.compute_gen_matching(int(parameters.loc["forward missmatches", "value"]), int(parameters.loc["reverse missmatches", "value"]), 
                                          primer_pairs, gen_record, parameters.loc["Nend miss.", "value"], hanging_primers=parameters.loc["hanging primers", "value"])
    return template, gen_record, primer_pairs

def save_template_primer_missmatches(output_file, template, header=None):
    m.store_matching_results(output_file, template, header)
    return

def load_gen_record(gen_file, check_integrity, check_uppercase, hanging_primers):
    if(check_integrity): print("INTEGRITY CHECKED")
    gen_record = None
    try:
        gen_record = ld.load_bio_files(gen_file, writable=check_integrity or check_uppercase, 
                                       check_uppercase=check_uppercase, file_format=None)
    except:
        print("Error at loading gen record")
    if(check_integrity): gen_record = ld.remove_bad_gens(gen_record)
    return gen_record

def load_primer_pairs(primer_pairs_file):
    primer_pairs = None
    try:
        primer_pairs = ld.load_csv_file(primer_pairs_file)
    except:
        print("Error at loading primer pairs file")
        
    return primer_pairs

def load_template(parameters):
    template = ld.load_template(parameters.loc["csv_template", "value"])
    
    gen_record = load_gen_record(parameters.loc["gen", "value"], parameters.loc["check_integrity", "value"], 
                                 parameters.loc["check_uppercase", "value"], parameters.loc["hanging primers", "value"])
    primer_pairs = load_primer_pairs(parameters.loc["primer_pairs", "value"])
    if(gen_record!=None and primer_pairs!=None):
        template = ld.restore_template(template, gen_record, primer_pairs)
        
    return template
        

def simulate():
    return