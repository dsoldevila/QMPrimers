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
        template, discarded, raw_stats, cooked_stats = m.compute_gen_matching(int(parameters.loc["forward missmatches", "value"]), int(parameters.loc["reverse missmatches", "value"]), 
                                          primer_pairs, gen_record, parameters.loc["Nend miss.", "value"], hanging_primers=parameters.loc["hanging primers", "value"])
    return template, discarded, gen_record, primer_pairs, raw_stats, cooked_stats

def save_matching_info(output_file, template, discarded, raw_stats, cooked_stats, header=None):
    m.store_matching_results(output_file+"_positive.csv", template, header)
    
    if(not discarded.empty):
        m.store_discarded(output_file+"_negative.csv", discarded)
    
    m.store_stats(output_file+"_stats.txt", raw_stats, cooked_stats)
    
    return

def load_gen_record(gen_file, check_integrity, check_uppercase, hanging_primers):
    if(check_integrity): print("INTEGRITY CHECKED")
    gen_record = None
    try:
        gen_record = ld.load_bio_files(gen_file, writable=check_integrity or check_uppercase, 
                                       check_uppercase=check_uppercase, file_format=None)
        print("Genome record file loaded!")
    except:
        print("Error at loading gen record")
    if(check_integrity): gen_record = ld.remove_bad_gens(gen_record)
    return gen_record

def load_primer_pairs(primer_pairs_file):
    primer_pairs = None
    try:
        primer_pairs = ld.load_csv_file(primer_pairs_file)
        print("Primer pairs file loaded!")
    except:
        print("Error at loading primer pairs file")
        
    return primer_pairs

def load_template(parameters):
    gen_record = load_gen_record(parameters.loc["gen", "value"], parameters.loc["check_integrity", "value"], 
                                 parameters.loc["check_uppercase", "value"], parameters.loc["hanging primers", "value"])
    primer_pairs = load_primer_pairs(parameters.loc["primer_pairs", "value"])
    try:
        template = ld.load_template(parameters.loc["csv_template", "value"])
        max_misses = int(parameters.loc["forward missmatches", "value"]) + int(parameters.loc["reverse missmatches", "value"])
        template, discarded, raw_stats, cooked_stats = ld.restore_template(template, gen_record, primer_pairs, max_misses)
        print("Template file restored!")
    except:
        print("Error at restoring template")
    return template, discarded, gen_record, primer_pairs, raw_stats, cooked_stats
        

def simulate():
    return