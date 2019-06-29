#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 13:31:13 2019

@author: david
@brief: This document was made to make the front end a little cleaner and to separate the front end from the backend, to not have to modify both
command line mode and gui mode when making changes on the back end. Now that the project is almost complete, maybe this file has no sense anymore. 
Especially if we take into account that the simulation backend bypasses this file.

EDIT: Now handles the Queue to communicate the TKINTER GUI with the backend
"""

import load_data as ld
import matching as m
from common import *
import simulation as s


"""
Matching Stuff
"""
def compute(parameters, thread_q=None):
    """
    @Brief calls the matching algorithms for both GUI and command line modes
    """
    template = None
    gen_record = load_gen_record(parameters.loc["gen", "value"], parameters.loc["check_integrity", "value"], 
                                 parameters.loc["check_uppercase", "value"], parameters.loc["hanging primers", "value"])
    primer_pairs = load_primer_pairs(parameters.loc["primer_pairs", "value"])
    if(gen_record!=None and primer_pairs!=None):
        template, discarded, raw_stats, cooked_stats = m.compute_gen_matching(int(parameters.loc["forward missmatches", "value"]), int(parameters.loc["reverse missmatches", "value"]), 
                                          primer_pairs, gen_record, parameters.loc["output_file", "value"], hanging_primers=parameters.loc["hanging primers", "value"])
        if(thread_q!=None): #TODO better queue this data in load_gen and load_primer_pairs respectively
            thread_q.put(template)
            thread_q.put(discarded)
            thread_q.put(raw_stats)
            thread_q.put(cooked_stats)
            thread_q.put(gen_record)
            thread_q.put(primer_pairs)
    return template, discarded, gen_record, primer_pairs, raw_stats, cooked_stats

def save_matching_info(output_file, template, header, discarded, raw_stats, cooked_stats):
    
    #try:
    m.store_matching_results(output_file+"_positive.csv", template, header)
    print("Template saved")
    
    if(not discarded.empty):
        m.store_discarded(output_file+"_negative.csv", discarded)
        print("Negatives saved")
    else:
        print("Negatives not saved")
        
    if(not raw_stats.empty and not cooked_stats.empty):
        m.store_stats(output_file+"_stats.txt", raw_stats, cooked_stats)
        print("Statistics saved")
    else:
        print("Statistics not saved")
    """       
    except:
        print("Unable to save files")
    """
    
    return

"""
Loading input stuff
"""

def load_gen_record(gen_file, check_integrity, check_uppercase, hanging_primers):
    gen_record = None
    try:
        gen_record = ld.load_bio_files(gen_file, writable=check_integrity or check_uppercase, 
                                       check_uppercase=check_uppercase)
        logging.info("Genome record file loaded!")
    except:
        logging.critical("Error at loading gen record")
    if(check_integrity): gen_record = ld.remove_bad_gens(gen_record)
    return gen_record

def load_primer_pairs(primer_pairs_file):
    primer_pairs = None
    try:
        primer_pairs = ld.load_csv_file(primer_pairs_file)
        logging.info("Primer pairs file loaded!")
    except:
        logging.critical("Error at loading primer pairs file")
        
    return primer_pairs

"""
Loading template stuff
"""
def load_template(parameters, thread_q=None):
    
    gen_record = parameters.loc["gen", "value"]
    primer_pairs = parameters.loc["primer_pairs", "value"]
                        
    if(gen_record  and primer_pairs):
        gen_record = load_gen_record(gen_record , parameters.loc["check_integrity", "value"], 
                                     parameters.loc["check_uppercase", "value"], parameters.loc["hanging primers", "value"])
        primer_pairs = load_primer_pairs(primer_pairs)
    try:
        template = ld.load_template(parameters.loc["csv_template", "value"])
        max_misses = int(parameters.loc["forward missmatches", "value"]) + int(parameters.loc["reverse missmatches", "value"])
        template, discarded, raw_stats, cooked_stats = ld.restore_template(template, gen_record, primer_pairs, max_misses)
        print("Template file restored!")
        
        if(thread_q!=None):
            thread_q.put(template)
            thread_q.put(discarded)
            thread_q.put(raw_stats)
            thread_q.put(cooked_stats)
            thread_q.put(gen_record)
            thread_q.put(primer_pairs)
        return template, discarded, gen_record, primer_pairs, raw_stats, cooked_stats

    except:
        logging.error("Unable to restore template")
        if(thread_q!=None): #TODO Ugly code, if repeated twice
            thread_q.put(pd.DataFrame())
            thread_q.put(pd.DataFrame())
            thread_q.put(pd.DataFrame())
            thread_q.put(pd.DataFrame())
            thread_q.put(gen_record)
            thread_q.put(primer_pairs)
        return pd.DataFrame(), pd.DataFrame(), gen_record, primer_pairs, pd.DataFrame(), pd.DataFrame()

def load_template_only(template_file):
    template = pd.DataFrame();
    try:
        template = ld.load_template(template_file)
        print("Template loaded!")
    except:
        logging.error("Unable to load template")
    return template

def get_Nend_match(template, nend, max_misses, thread_q=None):
    nend_template, raw_stats, cooked_stats =  m.get_Nend_template(template, nend, max_misses)
    if(thread_q!=None):
            thread_q.put(nend_template)
            thread_q.put(raw_stats)
            thread_q.put(cooked_stats)
    
    return nend_template, raw_stats, cooked_stats
        
"""
Simulation stuff
"""
"""
def simulate(template, sample_size, k, B, N):
    sim = s.Simulation(template, sample_size)
    raw, cooked = sim.simulate(k, B, N)
    return raw, cooked
"""

