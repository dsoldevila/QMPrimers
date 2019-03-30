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
    gen_record = load_gen_record(parameters["gen"], writable=parameters["hanging primers"], check_integrity=True)
    primer_pairs = load_primer_pairs(parameters["primer_pairs"])
    if(gen_record!=None and primer_pairs!=None):
        template = m.compute_gen_matching(int(parameters["forward missmatches"]), int(parameters["reverse missmatches"]), primer_pairs, gen_record, hanging_primers=parameters["hanging primers"])

    return template

def save_template_primer_missmatches(output_file, template, header=None):
    m.store_matching_results(output_file, template, header)
    return

def load_gen_record(gen_file, writable=False, check_integrity=False, file_format=None):
    gen_record = None
    try:
        gen_record = ld.load_bio_files(gen_file, writable=writable, file_format=file_format)
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

def simulate():
    return