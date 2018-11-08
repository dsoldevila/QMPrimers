#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 15:34:21 2018

@author: david
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def _pcr_amp(max_miss_f, max_miss_r, primer_tuple, gen_seq, model):
    
    pass

def compute_matching(max_miss_f, max_miss_r, primer_pairs, gen_seqs, model):
    """
    @parameters:
        max_miss_f: Maximum number of missmatches allowed in forward primer 
        max_miss_r: Maximum number of missmatches allowed in reverse primer
        primer_pairs: list of primer pairs
        gen_seqs: genomic sequences
        model: Using model 1 or 2
        amplicon: required distance between primers
    @return:
    """
    for gen in gen_seqs:
        for primer_tuple in primer_pairs:
            _pcr_amp(max_miss_f, max_miss_r, primer_tuple, gen_seqs.get(gen).seq, model)
    
    print(gen_seqs.get("ACEA1016-14_Aphis_spiraecola_BOLD"))
            
    return
