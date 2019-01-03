#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 17:02:08 2018

@author: david
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from common import *

import numpy as np
import pandas as pd

import load_data as ld

SCORE_TABLE = np.array([[1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0],
                        [0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0],
                        [0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0],
                        [0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0],
                        [1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                        [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                        [1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                        [1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                        [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype='uint8')
MATCH_TABLE = pd.DataFrame(SCORE_TABLE, index=list("ACGTWSMKRYBDHVNIZ"), columns=list("ACGTWSMKRYBDHVNIZ"))

def is_valid(score, min_score):
    score = np.asarray(score)
    min_score = np.asarray(min_score)
    return np.where(score >= min_score, score, 0)

def string2numpy(string):
   return np.array(list(string))

def _compute_primer_matching(max_misses, primer, len_primer, gen):
    result_matrix = MATCH_TABLE.loc[primer,gen] #get match table
    result_max_len = len(gen)-len_primer+1
    result_raw = np.zeros(result_max_len, dtype='i4') #TODO integer 32 too much?
    """
    offset=0
    dooooooo
    odoooooo
    oodooooo
    ooodoooo...
    """
    #TODO: Vectorize this for
    for i in range(result_max_len):
        result_raw[i] = np.sum((np.diagonal(result_matrix, offset=i))) #get scores of all diagonals
        
    is_score_valid = (result_raw>=len_primer-max_misses)
    n_results = is_score_valid.sum()
    result_raw = is_score_valid*result_raw #if score valid =score else =0
    
    result = np.zeros(n_results, dtype=[('score', '>i4'), ('start', '>i4'), ('end', '>i4')]) #TODO integer 32 too much?
    
    j = 0
    for i in range(result_max_len):
        #result[j] = (result_raw[i], i, i+len_primer) if result_raw[i] else result[j]
        if result_raw[i]:
            result[j] = (result_raw[i], i, i+len_primer)
            j=j+1
    return result

def compute_primer_pair_best_alignment(max_miss_f, max_miss_r, primer, gen):
    search_limit = primer.rlen+primer.max_amplicon
    if(primer.flen+search_limit>len(gen)): #If primer pair (plus amplicon) is larger than the genomic sequence, abort
        return None
    
    best_score = 0
    
    alignments = []
    result = MatchingList(gen)
    
    forward = string2numpy(primer.f)
    forward_matchings = _compute_primer_matching(max_miss_f, forward, primer.flen, gen[0:-search_limit]) #compute forward primer matching
    
    reverse = string2numpy(primer.r)
    for fm in forward_matchings: #for each matching with forward primer, compute reverse matchings
        ralignments = []
        start = fm[2]+primer.min_amplicon #forward match start + len(forward) + min amplicon
        end = fm[2]+primer.max_amplicon+primer.rlen #f match start + len(f) + max amplicon + len(r)
        reverse_matchings = _compute_primer_matching(max_miss_r, reverse, primer.rlen, gen[start:end])
        for rm in reverse_matchings: #get the best or bests matche(s) with this primer pair (alingments)
            score = fm[0] + rm[0]
            if (score > best_score): #if the score is better, erase the previous bests results
                alignments = [(fm, rm)]
                best_score = score
            elif (score == best_score):
                ralignments.append((fm, rm))
                
                
        for al in alignments:
            fm = al[0]
            rm= al[1]
            amplicon = primer.min_amplicon+rm[0]
            result.append(Matching(gen, primer, fm[1], start+rm[1], primer.flen-fm[0], primer.rlen-rm[0], amplicon, MATCH_TABLE))
            
    return result

def compute_gen_alignments(max_miss_f, max_miss_r, primer_pairs, gen):
    result = MatchingList(gen)
    for primer in primer_pairs:
        result.append(compute_primer_pair_best_alignment(max_miss_f, max_miss_r, primer, gen))
    pass

if(__name__=="__main__"):
    gen_record = ld.load_bio_file("Data/species_bold_own_genbank.fasta")
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")
    gen = gen_record.get("ACEA563-14_Aphis_gossypii_BOLD")
    primer = primer_pairs[4]
    result = compute_primer_pair_best_alignment(5, 5, primer, gen)
    
    print(result)