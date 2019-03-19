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

#This matrix tells the algorithm whether 2 nucleotides match or don't
SCORE_TABLE = np.array([[1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0],
                        [0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0],
                        [0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0],
                        [0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0],
                        [1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0],
                        [0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0],
                        [0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0],
                        [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
                        [1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0],
                        [1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0],
                        [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype='uint8')
MATCH_TABLE = pd.DataFrame(SCORE_TABLE, index=list("ACGTWSMKRYBDHVNIZ"), columns=list("ACGTWSMKRYBDHVNIZ"))

def _compute_matching_chunk(primer, len_primer, gen, start, end, matching_matrix, result_raw):
    """
    @brief: Computes matching matrix[start:end], the matching is done in chunks to improve memory efficiency
    """
    rstart = start+len_primer-1
    rend = end+len_primer-1

    matching_matrix[0:len_primer, 0:end-start] = MATCH_TABLE.loc[primer, gen[start:end]].values
    for i in range(len_primer):
        #matching_matrix[i, 0:end-start] = gen[start:end]==primer[i] #TODO
        result_raw[rstart-i:rend-i] = np.add(result_raw[rstart-i:rend-i], matching_matrix[i, 0:end-start])
        
    return

def _compute_matching_chunk_debug(primer, len_primer, gen, start, end, matching_matrix, result_raw):
    """
    @brief: Computes matching matrix[start:end], the matching is done in chunks to improve memory efficiency
    """
    rstart = start+len_primer-1
    rend = end+len_primer-1

    matching_matrix[0:len_primer, 0:end-start] = MATCH_TABLE.loc[primer, gen[start:end]].values
    for i in range(len_primer):
        #matching_matrix[i, 0:end-start] = gen[start:end]==primer[i] #TODO
        print("a")
        a= result_raw[rstart-i:rend-i]
        b=  matching_matrix[i, 0:end-start]
        result_raw[rstart-i:rend-i] = np.add(result_raw[rstart-i:rend-i], matching_matrix[i, 0:end-start])
        
    return

def _compute_forward_matching(max_misses, primer_pairs, gen, chunk_size, max_row_size, data_type='uint8'):
    """
    @brief: Computes all forward primers at once for a given genomic sequence
    """
    start = 0
    end = chunk_size
    #len_gen = gen.shape[0]
    len_gen = len(gen)
    if end > len_gen or end < start:
        end = len_gen
        
        
    matching_matrix =  np.zeros((max_row_size, chunk_size), dtype=data_type) #matrix buffer to compute matches
    rr_dict = {} #raw_result dictionary, raw_result is a buffer that stores a matching score for every fprimer position
    
    for pp in primer_pairs.values():
        rr_dict[pp.id]= (np.zeros((len_gen+pp.flen-1,), dtype=data_type))
    
    while(end <= len_gen):
        for pp in primer_pairs.values(): 
            _compute_matching_chunk(pp.f, pp.flen, gen, start, end, matching_matrix, rr_dict[pp.id])
        start = end
        end += chunk_size
    
    if(len_gen % chunk_size): #if there is "gen left", compute the rest
        for pp in primer_pairs.values():
            _compute_matching_chunk(pp.f, pp.flen, gen, start, len_gen, matching_matrix, rr_dict[pp.id])
     
    return _process_forward_matching(rr_dict, primer_pairs, max_misses, data_type)

def _process_forward_matching(rr_dict, primer_pairs, max_misses, data_type):
    """
    @brief: Transforms _compute_forward_matching result into something readable
    """
    pr_dict = {} # processed result dictionary
    
    for pp in primer_pairs.values():
        flen = pp.flen
        raw_result = rr_dict[pp.id]
        is_result_valid = raw_result >= (flen-max_misses)
        n_results = is_result_valid.sum()
        raw_result = is_result_valid*raw_result
        pro_result = np.empty(n_results, dtype=[('fscore', data_type), ('start', '>i4'), ('end', '>i4'), ('reverse_search', 'object')]) #TODO Use a pandas Dataframe instead?
        
        j=0
        for i in range(len(raw_result)):        
            if raw_result[i]:
                pro_result[j] = (raw_result[i], i-flen+1, i, None)
                j=j+1
        pr_dict[pp.id] = pro_result
        
    return  pr_dict


def _compute_reverse_matching(max_misses, primer_pairs, gen, max_row_size, pr_dict, data_type='uint8'):
    """
    @brief: Computes all forward primers at once for a given genomic sequence, given forward results
    """
    
    #Set Job queue
    #For each fmatch a reverse matching job is set and appended to it. Matching_request is an ordered list of all
    #reverse matching jobs to improve memory efficiency when computing
    matching_request = []
    max_diff = 0
    for pp in primer_pairs.values(): # for each primer_pair in primer_pairs
        min_range = pp.min_amplicon #range to search relative to the forward primer
        max_range = pp.max_amplicon+pp.rlen
        diff = max_range-min_range
        max_diff = diff if diff > max_diff else max_diff
        for fmatch in pr_dict[pp.id]: #for each processed result with primer_pair i
            start = min_range + fmatch[2] #range to search in gen
            end = max_range + fmatch[2]
            result_raw = np.zeros(((end-start)+pp.rlen-1,), dtype=data_type)
            fmatch[3] = (pp.id, start, end, result_raw)
            matching_request.append(fmatch[3])
    
    matching_request = sorted(matching_request, key=lambda req: req[1])
    
    matching_matrix =  np.zeros((max_row_size, max_diff), dtype=data_type) #matrix buffer to compute matches
    
    #Compute
    for mr in matching_request:
        pp = primer_pairs[mr[0]]
        _compute_matching_chunk_debug(pp.r, pp.rlen, gen, mr[1], mr[2], matching_matrix, mr[3])
        
    for pp in primer_pairs.values():
        for fmatch in pr_dict[pp.id]:
            rresult = fmatch[3]
            rresult = _process_reverse_match(rresult, pp, max_misses, data_type)
    
    return pr_dict

def _process_reverse_match(rresult, pp, max_misses, data_type):
    """
    @brief: Transforms a reverse (matching) result into something human readable
    """
    # rraw = (pp.id, start, end, result_raw)
    raw_result = rresult[3]
    start = rresult[1]
    rlen = pp.rlen
    is_result_valid = raw_result >= (rlen-max_misses)
    n_results = is_result_valid.sum()
    raw_result = is_result_valid*raw_result
    pro_result = np.empty(n_results, dtype=[('rscore', data_type), ('start', '>i4'), ('end', '>i4')]) #TODO Use a pandas Dataframe instead?
    
    j=0
    for i in range(len(raw_result)):        
        if raw_result[i]:
            pro_result[j] = (raw_result[i], i-rlen+1+start, i+start)
            j=j+1
    rresult[3] = pro_result
    
    return rresult

def get_gen_alignment(primer_pairs, gen, pr_dict):
    
    gen_alignment = GenAlignment(gen)
    for pp in primer_pairs.values():
        best_score = 0
        alignments = []
        for fmatch in pr_dict[pp.id]:
            reverse_matches = fmatch[3]
            for rmatch in reverse_matches:
                score = fmatch[0] + rmatch[0]
                if (score > best_score): #if the score is better, erase the previous bests results
                    alignments = [(fmatch, rmatch)]
                    best_score = score
                elif (score == best_score): #elif the score is equaly good, get this alignment too
                    alignments.append((fmatch, rmatch))
        
        primer_alignment = []          
        for al in alignments:
            fmatch = al[0]
            rmatch= al[1]
            amplicon = pp.min_amplicon+rmatch[1]
            primer_alignment.append(Alignment(gen, pp, fmatch[1], fmatch[1], rmatch[1]+amplicon+fmatch[2], rmatch[1]+amplicon+fmatch[2],
                                    pp.flen-fmatch[0], pp.rlen-rmatch[0], amplicon, MATCH_TABLE))
        gen_alignment.append(primer_alignment)
            
    return gen_alignment
                

def compute_matching(max_misses_f, max_misses_r, primer_pairs, gen_record, chunk_size, data_type='uint8'):
    gen_matching_list = []
    for gkey in gen_record:
        gen = gen_record[gkey]
        forward_matches = _compute_forward_matching(max_misses_f, primer_pairs, gen, chunk_size,30, data_type='uint8')
        full_matches = _compute_reverse_matching(max_misses_r, primer_pairs, gen, 30, forward_matches, data_type='uint8')
        gen_alignment = get_gen_alignment(primer_pairs, gen, full_matches)
        gen_matching_list.append(gen_alignment)
    
    return

if(__name__=="__main__"):
    
    gen_record = ld.load_bio_files(["Data/species_bold_own_genbank.fasta"])
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")

    #result = compute_gen_matching(5, 5, primer_pairs, gen_record)
    
    import time 
    time1 = time.time()
    result = compute_matching(5, 5, primer_pairs, gen_record)
    elapsedTime = ((time.time()-time1))
    print(int(elapsedTime))