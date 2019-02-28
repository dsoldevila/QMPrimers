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


def compute_forward_matching(max_misses, primer_pairs, gen, chunk_size, max_row_size, data_type='uint8'):
    start = 0
    end = chunk_size
    len_gen = gen.shape[0]
    if end > len_gen or end < start:
        end = len_gen
        
        
    matching_matrix =  np.zeros((max_row_size, chunk_size), dtype=data_type) #matrix buffer to compute matches
    rr_dict = {} #raw_result dictionary
    
    for pp in primer_pairs:
        rr_dict[pp.id]= (np.zeros((len_gen+pp.flen-1,), dtype=data_type))
    
    while(end <= len_gen):
        for pp in primer_pairs:
            _compute_matching_chunk(pp.f, pp.flen, gen, start, end, matching_matrix, rr_dict[pp.id])
        start = end
        end += chunk_size
    
    if(len_gen % chunk_size): #if there is "gen left", compute the rest
        for pp in primer_pairs:
            _compute_matching_chunk(pp.f, pp.flen, gen, start, len_gen, matching_matrix, rr_dict[pp.id])
     
    return process_matching(rr_dict, primer_pairs, max_misses, data_type)

def compute_reverse_matching(max_misses, primer_pairs, gen, chunk_size, max_row_size, pr_dict, data_type='uint8'):
    start = 0
    end = chunk_size
    len_gen = gen.shape[0]
    if end > len_gen or end < start:
        end = len_gen
        
    
    return

def compute_matching(max_misses_f, max_misses_r, primer_pairs, gen_record, chunk_size, pr_dict, data_type='uint8'):
    
    return

def process_matching(rr_dict, primer_pairs, max_misses, data_type):
    
    pr_dict = {} # processed result dictionary
    
    for pp in primer_pairs:
        flen = pp.flen
        raw_result = rr_dict[pp.id]
        raw_result = raw_result
        is_result_valid = raw_result >= (flen-max_misses)
        n_results = is_result_valid.sum()
        raw_result = is_result_valid*raw_result
        pro_result = np.empty(n_results, dtype=[('misses', data_type), ('start', '>i4'), ('end', '>i4')])
        
        j=0
        for i in range(len(raw_result)):        
            if raw_result[i]:
                pro_result[j] = (raw_result[i], i-flen+1, i)
                j=j+1
        pr_dict[pp.id] = pro_result
        
    return  pr_dict
    

def compute_primer_pair_best_alignment(max_miss_f, max_miss_r, primer, bio_gen, gen, hanging_primers):
    """
    Returns the best alignments between a genome and a primer pair
    @returns: PrimerAlignment instance
    """
    result = PrimerAlignment(primer, gen)
    max_amplicon = primer.max_amplicon
    search_limit = primer.rlen+max_amplicon
    len_gen = len(gen)
    if(primer.flen+search_limit>len_gen): #If primer pair plus max_amplicon is larger than the genomic sequence, check if with min_amplicon the same happens
        if(primer.flen+primer.rlen+primer.min_amplicon>len_gen): #If primer pair plus min_amplicon is larger than the genomic sequence, abort
            return result
        else: #else modify max_amplicon to keep the primer_pair within the limits
            max_amplicon = len_gen - (primer.flen + primer.rlen)
    
    best_score = 0
    
    alignments = []
    
    forward_matchings = _compute_primer_matching(max_miss_f, primer.f, primer.flen, gen[0:-search_limit]) #compute forward primer best matches
    for fm in forward_matchings: #for each match with forward primer, compute reverse matchings
        start = fm[2]+primer.min_amplicon #forward match start + len(forward) + min amplicon
        end = fm[2]+max_amplicon+primer.rlen #f match start + len(f) + max amplicon + len(r)
        reverse_matchings = _compute_primer_matching(max_miss_r, primer.r, primer.rlen, gen[start:end])
        for rm in reverse_matchings: #get the best or bests matche(s) with this primer pair (alingments)
            score = fm[0] + rm[0]
            if (score > best_score): #if the score is better, erase the previous bests results
                alignments = [(fm, rm)]
                best_score = score
            elif (score == best_score): #elif the score is equaly good, get this alignment too
                alignments.append((fm, rm))
                
                
    for al in alignments:
        fm = al[0]
        rm= al[1]
        amplicon = primer.min_amplicon+rm[1]
        result.append(Alignment(bio_gen, primer, fm[1], fm[1]-max_miss_f*hanging_primers, rm[1]+amplicon+fm[2], rm[1]+amplicon+fm[2]-max_miss_f*hanging_primers,
                                primer.flen-fm[0], primer.rlen-rm[0], amplicon, MATCH_TABLE))
            
    return result

def compute_gen_matching(max_miss_f, max_miss_r, primer_pairs, gen_record, hanging_primers=False):
    """
    Computes the best alignments between each genome and each primer
    @returns: List of GenAlignment instances
    """
    if(hanging_primers):
        gen_record = append_zeros(gen_record, max_miss_f, max_miss_r)
        
    gen_alignment_list = []
    size = len(gen_record)
    i = 0
    for gen_key in gen_record:
        print(gen_key, "{0:.2f}".format(i/size*100)+"%")
        i +=1
        gen = gen_record[gen_key]
        gen_alignment = GenAlignment(gen)
        numpy_gen = np.array(list(gen.seq))
        for primer in primer_pairs:
            alignment_list = compute_primer_pair_best_alignment(max_miss_f, max_miss_r, primer, gen, numpy_gen, hanging_primers)
            gen_alignment.append(alignment_list)
        gen_alignment_list.append(gen_alignment)
    return gen_alignment_list

if(__name__=="__main__"):
    
    gen_record = ld.load_bio_files(["Data/species_bold_own_genbank.fasta"])
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")

    #result = compute_gen_matching(5, 5, primer_pairs, gen_record)
    
    import time 
    time1 = time.time()
    result = compute_gen_matching(5, 5, primer_pairs, gen_record)
    elapsedTime = ((time.time()-time1))
    print(int(elapsedTime))