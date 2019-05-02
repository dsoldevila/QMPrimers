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


#TODO this func is currently unused, ~line 68
def is_valid(score, min_score):
    score = np.asarray(score)
    min_score = np.asarray(min_score)
    return np.where(score >= min_score, score, 0)

def append_zeros(gen_record, max_miss_f, max_miss_r):
    """
    Adds zeros ('Z') at the beginning and end of the genome to allow the primer to pass the genome's limits
    Example:
        genome AGATTCATT, primer TCTAGA scores 2 at pos 1
        genome ZZZAGATTCATTZZZ, primer TCTAGA scores 3 at pos -3
    """
    for gen_key in gen_record:
        gen_record[gen_key].seq = Seq("Z"*max_miss_f+str(gen_record[gen_key].seq)+"Z"*max_miss_r)
    return gen_record

def _compute_primer_matching(max_misses, primer, len_primer, gen):
    """
    Computes the best matches between a genome and a primer.
    @returns: Numpy matrix of arrays (score, start_pos, end_pos).
    """
    result_matrix = MATCH_TABLE.loc[primer, gen] #get match table

    result_max_len = len(gen)-len_primer+1
    result_raw = np.zeros(result_max_len, dtype='uint8') #TODO 0-255 should be enough, but better to not hardcode this

    result_matrix = result_matrix.values #get numpy matrix, raw indexes are faster than labels
    for i in range(len_primer):
        result_raw = np.add(result_raw, result_matrix[i,i:result_max_len+i]) #Speedup try n2, SpeedUp = 168/26 = 6,4x
    
    is_score_valid = (result_raw>=len_primer-max_misses)
    n_results = is_score_valid.sum()
    result_raw = is_score_valid*result_raw #if score valid =score else =0
    
    result = np.zeros(n_results, dtype=[('score', 'uint8'), ('start', '>i4'), ('end', '>i4')]) #TODO integer 32 too much /  integer 8 enough?
    
    j = 0
    for i in range(result_max_len):        
        if result_raw[i]:
            result[j] = (result_raw[i], i, i+len_primer)
            j=j+1

    return result

def compute_primer_pair_best_alignment(max_miss_f, max_miss_r, primer, gen, hanging_primers, template, discarded, Nend, alignment_processor):
    """
    Returns the best alignments between a genome and a primer pair
    @returns: PrimerAlignment instance
    """
    
    max_amplicon = primer.max_amplicon
    search_limit = primer.rlen+max_amplicon
    len_gen = len(gen)
    if(primer.flen+search_limit>len_gen): #If primer pair plus max_amplicon is larger than the genomic sequence, check if with min_amplicon the same happens
        if(primer.flen+primer.rlen+primer.min_amplicon>len_gen): #If primer pair plus min_amplicon is larger than the genomic sequence, abort
            print("Warning: Skipping gen "+gen.id+" primer pair "+str(primer.id))
            discarded.loc[discarded.shape[0]] = [primer.id, gen.id]
            alignment_processor.add_negative_2_stats(primer.id)
            return template, discarded
        else: #else modify max_amplicon to keep the primer_pair within the limits
            max_amplicon = len_gen - (primer.flen + primer.rlen)
    
    best_score = 0
    alignments = []

    forward_matchings = _compute_primer_matching(max_miss_f, primer.f.seq, primer.flen, gen.seq[0:-search_limit]) #compute forward primer best matches
    for fm in forward_matchings: #for each match with forward primer, compute reverse matchings
        start = fm[2]+primer.min_amplicon #forward match start + len(forward) + min amplicon
        end = fm[2]+max_amplicon+primer.rlen #f match start + len(f) + max amplicon + len(r)
        reverse_matchings = _compute_primer_matching(max_miss_r, primer.r.seq, primer.rlen, gen.seq[start:end])
        for rm in reverse_matchings: #get the best or bests matche(s) with this primer pair (alingments)
            score = fm[0] + rm[0]
            if (score > best_score): #if the score is better, erase the previous bests results
                alignments = [(fm, rm)]
                best_score = score
            elif (score == best_score): #elif the score is equaly good, get this alignment too
                alignments.append((fm, rm))
                
    if(alignments==[]):
        discarded.loc[discarded.shape[0]] = [primer.id, gen.id]
        alignment_processor.add_negative_2_stats(primer.id)
    for al in alignments:
        fm = al[0]
        rm= al[1]
        amplicon = primer.min_amplicon+rm[1]
        alignment_processor.get(gen, primer, fm[1], fm[1]-max_miss_f*hanging_primers, rm[1]+amplicon+fm[2], rm[1]+amplicon+fm[2]-max_miss_f*hanging_primers,
                                primer.flen-fm[0], primer.rlen-rm[0], amplicon, Nend) #TODO, creating a temp class overkill?
        template.loc[template.shape[0]] = alignment_processor.get_csv()
        
    return template, discarded

def compute_gen_matching(max_miss_f, max_miss_r, primer_pairs, gen_record, Nend, hanging_primers=False):
    """
    Computes the best alignments between each genome and each primer
    @returns: List of GenAlignment instances
    """
    if(hanging_primers):
        gen_record = append_zeros(gen_record, max_miss_f, max_miss_r)
    
    for pkey in primer_pairs:
        pp = primer_pairs[pkey]
        pp.f.seq = np.array(pp.f)
        pp.r.seq = np.array(pp.r)
        
    size = len(gen_record)
    i = 0       
    template_header = TEMPLATE_HEADER
    if(Nend):
        template_header.extend(["mismFN"+str(Nend), "mismRN"+str(Nend)])
    template = pd.DataFrame(columns=template_header)
    
    discarded = pd.DataFrame(columns=TEMPLATE_HEADER[0:2])

    alignment_processor = Alignment(max_miss_f + max_miss_r)
    
    for gen_key in gen_record:
        print(gen_key, "{0:.2f}".format(i/size*100)+"%")
        i +=1
        gen = gen_record[gen_key]
        gen.seq = np.array(gen.seq)
        for pkey in primer_pairs:
            try:
                template, discarded = compute_primer_pair_best_alignment(max_miss_f, max_miss_r, primer_pairs[pkey], gen, hanging_primers, template, discarded, Nend, alignment_processor)
            except:
                raise
                print("Error: Skipping gen "+gen.id+" primer pair "+str(pkey))

    raw_stats, cooked_stats = alignment_processor.get_stats()
    
    return template, discarded, raw_stats, cooked_stats

def store_matching_results(output_file, template, header=None):
    """
    Stores alignment results
    @param gen_matching_list list of GenMatching instances
    @return None
    """
    template.to_csv(output_file, columns=header, index_label="id");
    
    return

def store_stats(output_file, raw_stats, cooked_stats):
    with open(output_file,'w') as outfile:
        raw_stats.to_string(outfile)
        outfile.write("\n\n")
        cooked_stats.to_string(outfile)
    return

def store_discarded(output_file, discarded):
    discarded.to_csv(output_file, index_label="id")
    return

if(__name__=="__main__"):
    
    gen_record = ld.load_bio_files(["Data/species_bold_own_genbank.fasta"])
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")


    import time 
    time1 = time.time()
    result = compute_gen_matching(5, 5, primer_pairs, gen_record)
    elapsedTime = ((time.time()-time1))
    print(int(elapsedTime))

    