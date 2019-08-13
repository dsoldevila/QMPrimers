#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 17:02:08 2018

@author: david
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from common import *

import numpy as np
import pandas as pd
import load_data as ld


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

def compute_primer_pair_best_alignment(max_miss_f, max_miss_r, primer, gen, hanging_primers, template, discarded, alignment_processor):
    """
    Returns the best alignments between a genome and a primer pair
    @returns: PrimerAlignment instance
    """
    logging.debug("pp "+str(primer.id))
    
    max_amplicon = primer.max_amplicon
    search_limit = primer.rlen+max_amplicon
    len_gen = len(gen)
    if(primer.flen+search_limit>len_gen): #If primer pair plus max_amplicon is larger than the genomic sequence, check if with min_amplicon the same happens
        if(primer.flen+primer.rlen+primer.min_amplicon>len_gen): #If primer pair plus min_amplicon is larger than the genomic sequence, abort
            logging.warning("Skipping gen "+gen.id+" primer pair "+str(primer.id))
            discarded.loc[discarded.shape[0]] = [primer.id, gen.id]
            alignment_processor.add_negative_2_stats(primer.id)
            return template, discarded
        else: #else modify max_amplicon to keep the primer_pair within the limits
            max_amplicon = len_gen - (primer.flen + primer.rlen)
            search_limit = primer.rlen+max_amplicon
            logging.debug("Cutting max amplicon from "+str(primer.max_amplicon)+" to "+str(max_amplicon))
    
    best_score = 0
    alignments = []

    forward_matchings = _compute_primer_matching(max_miss_f, primer.f.seq, primer.flen, gen.seq[0:-search_limit]) #compute forward primer best matches
    for fm in forward_matchings: #for each match with forward primer, compute reverse matchings
        logging.debug("forward matching: "+str(fm))
        start = fm[2]+primer.min_amplicon #forward match start + len(forward) + min amplicon
        end = fm[2]+max_amplicon+primer.rlen #f match start + len(f) + max amplicon + len(r)
        logging.debug("Reverse match search from "+str(start)+" to "+str(end))
        reverse_matchings = _compute_primer_matching(max_miss_r, primer.r.seq, primer.rlen, gen.seq[start:end])
        for rm in reverse_matchings: #get the best or bests matche(s) with this primer pair (alingments)
            score = fm[0] + rm[0]
            if (score > best_score): #if the score is better, erase the previous bests results
                alignments = [(fm, rm)]
                best_score = score
            elif (score == best_score): #elif the score is equaly good, get this alignment too
                alignments.append((fm, rm))
                
    if(alignments==[]):
        logging.debug(" No alignment")
        discarded.loc[discarded.shape[0]] = [primer.id, gen.id]
        alignment_processor.add_negative_2_stats(primer.id)
    for al in alignments:
        fm = al[0]
        rm= al[1]
        amplicon = primer.min_amplicon+rm[1]
        alignment_processor.get(gen, primer, fm[1], fm[1]-max_miss_f*hanging_primers, amplicon+fm[2], amplicon+fm[2]-max_miss_f*hanging_primers,
                                primer.flen-fm[0], primer.rlen-rm[0], amplicon) #TODO, creating a temp class overkill?
        template.loc[template.shape[0]] = alignment_processor.get_csv()
        
    return template, discarded

def compute_gen_matching(max_miss_f, max_miss_r, primer_pairs, gen_record, output_file, hanging_primers=False):
    """
    Computes the best alignments between each genome and each primer
    @returns: List of GenAlignment instances
    """
    try:
        assert(max_miss_f>0), "Max forward misses must be greater than 0"
        assert(max_miss_r>0), "Max reverse misses must be greater than 0"
        #assert other params or let it crash?
        os.makedirs(os.path.dirname(output_file), exist_ok=True) #make sure output file exists
    except(Exception) as e:
        logging.error(e)
        return pd.DataFrame(), pd.DataFrame(),pd.DataFrame(), pd.DataFrame()
    logging.debug("Testing debug")
    if(hanging_primers):
        gen_record = append_zeros(gen_record, max_miss_f, max_miss_r)
    
    primerPair_list = [] #used to sort template by primer pairs
    for pkey in primer_pairs:
        pp = primer_pairs[pkey]
        pp.f.seq = np.array(pp.f)
        pp.r.seq = np.array(pp.r)
        primerPair_list.append(pkey)
    
    size = len(gen_record)
    i = 0       
    template_header = TEMPLATE_HEADER

    template = pd.DataFrame(columns=template_header)
    
    discarded = pd.DataFrame(columns=TEMPLATE_HEADER[0:2])

    alignment_processor = Alignment(max_miss_f + max_miss_r)
    
    template.to_csv(output_file+"_positive.csv", index_label="id")
    discarded.to_csv(output_file+"_negative.csv", index_label="id")
    
    t1 = 0
    d1 = 0
    for gen_key in gen_record:
        print(gen_key, "{0:.2f}".format(i/size*100)+"%")
        i +=1
        gen = gen_record[gen_key]
        gen.seq = np.array(gen.seq)
        for pkey in primer_pairs:
            try:
                t2 = t1
                d2 = d1
                template, discarded = compute_primer_pair_best_alignment(max_miss_f, max_miss_r, primer_pairs[pkey], gen, hanging_primers, template, discarded, alignment_processor)
                t1 = template.shape[0]
                d1 = discarded.shape[0]
                template.loc[template.index[t2:t1]].to_csv(output_file+"_positive.csv", mode='a', index_label="id", header=None)
                discarded.loc[discarded.index[d2:d1]].to_csv(output_file+"_negative.csv", mode='a', index_label="id", header=None)
            except:
                raise
                logging.error("Skipping gen "+gen.id+" primer pair "+str(pkey))
   
    raw_stats, cooked_stats = alignment_processor.get_stats()
    #store_stats("Temp file", output_file+"_stats.txt", raw_stats, cooked_stats)
    #print("Template, negative and statistics saved")
    print("Matching done!")
    
    template["primerPair"] = pd.Categorical(template["primerPair"], categories=primerPair_list, ordered=True)
    template.sort_values(['primerPair', 'fastaid'], inplace=True)
    template.reset_index(drop=True, inplace=True)
    
    discarded["primerPair"] = pd.Categorical(discarded["primerPair"], categories=primerPair_list, ordered=True)
    discarded.sort_values(['primerPair', 'fastaid'], inplace=True)
    discarded.reset_index(drop=True, inplace=True)
    return template, discarded, raw_stats, cooked_stats

"""
OTHER FUNCTIONS
"""

def store_matching_results(input_files, output_file, template, header):
    """
    Stores alignment results
    @param gen_matching_list list of GenMatching instances
    @return None
    """
    try:
        columns = template.columns.values
        h = []
        for i in range(len(header)):
            h.append(columns[header[i]])
        with open(output_file,'w') as outfile:
                outfile.write(input_files+"\n")
                outfile.write(str(datetime.datetime.now())+"\n")
                template.to_csv(outfile, index_label="id", columns=h)
        print("Template saved!")
    except(Exception) as e:
        logging.error(e)
    
    return

def store_stats(input_files, output_file, raw_stats, cooked_stats):
    try:
        with open(output_file,'w') as outfile:
            outfile.write(input_files+"\n")
            outfile.write(str(datetime.datetime.now())+"\n")
            raw_stats.to_string(outfile)
            outfile.write("\n\n")
            cooked_stats.to_string(outfile)
            print("Statistics saved")
    except(Exception) as e:
        logging.error(e)
    return

def store_discarded(input_files, output_file, discarded):
    try:
        with open(output_file,'w') as outfile:
            outfile.write(input_files+"\n")
            outfile.write(str(datetime.datetime.now())+"\n")
            discarded.to_csv(outfile, index_label="id")
            print("Negatives saved")
    except(Exception) as e:
        logging.error(e)
    return


def get_Nend_template(template, nend):
    """
    @brief With the computed template, generate a template but with Nend mismatches
    @Return new template, raw stats, cooked_stats
    """
    try:
        assert(not template.empty), "Bad template"
        assert(nend>0), "N-end must be greater than 0"
    except(Exception) as e:
        logging.error(e)
        return
    mismF = get_missmatch_column_name(template.columns.values, primer="f")
    mismR = get_missmatch_column_name(template.columns.values, primer="r")
    
    mismF = template[mismF].max()
    mismR = template[mismR].max()
    
    alignment = Alignment(mismF+mismR);
    
    #TODO this dupplicates the TEMPLATE_HEADER. This header should be generated from TEMPLATE_HEADER.
    header = ["primerPair","fastaid","primerF","primerR","mismFN"+str(nend),"mismRN"+str(nend),"ampliconLen", "F_pos", 
              "mismFN"+str(nend)+"_loc", "mismFN"+str(nend)+"_type", "mismFN"+str(nend)+"_base", "R_pos", "mismRN"+str(nend)+"_loc",
              "mismRN"+str(nend)+"_type", "mismRN"+str(nend)+"_base", "amplicon"]
    
    nend_template = pd.DataFrame(columns=header)
    size = template.shape[0]
    for i in range(template.shape[0]):
        nend_template.loc[i] = alignment.get_Nend(*template.loc[i], nend)
        print("Getting Nend "+"{0:.2f}".format(i/size*100)+"%")
    raw_stats, cooked_stats = alignment.get_stats()
    
    return nend_template, raw_stats, cooked_stats


def debug_matching(gen, primer_pair, mf, mr, output_file, hanging_primers=False):
    """
    This function computes and displays a single alignment. Used for debugging purposes
    """
    try:
        assert(len(gen)==1), "Multiple gen sequences detected"
        assert(len(primer_pair)==1), "Multiple primer pairs detected"
    except(Exception) as e:
        logging.error(e)
        return
    
    template, discarded, raw_stats, cooked_stats = compute_gen_matching(mf, mr, primer_pair, gen, output_file, hanging_primers=hanging_primers)
    
    if(template.empty):
        logging.warning("No result")
        return
    
    match_result = template.loc[0]
    pp = primer_pair[next(iter(primer_pair))]
    gen = gen[next(iter(gen))]
    
    fpos = match_result.at['F_pos'] -1
    
    rem_len = len(gen)-(fpos+pp.flen+match_result.at['ampliconLen']+pp.rlen)
    
    pp.f.seq = Seq(''.join(pp.f.seq))
    pp.r.seq = Seq(''.join(pp.r.seq))
    
    pp_aligned = '-'*fpos+pp.f.seq+'-'*match_result.at['ampliconLen']+pp.r.seq+'-'*rem_len
    pp_aligned=SeqRecord(pp_aligned)
    pp_aligned.id = pp.id
    align = MultipleSeqAlignment([gen, pp_aligned])
    
    print(align.format("clustal"))
    
    try:
        with open(output_file+".txt",'w') as outfile:
            outfile.write(align.format("clustal"))
            print("Debug saved")
    except(Exception) as e:
        logging.error(e)
        
    return
        
    

if(__name__=="__main__"):
    
    gen_record = ld.load_bio_files(["Data/species_bold_own_genbank.fasta"])
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")


    import time 
    time1 = time.time()
    result = compute_gen_matching(5, 5, primer_pairs, gen_record)
    elapsedTime = ((time.time()-time1))
    print(int(elapsedTime))

    