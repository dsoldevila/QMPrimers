#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#DEPRECATED
"""
Created on Thu Nov  8 15:34:21 2018

@author: David Soldevila
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from common import *

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Alphabet import IUPAC

import load_data as ld
import custom_pairwise2

"""
Important reminder:
    THe first pos of a pairwise2 alignment is the first with match:
        AAG and TAG will return a score of 2 starting at 1
    The last pos of parwise2 alignment CAN be the last with a match
        TAT and TAG will return two identical matches, with ending positions of 2 and 1, it seems that the one with higher end position has order priority.
"""
#IUPAC Notation
MATCH_TABLE = {('A', 'A'): 1, ('A', 'C'): 0, ('A', 'G'): 0, ('A', 'T'): 0, ('A', 'W'): 1, ('A', 'S'): 0, ('A', 'M'): 1, ('A', 'I'): 1,
               ('A', 'K'): 0, ('A', 'R'): 1, ('A', 'Y'): 0, ('A', 'B'): 0, ('A', 'D'): 1, ('A', 'H'): 1, ('A', 'V'): 1, ('A', 'N'): 0, ('A', 'Z'): 0,
               ('C', 'A'): 0, ('C', 'C'): 1, ('C', 'G'): 0, ('C', 'T'): 0, ('C', 'W'): 0, ('C', 'S'): 1, ('C', 'M'): 1, ('C', 'I'): 1,
               ('C', 'K'): 0, ('C', 'R'): 0, ('C', 'Y'): 1, ('C', 'B'): 1, ('C', 'D'): 0, ('C', 'H'): 1, ('C', 'V'): 1, ('C', 'N'): 0, ('C', 'Z'): 0,
               ('G', 'A'): 0, ('G', 'C'): 0, ('G', 'G'): 1, ('G', 'T'): 0, ('G', 'W'): 0, ('G', 'S'): 1, ('G', 'M'): 0, ('G', 'I'): 1, 
               ('G', 'K'): 1, ('G', 'R'): 1, ('G', 'Y'): 0, ('G', 'B'): 1, ('G', 'D'): 1, ('G', 'H'): 0, ('G', 'V'): 1, ('G', 'N'): 0, ('G', 'Z'): 0,
               ('T', 'A'): 0, ('T', 'C'): 0, ('T', 'G'): 0, ('T', 'T'): 1, ('T', 'W'): 1, ('T', 'S'): 0, ('T', 'M'): 0, ('T', 'I'): 1, 
               ('T', 'K'): 1, ('T', 'R'): 0, ('T', 'Y'): 1, ('T', 'B'): 1, ('T', 'D'): 1, ('T', 'H'): 1, ('T', 'V'): 0, ('T', 'N'): 0, ('T', 'Z'): 0,
               ('Y', 'A'): 0, ('Y', 'C'): 1, ('Y', 'G'): 0, ('Y', 'T'): 1, ('Y', 'W'): 0, ('Y', 'S'): 0, ('Y', 'M'): 0, ('Y', 'I'): 0, 
               ('Y', 'K'): 0, ('Y', 'R'): 0, ('Y', 'Y'): 1, ('Y', 'B'): 0, ('Y', 'D'): 0, ('Y', 'H'): 0, ('Y', 'V'): 0, ('Y', 'N'): 0, ('Y', 'Z'): 0,
               ('N', 'A'): 0, ('N', 'C'): 0, ('N', 'G'): 0, ('N', 'T'): 0, ('N', 'W'): 0, ('N', 'S'): 0, ('N', 'M'): 0, ('N', 'I'): 0,
               ('N', 'K'): 0, ('N', 'R'): 0, ('N', 'Y'): 0, ('N', 'B'): 0, ('N', 'D'): 0, ('N', 'H'): 0, ('N', 'V'): 0, ('N', 'N'): 0, ('N', 'Z'): 0,
               ('W', 'A'): 1, ('W', 'B'): 0, ('W', 'C'): 0, ('W', 'D'): 0, ('W', 'G'): 0, ('W', 'H'): 0, ('W', 'I'): 0, ('W', 'K'): 0, ('W', 'M'): 0, 
               ('W', 'N'): 0, ('W', 'R'): 0, ('W', 'S'): 0, ('W', 'T'): 1, ('W', 'V'): 0, ('W', 'W'): 1, ('W', 'Y'): 0, ('W', 'Z'): 0,
               ('M', 'A'): 1, ('M', 'B'): 0, ('M', 'C'): 1, ('M', 'D'): 0, ('M', 'G'): 0, ('M', 'H'): 0, ('M', 'I'): 0, ('M', 'K'): 0, ('M', 'M'): 1, 
               ('M', 'N'): 0, ('M', 'R'): 0, ('M', 'S'): 0, ('M', 'T'): 0, ('M', 'V'): 0, ('M', 'W'): 0, ('M', 'Y'): 1, ('M', 'Z'): 0,
               ('R', 'A'): 1, ('R', 'B'): 0, ('R', 'C'): 0, ('R', 'D'): 0, ('R', 'G'): 1, ('R', 'H'): 0, ('R', 'I'): 0, ('R', 'K'): 0, ('R', 'M'): 0, 
               ('R', 'N'): 0, ('R', 'R'): 1, ('R', 'S'): 0, ('R', 'T'): 0, ('R', 'V'): 0, ('R', 'W'): 0, ('R', 'Y'): 0, ('R', 'Z'): 0,
               ('S', 'A'): 0, ('S', 'B'): 0, ('S', 'C'): 1, ('S', 'D'): 0, ('S', 'G'): 1, ('S', 'H'): 0, ('S', 'I'): 0, ('S', 'K'): 0, ('S', 'M'): 0, 
               ('S', 'N'): 0, ('S', 'R'): 0, ('S', 'S'): 1, ('S', 'T'): 0, ('S', 'V'): 0, ('S', 'W'): 0, ('S', 'Y'): 0, ('S', 'Z'): 0,
               ('K', 'A'): 0, ('K', 'B'): 0, ('K', 'C'): 0, ('K', 'D'): 0, ('K', 'G'): 1, ('K', 'H'): 0, ('K', 'I'): 0, ('K', 'K'): 1, ('K', 'M'): 0,
               ('K', 'N'): 0, ('K', 'R'): 0, ('K', 'S'): 0, ('K', 'T'): 1, ('K', 'V'): 0, ('K', 'W'): 0, ('K', 'Y'): 0, ('K', 'Z'): 0,
               ('B', 'A'): 0, ('B', 'B'): 1, ('B', 'C'): 1, ('B', 'D'): 0, ('B', 'G'): 1, ('B', 'H'): 0, ('B', 'I'): 0, ('B', 'K'): 0, ('B', 'M'): 0,
               ('B', 'N'): 0, ('B', 'R'): 0, ('B', 'S'): 0, ('B', 'T'): 1, ('B', 'V'): 0, ('B', 'W'): 0, ('B', 'Y'): 0, ('B', 'Z'): 0,
               ('D', 'A'): 1, ('D', 'B'): 0, ('D', 'C'): 0, ('D', 'D'): 1, ('D', 'G'): 1, ('D', 'H'): 0, ('D', 'I'): 0, ('D', 'K'): 0, ('D', 'M'): 0, 
               ('D', 'N'): 0, ('D', 'R'): 0, ('D', 'S'): 0, ('D', 'T'): 1, ('D', 'V'): 0, ('D', 'W'): 0, ('D', 'Y'): 0, ('D', 'Z'): 0,
               ('H', 'A'): 1, ('H', 'B'): 0, ('H', 'C'): 1, ('H', 'D'): 0, ('H', 'G'): 0, ('H', 'H'): 1, ('H', 'I'): 0, ('H', 'K'): 0, ('H', 'M'): 0, 
               ('H', 'N'): 0, ('H', 'R'): 0, ('H', 'S'): 0, ('H', 'T'): 1, ('H', 'V'): 0, ('H', 'W'): 0, ('H', 'Y'): 0, ('H', 'Z'): 0,
               ('V', 'A'): 1, ('V', 'B'): 0, ('V', 'C'): 1, ('V', 'D'): 0, ('V', 'G'): 1, ('V', 'H'): 0, ('V', 'I'): 0, ('V', 'K'): 0, ('V', 'M'): 0, 
               ('V', 'N'): 0, ('V', 'R'): 0, ('V', 'S'): 0, ('V', 'T'): 0, ('V', 'V'): 1, ('V', 'W'): 0, ('V', 'Y'): 0, ('V', 'Z'): 0,
               ('Z', 'A'): 0, ('Z', 'B'): 0, ('Z', 'C'): 0, ('Z', 'D'): 0, ('Z', 'G'): 0, ('Z', 'H'): 0, ('Z', 'I'): 0, ('Z', 'K'): 0, ('Z', 'M'): 0, 
               ('Z', 'N'): 0, ('Z', 'R'): 0, ('Z', 'S'): 0, ('Z', 'T'): 0, ('Z', 'V'): 0, ('Z', 'W'): 0, ('Z', 'Y'): 0, ('Z', 'Z'): 0,
               ('I', 'A'): 1, ('I', 'B'): 1, ('I', 'C'): 1, ('I', 'D'): 1, ('I', 'G'): 1, ('I', 'H'): 1, ('I', 'I'): 1, ('I', 'K'): 1, ('I', 'M'): 1, 
               ('I', 'N'): 1, ('I', 'R'): 1, ('I', 'S'): 1, ('I', 'T'): 1, ('I', 'V'): 1, ('I', 'W'): 1, ('I', 'Y'): 1, ('I', 'Z'): 1}
                
"""
char = 'A'
for key in MATCH_TABLE:
    c = key[0]
    if(c != char):
        print('')
        char = c
    else:
        print(' ', end='')
    print(str(key)+': '+str(MATCH_TABLE[key])+',', end='')
    char = c
"""

    
def _compute_matching_pairwise2(max_miss_f, max_miss_r, pp, gen):
    gen_seq = gen.seq
    
    search_limit = pp.rlen+pp.max_amplicon
    if(pp.flen+search_limit>len(gen_seq)): #If primer pair (plus amplicon) is larger than the genomic sequence, abort
        return None
   
    max_score = pp.rlen - max_miss_r + pp.flen - max_miss_f
    fpos = 0
    rpos = 0
    falignment = None
    ralignment = None
    
    forward_alignments = pairwise2.align.localds(gen_seq[0:-search_limit], pp.f.seq, MATCH_TABLE, -pp.flen, -pp.flen) #compute forward matches
    #forward_alignments = tuple(align1, align2(----sequence----), score, begin, end)
    for falign in forward_alignments: #for every forward match
        start = falign[4]+pp.min_amplicon #forward match start + len(forward) + min amplicon
        end = falign[4]+pp.max_amplicon+pp.rlen #f match start + len(f) + max amplicon + len(r)
        reverse_alignments = pairwise2.align.localds(gen_seq[start:end], pp.r.seq, MATCH_TABLE, -pp.rlen, -pp.rlen) #compute backward matches
        for ralign in reverse_alignments: #select best pair match
            if(falign[2]+ralign[2] > max_score):
                max_score = falign[2]+ralign[2]
                fpos = falign[4]-pp.flen #falign[3] #see important reminder
                rpos = ralign[4]-pp.rlen #ralign[3]
                falignment = falign
                ralignment = ralign
    
    if(falignment == None): #if there is no match return none
        return None
   
    amplicon = pp.min_amplicon+rpos
    rpos_abs = amplicon + fpos + pp.flen
    
    return Matching(gen, pp, fpos, rpos_abs, pp.flen-falignment[2], pp.rlen-ralignment[2], amplicon, MATCH_TABLE)
    

def compute_template_missmatches(fmaxm, rmaxm, primer_pairs, gen_record):
    result = []
    for gen_key in gen_record:
        print(gen_key)
        gen = gen_record[gen_key]
        gen_matches = MatchingList(gen)
        for pp in primer_pairs:
            match = _compute_matching_pairwise2(fmaxm, rmaxm, pp, gen)
            if(match): gen_matches.append(match)
            
        result.append(gen_matches)
    
    return result

def cutom_pw2_test(clean=True):
    
    """What should return
    AGCTAGCWTAGTTGCCA
    ---TAGCT---------	5.0
    -------TAGCT-----	2.0	
    --------TAGCT----	4.0
    -----------TAGCT-	3.0
    ------------TAGCT	2.0
    """
    gen = "AGCTAGCWTAGTTGCCA"
    primer = "TAGCT"

    print(pairwise2.align.localds(gen, primer, MATCH_TABLE, -20, -20))
    print("Custom")
    
    
    result = custom_pairwise2.align.localds(gen, primer, MATCH_TABLE, -20, -20)
    
    if(clean):
        clean_result = []
        starts = []
        for r in result:
            if r[3] not in starts:
                clean_result.append(r)
                starts.append(r[3])
        result = clean_result
            
    for r in result:
        print(r[0])
        print(r[1])
        print(r[2])
    

if (__name__=="__main__"):
    gen_record = ld.load_bio_file("Data/species_bold_own_genbank.fasta")
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")
    """
    result = compute_template_missmatches(10, 10, primer_pairs, gen_record)
    """
    """
    #Individual testing
    gen_name =  "FBARB265-11_Clubiona_leucaspis_BOLD" #"SPSLO262-12_Macaroeris_nidicolens_BOLD" "XXX-99_Xysticus_OWN"
    gen = gen_record.get(gen_name)
    gen = {gen_name:gen}
    primer_list = [primer_pairs[5-1]]
    result = compute_template_missmatches(10, 10, primer_list, gen)
    print(result[0])
    """
    cutom_pw2_test(clean=False)
        
    
        
   