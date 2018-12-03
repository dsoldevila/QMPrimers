#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 15:34:21 2018

@author: david
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from common import *

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Alphabet import IUPAC

import load_data as ld


#A,C,G,T,Y,N,-,W,M,R,
MATCH_TABLE = {('A', 'A'): 1, ('A', 'C'): 0, ('A', 'G'): 0, ('A', 'T'): 0, ('A', 'U'): 0, ('A', 'W'): 1, ('A', 'S'): 0, ('A', 'M'): 1, ('A', 'I'): 1,
               ('A', 'K'): 0, ('A', 'R'): 1, ('A', 'Y'): 0, ('A', 'B'): 0, ('A', 'D'): 1, ('A', 'H'): 1, ('A', 'V'): 1, ('A', 'N'): 0, ('A', 'Z'): 0,
               ('C', 'A'): 0, ('C', 'C'): 1, ('C', 'G'): 0, ('C', 'T'): 0, ('C', 'U'): 0, ('C', 'W'): 0, ('C', 'S'): 1, ('C', 'M'): 1, ('C', 'I'): 1, 
               ('C', 'K'): 0, ('C', 'R'): 0, ('C', 'Y'): 1, ('C', 'B'): 1, ('C', 'D'): 0, ('C', 'H'): 1, ('C', 'V'): 1, ('C', 'N'): 0, ('C', 'Z'): 0,
               ('G', 'A'): 0, ('G', 'C'): 0, ('G', 'G'): 1, ('G', 'T'): 0, ('G', 'U'): 0, ('G', 'W'): 0, ('G', 'S'): 1, ('G', 'M'): 0, ('G', 'I'): 1,
               ('G', 'K'): 1, ('G', 'R'): 1, ('G', 'Y'): 0, ('G', 'B'): 1, ('G', 'D'): 1, ('G', 'H'): 0, ('G', 'V'): 1, ('G', 'N'): 0, ('G', 'Z'): 0,
               ('T', 'A'): 0, ('T', 'C'): 0, ('T', 'G'): 0, ('T', 'T'): 1, ('T', 'U'): 0, ('T', 'W'): 1, ('T', 'S'): 0, ('T', 'M'): 0, ('T', 'I'): 1, 
               ('T', 'K'): 1, ('T', 'R'): 0, ('T', 'Y'): 1, ('T', 'B'): 1, ('T', 'D'): 1, ('T', 'H'): 1, ('T', 'V'): 0, ('T', 'N'): 0, ('T', 'Z'): 0,
               ('Y', 'A'): 0, ('Y', 'C'): 1, ('Y', 'G'): 0, ('Y', 'T'): 1, ('Y', 'U'): 0, ('Y', 'W'): 1, ('Y', 'S'): 1, ('Y', 'M'): 1, ('Y', 'I'): 1, 
               ('Y', 'K'): 1, ('Y', 'R'): 0, ('Y', 'Y'): 1, ('Y', 'B'): 1, ('Y', 'D'): 1, ('Y', 'H'): 1, ('Y', 'V'): 1, ('Y', 'N'): 0, ('Y', 'Z'): 0,
               ('N', 'A'): 0, ('N', 'C'): 0, ('N', 'G'): 0, ('N', 'T'): 0, ('N', 'U'): 0, ('N', 'W'): 0, ('N', 'S'): 0, ('N', 'M'): 0, ('N', 'I'): 0, 
               ('N', 'K'): 0, ('N', 'R'): 0, ('N', 'Y'): 0, ('N', 'B'): 0, ('N', 'D'): 0, ('N', 'H'): 0, ('N', 'V'): 0, ('N', 'N'): 0, ('N', 'Z'): 0,
               ('W', 'A'): 1, ('W', 'B'): 1, ('W', 'C'): 0, ('W', 'D'): 1, ('W', 'G'): 0, ('W', 'H'): 1, ('W', 'I'): 1, ('W', 'K'): 1, ('W', 'M'): 1, 
               ('W', 'N'): 0, ('W', 'R'): 1, ('W', 'S'): 0, ('W', 'T'): 1, ('W', 'U'): 0, ('W', 'V'): 1, ('W', 'W'): 1, ('W', 'Y'): 1, ('W', 'Z'): 0,
               ('M', 'A'): 1, ('M', 'B'): 0, ('M', 'C'): 1, ('M', 'D'): 0, ('M', 'G'): 0, ('M', 'H'): 0, ('M', 'I'): 0, ('M', 'K'): 0, ('M', 'M'): 1, 
               ('M', 'N'): 0, ('M', 'R'): 0, ('M', 'S'): 0, ('M', 'T'): 0, ('M', 'U'): 0, ('M', 'V'): 0, ('M', 'W'): 0, ('M', 'Y'): 1, ('M', 'Z'): 0,
               ('R', 'A'): 0, ('R', 'B'): 0, ('R', 'C'): 0, ('R', 'D'): 0, ('R', 'G'): 0, ('R', 'H'): 0, ('R', 'I'): 0, ('R', 'K'): 0, ('R', 'M'): 0, 
               ('R', 'N'): 0, ('R', 'R'): 1, ('R', 'S'): 0, ('R', 'T'): 0, ('R', 'U'): 0, ('R', 'V'): 0, ('R', 'W'): 0, ('R', 'Y'): 0, ('R', 'Z'): 0,
               ('-', 'A'): 0, ('-', 'C'): 0, ('-', 'G'): 0, ('-', 'T'): 0, ('-', 'U'): 0, ('-', 'W'): 0, ('-', 'S'): 0, ('-', 'M'): 0, ('-', 'I'): 0, 
               ('-', 'K'): 0, ('-', 'R'): 0, ('-', 'Y'): 0, ('-', 'B'): 0, ('-', 'D'): 0, ('-', 'H'): 0, ('-', 'V'): 0, ('-', 'N'): 0, ('-', 'Z'): 0}
                
"""
for key in MATCH_TABLE:
    c = key[0]
    if(c != char):
        print('')
        char = c
    else:
        print(' ', end='')
    print(str(key)+': '+str(MATCH_TABLE[key])+',', end='')"""
    
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
    
    forward_alignments = pairwise2.align.localds(gen_seq[0:-search_limit], pp.f.seq, MATCH_TABLE, -5, -5) #compute forward matches
    for falign in forward_alignments: #for every forward match
        start = falign[3]+pp.flen+pp.min_amplicon #forward match start + len(forward) + min amplicon
        end = falign[3]+pp.flen+pp.max_amplicon+pp.rlen #f match start + len(f) + max amplicon + len(r)
        reverse_alignments = pairwise2.align.localds(gen_seq[start:end], pp.r.seq, MATCH_TABLE, -5, -5) #compute backward matches
        for ralign in reverse_alignments: #select best pair match
            if(falign[2]+ralign[2] > max_score):
                max_score = falign[2]+ralign[2]
                fpos = falign[3]
                rpos = ralign[3]
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

if (__name__=="__main__"):
    gen_record = ld.load_bio_file("Data/species_bold_own_genbank.fasta")
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")
    result = compute_template_missmatches(10, 10, primer_pairs, gen_record)
    print(result[0])
        
    """gen = gen_record.get("GBCH0808-06_Drassodes_lapidosus_BOLD")
    pp = primer_pairs[0]
    match = _compute_matching_pairwise2(11, 11, pp, gen)
    print(match)"""

    """
    alignments = pairwise2.align.localms(gen_record.get("ACEA1016-14_Aphis_spiraecola_BOLD").seq, primer_pairs[1].f.seq, 1, 0, -5, -5)
    for a in alignments:
        #a = tuple(align1, align2(----sequence----), score, begin, end)
        #print(format_alignment(*a))
        print(a[4])"""