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


#IUPACAmbiguousDNA
MATCH_TABLE = {('A', 'A'): 1, ('A', 'C'): 0, ('A', 'G'): 0, ('A', 'T'): 0, ('A', 'U'): 0, ('A', 'W'): 1, ('A', 'S'): 0, ('A', 'M'): 1, ('A', 'I'): 1,
               ('A', 'K'): 0, ('A', 'R'): 1, ('A', 'Y'): 0, ('A', 'B'): 0, ('A', 'D'): 1, ('A', 'H'): 1, ('A', 'V'): 1, ('A', 'N'): 1, ('A', 'Z'): 0,
               ('C', 'A'): 0, ('C', 'C'): 1, ('C', 'G'): 0, ('C', 'T'): 0, ('C', 'U'): 0, ('C', 'W'): 0, ('C', 'S'): 1, ('C', 'M'): 1, ('C', 'I'): 1,
               ('C', 'K'): 0, ('C', 'R'): 0, ('C', 'Y'): 1, ('C', 'B'): 1, ('C', 'D'): 0, ('C', 'H'): 1, ('C', 'V'): 1, ('C', 'N'): 1, ('C', 'Z'): 0, 
               ('G', 'A'): 0, ('G', 'C'): 0, ('G', 'G'): 1, ('G', 'T'): 0, ('G', 'U'): 0, ('G', 'W'): 0, ('G', 'S'): 1, ('G', 'M'): 0, ('G', 'I'): 1,
               ('G', 'K'): 1, ('G', 'R'): 1, ('G', 'Y'): 0, ('G', 'B'): 1, ('G', 'D'): 1, ('G', 'H'): 0, ('G', 'V'): 1, ('G', 'N'): 1, ('G', 'Z'): 0,
               ('T', 'A'): 0, ('T', 'C'): 0, ('T', 'G'): 0, ('T', 'T'): 1, ('T', 'U'): 0, ('T', 'W'): 1, ('T', 'S'): 0, ('T', 'M'): 0, ('T', 'I'): 1,
               ('T', 'K'): 1, ('T', 'R'): 0, ('T', 'Y'): 1, ('T', 'B'): 1, ('T', 'D'): 1, ('T', 'H'): 1, ('T', 'V'): 0, ('T', 'N'): 1, ('T', 'Z'): 0,
               ('Y', 'A'): 0, ('Y', 'C'): 1, ('Y', 'G'): 0, ('Y', 'T'): 1, ('Y', 'U'): 0, ('Y', 'W'): 1, ('Y', 'S'): 1, ('Y', 'M'): 1, ('Y', 'I'): 1,
               ('Y', 'K'): 1, ('Y', 'R'): 0, ('Y', 'Y'): 1, ('Y', 'B'): 1, ('Y', 'D'): 1, ('Y', 'H'): 1, ('Y', 'V'): 1, ('Y', 'N'): 1, ('Y', 'Z'): 0,
               ('N', 'A'): -5, ('N', 'C'): -5, ('N', 'G'): -5, ('N', 'T'): -5, ('N', 'U'): -5, ('N', 'W'): -5, ('N', 'S'): -5, ('N', 'M'): -5, ('N', 'I'): -5,
               ('N', 'K'): -5, ('N', 'R'): -5, ('N', 'Y'): -5, ('N', 'B'): -5, ('N', 'D'): -5, ('N', 'H'): -5, ('N', 'V'): -5, ('N', 'N'): -5, ('N', 'Z'): -5}


    
def _compute_matching_pairwise2(max_miss_f, max_miss_r, pp, gen):
    gen_seq = gen.seq
    
    max_score = 0
    fpos = 0
    rpos = 0
    falignment = None
    ralignment = None
    
    search_limit = pp.rlen+pp.max_amplicon
    if(pp.flen+search_limit>len(gen_seq)): #If primer pair (plus amplicon) is larger than the genomic sequence, abort
        return None
    
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
    for gen_matching in result:
        print(result)
        
    """gen = gen_record.get("AGB001-11_Salticus_scenicus_BOLD")
    pp = primer_pairs[14]
    match = _compute_matching_pairwise2(11, 11, pp, gen)
    print(match)"""

    """
    alignments = pairwise2.align.localms(gen_record.get("ACEA1016-14_Aphis_spiraecola_BOLD").seq, primer_pairs[1].f.seq, 1, 0, -5, -5)
    for a in alignments:
        #a = tuple(align1, align2(----sequence----), score, begin, end)
        #print(format_alignment(*a))
        print(a[4])"""