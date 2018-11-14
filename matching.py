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

import load_data as ld


#IUPACAmbiguousDNA
MATCH_TABLE = {"A":{"A":1, "C":0, "G":0, "T":0, "U":0, "W":1, "S":0,"M":1,"K":0,"R":1,"Y":0,"B":0,"D":1,"H":1,"V":1,"N":1,"Z":0}, 
               "C":{"A":0, "C":1, "G":0, "T":0, "U":0, "W":0, "S":1,"M":1,"K":0,"R":0,"Y":1,"B":1,"D":0,"H":1,"V":1,"N":1,"Z":0}, 
               "G":{"A":0, "C":0, "G":1, "T":0, "U":0, "W":0, "S":1,"M":0,"K":1,"R":1,"Y":0,"B":1,"D":1,"H":0,"V":1,"N":1,"Z":0},
               "T":{"A":0, "C":0, "G":0, "T":1, "U":0, "W":1, "S":0,"M":0,"K":1,"R":0,"Y":1,"B":1,"D":1,"H":1,"V":0,"N":1,"Z":0}}
        
"""def _pcr_amp(max_miss_f, max_miss_r, pp, gen_seq, model):
    
    len_f = len(pp.f.seq) #es té en compte que els dos primers tenen la mateixa llargada
    len_r = len(pp.r.seq)
    amplicon = pp.min_amplicon
    lim = len_f + len_r + amplicon
    
    min_miss_f = max_miss_f
    min_miss_r = max_miss_r
    
    i_f = 0
    i_r = 0
    
    
    for i in range(len(gen_seq.seq)-lim-1):
        miss_f = 0
        miss_r = 0
        j = 0
        for j in range(len_f-1):
            miss_f += MATCH_TABLE[gen_seq.seq[i+j]][pp.f.seq[j]]
            miss_r += MATCH_TABLE[gen_seq.seq[i+lim-len_r+j]][pp.f.seq[j]]
        
        con = miss_f<min_miss_f and miss_r<min_miss_r
        
        min_miss_f = miss_f if con else min_miss_f
        min_miss_r = miss_r if con else min_miss_r
        i_f = i if con else i_f
        i_r = i+lim-len_r if con else i_r    
        
    return Matching(gen_seq, pp, min_miss_f, min_miss_r, i_f, i_r, amplicon)
"""        

def compute_matching(max_miss_f, max_miss_r, primer_pairs, gen_record, model=1):
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
    """for gen in gen_seqs:
        for primer_pair in primer_pairs:"""
    
    pp = primer_pairs[0]
    gen_seq = gen_record.get("ACEA1016-14_Aphis_spiraecola_BOLD")
    alignments = pairwise2.align.globalxs(gen_seq.seq[0:int(len(gen_seq.seq)/8)], pp.f.seq, -1000, -1000)
    print('Score',alignments[0][2])
    print(format_alignment(*alignments[0]))
    
    pass

if "__main__":
    gen_record = ld.load_bio_file("species_bold_own_genbank.fasta")
    primer_pairs = ld.load_csv_file("P&PP.csv")
    
    compute_matching(5, 5, primer_pairs, gen_record)