#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#DEPRECATED
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
MATCH_TABLE = {"A":{"A":1, "C":0, "G":0, "T":0, "U":0, "W":1, "S":0,"M":1,"K":0,"R":1,"Y":0,"B":0,"D":1,"H":1,"V":1,"N":1,"Z":0}, 
               "C":{"A":0, "C":1, "G":0, "T":0, "U":0, "W":0, "S":1,"M":1,"K":0,"R":0,"Y":1,"B":1,"D":0,"H":1,"V":1,"N":1,"Z":0}, 
               "G":{"A":0, "C":0, "G":1, "T":0, "U":0, "W":0, "S":1,"M":0,"K":1,"R":1,"Y":0,"B":1,"D":1,"H":0,"V":1,"N":1,"Z":0},
               "T":{"A":0, "C":0, "G":0, "T":1, "U":0, "W":1, "S":0,"M":0,"K":1,"R":0,"Y":1,"B":1,"D":1,"H":1,"V":0,"N":1,"Z":0}}
            

class Matching_old:
    def __init__(self, gen, primer_pair, range_f, range_r, misses_f, misses_f_loc, misses_r, misses_r_loc, amplicon):
        self.gen = gen
        self.primer_pair = primer_pair
        self.range_f = range_f
        self.range_r = range_r
        self.misses_f = misses_f
        self.misses_f_loc = self._misses2str(misses_f_loc)
        self.misses_r = misses_r
        self.misses_r_loc = self._misses2str(misses_r_loc)
        self.amplicon = amplicon
        
    def _misses2str(self, misses_loc):
        array = " " if misses_loc[0] else "|"
        for i in range(1, len(misses_loc)):
            char = " " if misses_loc[i] else "|"
            array +=char
        return array
        
    def __str__(self):
        info = ("PRIME PAIR "+str(self.primer_pair.id)+"\n"+
              "Forward's missmatches: "+str(self.misses_f)+"\n"+
              "Backward's missmatches: "+str(self.misses_r)+"\n"+
              "Pair's amplicon: "+str(self.amplicon)+"\n"+
              "Forward's match "+ str(self.range_f)+"\n"+
              " "+str(self.gen.seq[self.range_f[0]:self.range_f[1]])+"\n"+
              " "+self.misses_f_loc+"\n"+
              " "+str(self.primer_pair.f.seq) +"\n"+
              "Backwards's match "+ str(self.range_r)+"\n"+
              " "+str(self.gen.seq[self.range_r[0]:self.range_r[1]])+"\n"+
              " "+self.misses_r_loc+"\n"+
              " "+str(self.primer_pair.r.seq) +"\n")
        return info
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
    
    pp = primer_pairs[1]
    gen_seq = gen_record.get("ACEA1016-14_Aphis_spiraecola_BOLD")
    seq = gen_seq.seq
    

    
    min_mf = pp.flen
    min_mr = pp.rlen
    i_f = 0
    i_r = 0
    range_max = len(seq)-pp.flen-pp.rlen-pp.max_amplicon-1
    
    missf_array = [0]*pp.flen
    missr_array = [0]*pp.rlen
    missf_loc = [0]*pp.flen
    missr_loc = [0]*pp.rlen
    
    for i in range(0, range_max): #max and min amplicon are equal
        missf = 0
        missr = 0
        for j in range(0, pp.flen): #primers don't have the same lenght
            r = not(MATCH_TABLE[seq[i+j]][pp.f.seq[j]])
            missf_array[j] = 1 if r else 0
            missf += r
        for j in range(0, pp.rlen):
            r = not(MATCH_TABLE[seq[i+j+pp.max_amplicon+pp.flen]][pp.r.seq[j]])
            missr_array[j] = 1 if r else 0
            missr += r
            
        con = missf <= min_mf and missr <= min_mr
        
        min_mf = missf if con else min_mf
        min_mr = missr if con else min_mr
        i_f = i if con else i_f
        i_r = i+pp.max_amplicon+pp.flen if con else i_r
        missf_loc = missf_array.copy() if con else missf_loc
        missr_loc = missr_array.copy() if con else missr_loc
    
    result = Matching_old(gen_seq, pp, (i_f, i_f+pp.flen), (i_r, i_r+pp.rlen), min_mf, missf_loc, min_mr, missr_loc, pp.max_amplicon)
    print(result)
    return

if (__name__=="__main__"):
    gen_record = ld.load_bio_file("Data/species_bold_own_genbank.fasta")
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")
    compute_matching(10, 10, primer_pairs, gen_record)

    """
    alignments = pairwise2.align.localms(gen_record.get("ACEA1016-14_Aphis_spiraecola_BOLD").seq, primer_pairs[1].f.seq, 1, 0, -5, -5)
    for a in alignments:
        #a = tuple(align1, align2(----sequence----), score, begin, end)
        #print(format_alignment(*a))
        print(a[4])"""