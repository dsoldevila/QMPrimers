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
MATCH_TABLE_old = {"A":{"A":1, "C":0, "G":0, "T":0, "U":0, "W":1, "S":0,"M":1,"K":0,"R":1,"Y":0,"B":0,"D":1,"H":1,"V":1,"N":1,"Z":0}, 
               "C":{"A":0, "C":1, "G":0, "T":0, "U":0, "W":0, "S":1,"M":1,"K":0,"R":0,"Y":1,"B":1,"D":0,"H":1,"V":1,"N":1,"Z":0}, 
               "G":{"A":0, "C":0, "G":1, "T":0, "U":0, "W":0, "S":1,"M":0,"K":1,"R":1,"Y":0,"B":1,"D":1,"H":0,"V":1,"N":1,"Z":0},
               "T":{"A":0, "C":0, "G":0, "T":1, "U":0, "W":1, "S":0,"M":0,"K":1,"R":0,"Y":1,"B":1,"D":1,"H":1,"V":0,"N":1,"Z":0}}

MATCH_TABLE = {('A', 'A'): 1, ('A', 'C'): 0, ('A', 'G'): 0, ('A', 'T'): 0, ('A', 'U'): 0, ('A', 'W'): 1, ('A', 'S'): 0, ('A', 'M'): 1, 
               ('A', 'K'): 0, ('A', 'R'): 1, ('A', 'Y'): 0, ('A', 'B'): 0, ('A', 'D'): 1, ('A', 'H'): 1, ('A', 'V'): 1, ('A', 'N'): 1, ('A', 'Z'): 0,
               ('C', 'A'): 0, ('C', 'C'): 1, ('C', 'G'): 0, ('C', 'T'): 0, ('C', 'U'): 0, ('C', 'W'): 0, ('C', 'S'): 1, ('C', 'M'): 1, 
               ('C', 'K'): 0, ('C', 'R'): 0, ('C', 'Y'): 1, ('C', 'B'): 1, ('C', 'D'): 0, ('C', 'H'): 1, ('C', 'V'): 1, ('C', 'N'): 1, ('C', 'Z'): 0, 
               ('G', 'A'): 0, ('G', 'C'): 0, ('G', 'G'): 1, ('G', 'T'): 0, ('G', 'U'): 0, ('G', 'W'): 0, ('G', 'S'): 1, ('G', 'M'): 0,
               ('G', 'K'): 1, ('G', 'R'): 1, ('G', 'Y'): 0, ('G', 'B'): 1, ('G', 'D'): 1, ('G', 'H'): 0, ('G', 'V'): 1, ('G', 'N'): 1, ('G', 'Z'): 0,
               ('T', 'A'): 0, ('T', 'C'): 0, ('T', 'G'): 0, ('T', 'T'): 1, ('T', 'U'): 0, ('T', 'W'): 1, ('T', 'S'): 0, ('T', 'M'): 0,
               ('T', 'K'): 1, ('T', 'R'): 0, ('T', 'Y'): 1, ('T', 'B'): 1, ('T', 'D'): 1, ('T', 'H'): 1, ('T', 'V'): 0, ('T', 'N'): 1, ('T', 'Z'): 0}

    
def _compute_matching_pairwise2(max_miss_f, max_miss_r, pp, gen):
    
    #for every genome
    #for every pair
    #compute forward matches
    #compute backward matches
    #select best pair match
    gen = gen.seq
    
    max_score = 0
    fpos = 0
    rpos = 0
    falignment = None
    ralignment = None
    
    forward_alignments = pairwise2.align.localds(gen, pp.f.seq, MATCH_TABLE, -5, -5)
    for falign in forward_alignments:
        reverse_alignments = pairwise2.align.localds(gen[falign[3]+pp.flen+pp.min_amplicon:falign[3]+pp.flen+pp.max_amplicon+pp.rlen], pp.r.seq, MATCH_TABLE, -5, -5)
        for ralign in reverse_alignments:
            if(falign[2]+ralign[2] > max_score):
                max_score = falign[2]+ralign[2]
                fpos = falign[3]
                rpos = ralign[3]
                falignment = falign
                ralignment = ralign
    
            #si forward + backward score < mínim, guarda el resultat
        #a = tuple(align1, align2(----sequence----), score, begin, end)
        #print(format_alignment(*a))
    print("Forward's missmatches: %d at %d" % (pp.flen-falignment[2], fpos))
    print("Reverse's missmatches: %d at %d" % (pp.rlen-ralignment[2], rpos))
    print(format_alignment(*falignment))
    print(format_alignment(*ralignment))
    
    amplicon = pp.min_amplicon+rpos
    rpos_abs = amplicon + fpos + pp.flen
    result = Matching(gen_record.get("ACEA1016-14_Aphis_spiraecola_BOLD"), pp, fpos, rpos_abs, pp.flen-falignment[2], pp.rlen-ralignment[2], amplicon)
    print(result)
    
    return

def compute_template_missmatches(fmaxm, rmaxm, primer_pairs, gen_record):
    result = []
    for gen in gen_record:
        gen_matches = MatchingList(gen)
        for pp in primer_pairs:
            gen_matches.append(_compute_matching_pairwise2(fmaxm, rmaxm, pp, gen))
            
        result.append(gen_matches)
    
    return result

if (__name__=="__main__"):
    gen_record = ld.load_bio_file("Data/species_bold_own_genbank.fasta")
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")
    _compute_matching_pairwise2(10, 10, primer_pairs, gen_record)

    """
    alignments = pairwise2.align.localms(gen_record.get("ACEA1016-14_Aphis_spiraecola_BOLD").seq, primer_pairs[1].f.seq, 1, 0, -5, -5)
    for a in alignments:
        #a = tuple(align1, align2(----sequence----), score, begin, end)
        #print(format_alignment(*a))
        print(a[4])"""