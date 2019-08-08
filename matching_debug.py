#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 12:04:10 2019

@author: david
"""
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import load_data as ld
import matching as m
from common import *



def debug_matching(gen, primer_pair, mf, mr, output_file, hanging_primers=False):
    assert(len(gen)==1)
    assert(len(primer_pair)==1)
    
    template, discarded, raw_stats, cooked_stats = m.compute_gen_matching(mf, mr, primer_pair, gen, output_file, hanging_primers=hanging_primers)
    
    if(template.empty):
        logging.warning("No result")
        return
    
    match_result = template.loc[0]
    pp = primer_pair[next(iter(primer_pair))]
    gen = gen[next(iter(gen))]
    
    rem_len = len(gen)-(match_result.at['F_pos']+pp.flen+match_result.at['ampliconLen']+pp.rlen)
    
    pp.f.seq = Seq(''.join(pp.f.seq))
    pp.r.seq = Seq(''.join(pp.r.seq))
    
    pp_aligned = '-'*match_result.at['F_pos']+pp.f.seq+'-'*match_result.at['ampliconLen']+pp.r.seq+'-'*rem_len
    pp_aligned=SeqRecord(pp_aligned)
    pp_aligned.id = pp.id
    align = MultipleSeqAlignment([gen, pp_aligned])
    print(align.format("clustal"))
    
if(__name__=="__main__"):
    
    
    ppf = ld.load_csv_file("test_input/single_PP.csv")
    genr = ld.load_bio_files("test_input/test_single.fasta")
    """
    gen = genr["ACEA1016-14_Aphis_spiraecola_BOLD"]
    primer = ppf['1'].f
    leng = len(gen)
    lenp = len(primer)
    lenr = leng-lenp
    primer = primer+"-"*lenr
    align = MultipleSeqAlignment([gen, primer])
    print(align.format("clustal"))
    """
    debug_matching(genr, ppf, 10, 10, "./debug")



