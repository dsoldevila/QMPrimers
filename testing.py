#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 11:37:52 2019

@author: david

The purpose of this code is to check if this program is working properly
"""
import load_data as ld
from common import *
import matching as m
import interface as i

import cProfile
import pstats
import time
import csv
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import simulation as s

def split(gen_record, percent):
    
    partial_gen_record = {}
    max_gens=len(gen_record)*percent
    i = 0
    for gen in gen_record:
        partial_gen_record.update({gen: gen_record[gen]})
        i+=1
        if(i>max_gens):
            break
    return partial_gen_record
    

def test_all_pairs():
    #"primerPair","id","fastaid","organism","subgrup","primerF","primerR","mismFT","mismRT","mismTT","mismF3","mismR3","mismT3","long"
    trusted_results = pd.read_csv("Test_data/mismatches_allPrimers_allMitochondria.csv", sep=',')
    global_check = {"amplicon": 1, "missf": 1, "missr":1}
    check = {"amplicon": 1, "missf": 1, "missr":1}
    gen_record = ld.load_bio_files(["Data/mitochondrion.1.1.genomic.fna"], writable=True)
    gen_record = split(gen_record, 0.01)
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")
    
    result = m.compute_gen_matching(5, 5, primer_pairs, gen_record, 0, hanging_primers=True)
    gen_matching_table = result[0]
    
    header = TEMPLATE_HEADER
    correct_alignments = pd.DataFrame(columns=header)
    
    info = {"total_gens": len(gen_record), "matches_skipped":0, "alignments_processed": 0, "multiple_alignment_cases":0, "better_alignments":0}

    for index in gen_matching_table.index:
        genid = gen_matching_table.iloc[index].loc["fastaid"]
        pp = gen_matching_table.iloc[index].loc["primerPair"]
        target = trusted_results.loc[trusted_results['fastaid'].str.contains(genid[4:-1])]
        target = target.loc[target['primerPair'] == int(pp)]
        if(target.empty):
            #print("Target empty, skipping this gen...")
            info["matches_skipped"]+=1
            #Add them to not checked list
        else:
            pp = primer_pairs[pp]
            info["alignments_processed"]+=1
            #amplicon
            amplicon = gen_matching_table.iloc[index].loc["amplicon"]
            if(amplicon < pp.min_amplicon or amplicon > pp.max_amplicon):  
                print (genid)
                print("PRIMER PAIR: ", pp.id)
                print("Amplicon outside range")
                print(amplicon, pp.min_amplicon, pp.max_amplicon)
                global_check["amplicon"]=0
                check["amplicon"]=0
            if(amplicon != target['long'].iat[0]):
                print (genid)
                print("PRIMER PAIR: ", pp.id)
                print("Amplicon not matching")
                global_check["amplicon"]=0
                check["amplicon"]=0
                
            fm = target['mismFT'].iat[0]
            rm = target['mismRT'].iat[0]
            #fm
            if(gen_matching_table.iloc[index].loc["mismFT"] > fm):
                print (genid)
                print("PRIMER PAIR: ", pp.id)
                print("Bad forward matching")
                global_check["missf"]=0
                check["missf"]=0
            #rm
            if(gen_matching_table.iloc[index].loc["mismRT"] > rm):
                print (genid)
                print("PRIMER PAIR: ", pp.id)
                print("Bad reverse matching")
                global_check["missr"]=0
                check["missr"]=0
            if 0 not in check.values(): #if everything is correct check if the result found is better
                correct_alignments.loc[correct_alignments.shape[0]] = gen_matching_table.loc[index]
                if(gen_matching_table.iloc[index].loc["mismFT"]+gen_matching_table.iloc[index].loc["mismRT"] < fm+rm):
                    info["better_alignments"]+=1
                else:
                    info["better_alignments"]+=0
                    #correct_alignment_list.append(al)
                            
    for info_key in info:
        print(info_key, info[info_key])
        
    sum_check = 0
    for c in global_check:
        sum_check += check[c]
    if(sum_check == len(check)):
        print("SUCCESS!")
    else:
        print("TEST FAILED")
        
    correct_alignments.to_csv("Test_data/correct_alignments.csv", index_label="id")
    """
    store_results("Test_data/better_alignments.csv", better_alignment_list, header)
    store_results("Test_data/not_tested_alignments.csv", not_tested_alignment_list, header)
    ld.store_matching_results("Test_data/full_alignments.csv", gen_alignment_list, header)
    """
    
    return


def pandas_scalability_test():
    """
    @brief: check if time is linearly proportional to matrix size
    """
    gen_record_large = ld.load_bio_files(["Data/mitochondrion.1.1.genomic.fna"])
    gen_record = ld.load_bio_files(["Data/species_bold_own_genbank.fasta"])
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")
    
    max_len = 0
    key = None
    for gen_key in gen_record_large:
        leng = len(gen_record_large[gen_key])
        if (leng>max_len):
            max_len = leng
            key = gen_key
            
    primer = primer_pairs[4].f
    #largest_matrix len= 2M*primer_len
    gen = gen_record_large[key]
    time1 = time.time()
    result_matrix = m.MATCH_TABLE.loc[gen, primer]
    elapsedTime_l = ((time.time()-time1))
    
    #len = 658*primer_len
    gen = gen_record["ACEA563-14_Aphis_gossypii_BOLD"]
    len_s = len(gen)
    time1 = time.time()
    result_matrix = m.MATCH_TABLE.loc[gen, primer]
    elapsedTime_s = ((time.time()-time1))
    
    print(int(elapsedTime_l*(10**9)/max_len), int(elapsedTime_s*(10**9)/len_s))
    
    return


def check_if_multiple_alignments_are_frequent():
    gen_record = ld.load_bio_files(["Data/species_bold_own_genbank.fasta"]) 
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")

    gen_alignment_list = m.compute_gen_matching(5, 5, primer_pairs, gen_record) 
    for gm in gen_alignment_list:
       matching_list = gm.get_matching_list()
       for al_list in matching_list:
           alignments = al_list.get_list()
           if(len(alignments>1)):
               print(al_list.gen.id, al_list.primer_pair.id, len(alignments))
    return

def performance_test():
    cProfile.run('restore_template()', 'temp.profile')
    stats = pstats.Stats('temp.profile')
    stats.strip_dirs().sort_stats('cumtime').print_stats(10)

        
def check_uppercase():
    gen_record = ld.load_bio_files(["Data/sbog_test.fasta"], check_uppercase=True) 
    print(gen_record["ACEA1016-14_Aphis_spiraecola_BOLD"])
    
def restore_template():
    gen_record = ld.load_bio_files(["Data/species_bold_own_genbank.fasta"])
    #gen_record = {"AGB001-11_Salticus_scenicus_BOLD": gen_record["AGB001-11_Salticus_scenicus_BOLD"]};
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")
    #primer_pairs = [primer_pairs[5]]
    template, discarded, rs, cs = m.compute_gen_matching(5, 5, primer_pairs, gen_record, 0) 
    
    header = ["primerPair","fastaid","primerF","primerR","mismFT","mismRT","amplicon", "F_pos", "mismFT_loc", "mismFT_type", 
                                     "mismFT_base", "R_pos"]
    m.store_matching_results("Test_data/test1.csv", template, header=TEMPLATE_HEADER)
    templateR = ld.load_template("Test_data/test1.csv")
    templateR, discarded, rs, cs = ld.restore_template(templateR, gen_record, primer_pairs, 10)
    
    i.save_matching_info("Test_data/test2", templateR, discarded, rs, cs, header=TEMPLATE_HEADER)
    return

def simulate_whitebox():
    template = pd.read_csv("/home/david/Git/QMPrimers/Test_data/test1.csv")
    full_sample = template["fastaid"].unique()
    sample_size = 6
    primer_pair = 14
    k = 0.3
    B = 4
    sim = s.Simulation(template, sample_size)
    sample = sim.get_random_sample(full_sample, sample_size, k)
    
    sim.amplify(template, primer_pair, sample, B)
    print(sample)
    tmp = sample[["oprop", "fprop"]]
    tmp = tmp.astype('float64')
    tmp = tmp.corr()
    print(tmp)

def simulate():
    template = None
    sample_size = 10
    k = 0.3
    B = 4
    N=10
    ci = 0.90
    sim = s.Simulation()
    raw, cooked = sim.simulate(template, sample_size, k, B, N, ci)
    sim.store_data("sim_test", raw, cooked, "output_positive.csv", sample_size, k, B, N)
    
def matching_test():
    init_logger()
    set_verbosity(True)
    gen_record = ld.load_bio_files(["Data/plant_dna.fasta"])
    #gen_record = split(gen_record, 0.005)
    #gen_record = {"AGB001-11_Salticus_scenicus_BOLD": gen_record["AGB001-11_Salticus_scenicus_BOLD"]};
    primer_pairs = ld.load_csv_file("Data/PP_chl_one.csv")
    #primer_pairs = {"6":primer_pairs["6"]}
    
    output = "Test_data/output"
    header = [i for i in range(len(TEMPLATE_HEADER))]
    template, discarded, rs, cs = m.compute_gen_matching(5, 5, primer_pairs, gen_record, output) 
    
    #m.store_matching_results("Test_data/test1.csv", template, header=TEMPLATE_HEADER)
    i.save_matching_info("TEST", output, template, header, discarded, rs, cs)
    close_logger()
    print(gen_record["NC_018766.1"][57613:57613+30])
    print(''.join(primer_pairs["1"].r.seq))
    return

def get_nend_test():
    nend = 2
    template = pd.read_csv("/home/david/Git/QMPrimers/Test_data/test_positive.csv")
    template2, raw, cooked = m.get_Nend_template(template, nend)
    
    header = ["primerPair","fastaid","primerF","primerR","mismFN"+str(nend),"mismRN"+str(nend),"amplicon", "F_pos", 
              "mismFN"+str(nend)+"_loc", "mismFN"+str(nend)+"_type", "mismFN"+str(nend)+"_base", "R_pos", "mismRN"+str(nend)+"_loc",
              "mismRN"+str(nend)+"_type", "mismRN"+str(nend)+"_base"]
    
    m.store_matching_results("Test_data/test_nend.csv", template2, header=header)
    return

parameters = [
        ["gen", "Data/sbog_test.fasta", "Genome file dir, no support for multiple files in cl", "-gf", "entry"],
        ["primer_pairs", "Data/PP.csv", "Primer pairs file dir. A particular header must be used in the file", "-pf", "entry"],
        ["output_file", os.path.join(os.getcwd(),"test"), "Location of the output files, no extension", "-o", "entry"],
        ["forward missmatches", 5, "Maximum number of missmatches allowed on forward primer", "-fm", "param"],
        ["reverse missmatches", 5, "Maximum number of missmatches allowed on reverse primer", "-rm", "param"], 
        ["Nend miss.", 0, "Missmatches in the last N positions on forward and in the first N pos. on reverse ", "-nend", "info"],
        ["hanging primers", False, "Primers allowed to match between [0-mf,len(genome)+mr] instead of just between genome's length", "--hanging", "param"],
        ["check_integrity", False, "Checks integrity of gen files, integrity of primer file is always checked", "--checki", "param"],
        ["check_uppercase", False, "Checks that all gens are in upper case, lower case gens will trigger an integrity file", "--checku", "param"],
        ["csv_template", "output_positive.csv", "Precomputed missmatching template", "-i", "entry"],
        ["verbose", False, "Outputs extra information", "--v", "cmd"]]                      
parameters = pd.DataFrame([x[1:] for x in parameters], index = [x[0] for x in parameters], columns=["value", "description", "flag", "type"])

def load_template_test():
    template, discarded, gen_record, primer_pairs, raw_stats, cooked_stats = i.load_template(parameters)
    print(template)
    
    nend=5
    nend_template, raw_stats, cooked_stats = i.get_Nend_match(template, nend, 10)
    
    header = [x for x in range(len(TEMPLATE_HEADER))]
    m.store_matching_results("test", "test.csv", template, header=header)
    
    
    

if(__name__=="__main__"):
    init_logger()
    time1 = time.time()
    simulate()
    elapsedTime = ((time.time()-time1))
    print(int(elapsedTime)/60)
    
    #performance_test()
    close_logger()
