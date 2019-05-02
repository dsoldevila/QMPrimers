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

"Name, forward starting pos, forward's missmatches loc, reverse's missmatches loc"
target1 = [("XXX-99_Xysticus_OWN", 11, [0,12,13,15],[2,5]),
        ("FBARB265-11_Clubiona_leucaspis_BOLD", 0,[0,12],[2,5,14,18]),
        ("ACEA902-14_Myzus_persicae_BOLD", 0,[13,15],[2,17,21]),
        ("ACEA563-14_Aphis_gossypii_BOLD", 0,[13,15],[5,21]),
        ("ACEA640-14_Aphis_craccivora_BOLD", 0,[13,15],[21]),
        ("ACEA833-14_Aphis_spiraecola_BOLD", 0,[13,15],[2,21]),
        ("ACEA589-14_Aphis_spiraecola_BOLD", 0,[13,15],[2,21]),
        ("GBBSP293-15_Philodromus_cespitum_BOLD", 0,[0,12],[2,5]),
        ("NLARA172-12_Philodromus_cespitum_BOLD", 0,[0,12],[2]),
        ("NLARA068-12_Philodromus_cespitum_BOLD", 0,[0,3,12],[2]),
        ("GBBSP2126-16_Philodromus_cespitum_BOLD", 0,[0,12],[2]),
        ("XXX-99_Pilophorus_perplexus_OWN", 11,[0,4,12,15],[2]),
        ("SPSLO262-12_Macaroeris_nidicolens_BOLD", 0,[2,12],[2]),
        ("GCOL11562-16_Scymnus_interruptus_BOLD", 0,[9,13,18],[5,14,23]),
        ("SPSLO172-12_Theridion_pinastri_BOLD", 0,[12,13,15],[11,20]),
        ("SPSLO351-13_Theridion_pinastri_BOLD", 0,[12,13,15],[11,20]),
        ("FBARB426-11_Theridion_varians_BOLD", 0,[0,6,12,13,15],[2]),
        ("SPSLO173-12_Theridion_varians_BOLD", 0,[0,6,12,13,15,18],[2,5]),
        ("XXX-99_Theridion_OWN", 18, [0,6,12,13,15,18], [2,5]),
        ("XXX-99_Chrysoperla_OWN", 8, [0,18], [5]),
        ("XXX-99_Rodolia_cardinalis_OWN", 11, [9,13], [2,14,23]),
        ("XXX-99_Thripidae_OWN", 11, [2], [5]),
        ("XXX-99_Scymnus_subvillosus_OWN", 11, [13], [23]),
        ("XXX-99_Trichopsocus_clarus_OWN", 18, [], [8]),
        ("XXX-99_Ceratitis_capitata_OWN", 11, [], []),
        ("GCOL10228-16_Adalia_decempunctata_BOLD", 0, [13,15], []),
        ("XXX-99_Adalia_decempunctata_OWN", 26, [13,15], []),
        ("GBBSP350-15_Platnickina_tincta_BOLD", 0, [0,12,13], [23]),
        ("TURAR1541-10_Platnickina_tincta_BOLD", 0, [0,12,13], [23]),
        ("NLARA078-12_Platnickina_tincta_BOLD", 0, [0,12,13], [23]),
        ("GBBSP1961-15_Platnickina_tincta_BOLD", 0, [0,12,13], [23])]

def validation_test1():
    gen_record = ld.load_bio_files(["Data/species_bold_own_genbank.fasta"], writable=True)
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")
    primer = primer_pairs[4] #Zeale
    
    gen_alignment_list = []
    
    test_passed = 1
    for t in target1:
        print(t[0])
        gen = gen_record.get(t[0])
        gen_alignment = GenAlignment(gen)

        alignment_list = m.compute_primer_pair_best_alignment(6, 4, primer, gen, hanging_primers=True)
        al= alignment_list.get_list()
        if (len(al)==1):
            al = al[0]
            check = 1
            if(t[1]!=al.fpos):
                print(t[1], al.fpos)
                check = 0
            if(t[2]!=al.fm_loc):
                print(t[2],al.fm_loc)
                check = 0
            if(t[1]+primer.min_amplicon+primer.flen!=al.rpos):
                print(t[1]+primer.min_amplicon+primer.flen,al.rpos)
                check = 0
            if(t[3]!=al.rm_loc):
                print(t[3], al.rm_loc)
                check = 0
            
            if(check==1):
                print("passed")
            else:
                test_passed = 0
                print("failed")
        elif(len(al)>1):
            test_passed = 0
            print("error")
            for a in al:
                print(a)
    if(test_passed):
        print("SUCCESS!")
    else:
        print("TEST FAILED")
    return

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
    
def single_test():
    gen_record = ld.load_bio_files(["Data/mitochondrion.1.1.genomic.fna"], writable=True)
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")
    gen_alignment_list = []
    with open("Test_data/multiple_alignment_cases.txt") as f:
        content = f.readlines()
        content = [x.strip() for x in content]
    i=0
    while i < len(content):
        gen = {content[i]:gen_record[content[i]]}
        pp = primer_pairs[int(content[i+1])-1]
        gen_alignment_list.extend(m.compute_gen_matching(5, 5, [pp], gen, hanging_primers=True))
        i = i+3
        
    multiple_alignment_list = []
    for gm in gen_alignment_list:
        matching_list = gm.get_matching_list()
        for al_list in matching_list: 
            al = al_list.get_list()
            for a in al:
                multiple_alignment_list.append(a)
                
    header = ["primerPair","fastaid","primerF","primerR","mismFT","mismRT","amplicon", "F_pos", "mismFT_type", "R_pos", "mismRT_type"]
    store_results("Test_data/multiple_alignments.csv", multiple_alignment_list, header)          
        
    return

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
            pp = primer_pairs[int(pp)-1]
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
        
    correct_alignments.to_csv("correct_alignments.csv", index_label="id")
    """
    store_results("Test_data/better_alignments.csv", better_alignment_list, header)
    store_results("Test_data/not_tested_alignments.csv", not_tested_alignment_list, header)
    ld.store_matching_results("Test_data/full_alignments.csv", gen_alignment_list, header)
    """
    
    return

def check_better_alignments():
    gen_keys = ["ref|NC_022449.1|","ref|NC_022472.1|","ref|NC_022671.1|","ref|NC_022680.1|","ref|NC_022682.1|","ref|NC_022710.1|",
                "ref|NC_022922.1|","ref|NC_023088.1|"]
    gen_record = ld.load_bio_files(["Data/mitochondrion.1.1.genomic.fna"])
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")
    pp = primer_pairs[4]
    
    gen = gen_record[gen_keys[0]]
    print(gen.seq[1496:1496+pp.flen])
    print(pp.f.seq)
    print(gen.seq[1496+pp.flen+pp.max_amplicon:1496+pp.flen+pp.max_amplicon+pp.rlen])
    print(pp.r.seq)
    
    return

def check_multiple_alignment():
    fprimer = Seq("ATCG")
    fprimer = SeqRecord(fprimer)
    rprimer = Seq("ACGT")
    rprimer = SeqRecord(rprimer)
    pp = PrimerPair(1, fprimer, rprimer, 4, 4)
    primer_pairs = [pp]
    
    gen = Seq("ATCGWACTACGTATCGWACTACGT")
    gen = SeqRecord(gen)
    gen.id = "gen"
    gen_record={"gen": gen}
    
    result = m.compute_gen_matching(1, 1, primer_pairs, gen_record)
    for gm in result:
        matching_list = gm.get_matching_list()
        for al_list in matching_list:
            print(al_list)
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

def store_results(output_file, alignment_list, header=None):
    """
    Stores alignment results
    @param gen_matching_list list of GenMatching instances
    @return None
    """
    with open(output_file, 'w', newline='') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
        if(header):
            filewriter.writerow(header)
        for a in alignment_list:
            filewriter.writerow(a.get_csv())   
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
    
def new_single_test():
    parameters = [
        ["gen", ["Data/sbog_test.fasta"], "Genome file dir, no support for multiple files in cl", "-gf", "entry"],
        ["primer_pairs", "Data/P&PP.csv", "Primer pairs file dir. A particular header must be used in the file", "-pf", "entry"],
        ["output_file", "output.csv", "Location of the output file", "-o", "entry"],
        ["forward missmatches", 5, "Maximum number of missmatches allowed on forward primer", "-fm", "param"],
        ["reverse missmatches", 5, "Maximum number of missmatches allowed on reverse primer", "-rm", "param"], 
        ["Nend miss.", 3, "Missmatches in the last N positions on forward and in the first N pos. on reverse ", "-nend", "info"],
        ["hanging primers", False, "Primers allowed to match between [0-mf,len(genome)+mr] instead of just between genome's length", "--hanging", "param"],
        ["check_integrity", False, "Checks integrity of gen files, integrity of primer file is always checked", "--checki", "param"],
        ["check_uppercase", False, "Checks that all gens are in upper case, lower case gens will trigger an integrity file", "--checku", "param"],
        ["csv_template", "<No Precomputed Template>", "Precomputed missmatching template", "-i", "entry"]]
    parameters = pd.DataFrame([x[1:] for x in parameters], index = [x[0] for x in parameters], columns=["value", "description", "flag", "type"])
    template = i.compute(parameters)
    try:
        print(template.head())
    except:
        print("Error")
        
def check_uppercase():
    gen_record = ld.load_bio_files(["Data/sbog_test.fasta"], check_uppercase=True) 
    print(gen_record["ACEA1016-14_Aphis_spiraecola_BOLD"])
    
def restore_template():
    gen_record = ld.load_bio_files(["Data/sbog_test.fasta"])
    #gen_record = {"AGB001-11_Salticus_scenicus_BOLD": gen_record["AGB001-11_Salticus_scenicus_BOLD"]};
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")
    #primer_pairs = [primer_pairs[5]]
    template, discarded, rs, cs = m.compute_gen_matching(5, 5, primer_pairs, gen_record, 0) 
    
    header = ["primerPair","fastaid","primerF","primerR","mismFT","mismRT","amplicon", "F_pos", "mismFT_loc", "mismFT_type", 
                                     "mismFT_base", "R_pos"]
    #m.store_matching_results("test1.csv", template, header=TEMPLATE_HEADER)
    #templateR = ld.load_template("test1.csv")
    #templateR, rs, cs = ld.restore_template(templateR, gen_record, primer_pairs, 10)
    
    i.save_matching_info("test", template, discarded, rs, cs, header=TEMPLATE_HEADER)
    return


if(__name__=="__main__"):
    
    time1 = time.time()
    test_all_pairs()
    elapsedTime = ((time.time()-time1))
    print(int(elapsedTime)/60)
    
    #performance_test()
