#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 11:37:52 2019

@author: david

The purpose of this code is to check if this program is working properly
"""
import load_data as ld
from common import *
from matching import *

import cProfile
import pstats

target1 = [("XXX-99_Xysticus_OWN", 11, [0,12,13,15],157,[2,5]),
        ("FBARB265-11_Clubiona_leucaspis_BOLD", 0,[0,12],157,[2,5,14,18]),
        ("ACEA902-14_Myzus_persicae_BOLD", 0,[13,15],157,[2,17,21]),
        ("ACEA563-14_Aphis_gossypii_BOLD", 0,[13,15],157,[5,21]),
        ("ACEA640-14_Aphis_craccivora_BOLD", 0,[13,15],157,[21]),
        ("ACEA833-14_Aphis_spiraecola_BOLD", 0,[13,15],157,[2,21]),
        ("ACEA589-14_Aphis_spiraecola_BOLD", 0,[13,15],157,[2,21]),
        ("GBBSP293-15_Philodromus_cespitum_BOLD", 0,[0,12],157,[2,5]),
        ("NLARA172-12_Philodromus_cespitum_BOLD", 0,[0,12],157,[2]),
        ("NLARA068-12_Philodromus_cespitum_BOLD", 0,[0,3,12],157,[2]),
        ("GBBSP2126-16_Philodromus_cespitum_BOLD", 0,[0,12],157,[2]),
        ("XXX-99_Pilophorus_perplexus_OWN", 11,[0,4,12,15],157,[2]),
        ("SPSLO262-12_Macaroeris_nidicolens_BOLD", 0,[2,12],157,[2]),
        ("GCOL11562-16_Scymnus_interruptus_BOLD", 0,[9,13,18],157,[5,14,23]),
        ("SPSLO172-12_Theridion_pinastri_BOLD", 0,[12,13,15],157,[11,20]),
        ("SPSLO351-13_Theridion_pinastri_BOLD", 0,[12,13,15],157,[11,20]),
        ("FBARB426-11_Theridion_varians_BOLD", 0,[0,6,12,13,15],157,[2]),
        ("SPSLO173-12_Theridion_varians_BOLD", 0,[0,6,12,13,15,18],157,[2,5]),
        ("XXX-99_Theridion_OWN", 18, [0,6,12,13,15,18], 157, [2,5]),
        ("XXX-99_Chrysoperla_OWN", 8, [0,18], 157, [5]),
        ("XXX-99_Rodolia_cardinalis_OWN", 11, [9,13], 157, [2,14,23]),
        ("XXX-99_Thripidae_OWN", 11, [2], 157, [5]),
        ("XXX-99_Scymnus_subvillosus_OWN", 11, [13], 157, [23]),
        ("XXX-99_Trichopsocus_clarus_OWN", 18, [], 157, [8]),
        ("XXX-99_Ceratitis_capitata_OWN", 11, [], 157, []),
        ("GCOL10228-16_Adalia_decempunctata_BOLD", 0, [13,15], 157, []),
        ("XXX-99_Adalia_decempunctata_OWN", 26, [13,15], 157, []),
        ("GBBSP350-15_Platnickina_tincta_BOLD", 0, [0,12,13], 157, [23]),
        ("TURAR1541-10_Platnickina_tincta_BOLD", 0, [0,12,13], 157, [23]),
        ("NLARA078-12_Platnickina_tincta_BOLD", 0, [0,12,13], 157, [23]),
        ("GBBSP1961-15_Platnickina_tincta_BOLD", 0, [0,12,13], 157, [23])]

def test1(gen_record, primer_pairs):
    primer = primer_pairs[4] #Zeale
    
    gen_matching_list = []
    
    test_passed = 1
    for t in target1:
        print(t[0])
        gen = gen_record.get(t[0])
        gen_matching = GenMatching(gen)

        alignment_list = compute_primer_pair_best_alignment(6, 4, primer, gen)
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
            if(t[1]+t[3]+21!=al.rpos):
                print(t[1]+t[3]+21,al.rpos)
                check = 0
            if(t[4]!=al.rm_loc):
                print(t[4], al.rm_loc)
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

def performance_test(gen_record, primer_pairs):
    cProfile.run('compute_gen_matching(5, 5, primer_pairs, gen_record)', 'myFunction.profile')
    stats = pstats.Stats('myFunction.profile')
    stats.strip_dirs().sort_stats('time').print_stats()

if(__name__=="__main__"):
    gen_record = ld.load_bio_file("Data/species_bold_own_genbank.fasta")
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")
    performance_test(gen_record, primer_pairs)
    #print(compute_primer_pair_best_alignment(6, 4, primer_pairs[4], gen_record.get("GBBSP1961-15_Platnickina_tincta_BOLD")))