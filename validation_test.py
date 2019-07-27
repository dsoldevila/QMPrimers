#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 16:58:39 2019

@author: david

This scripts aims to automatize the tests. Running this file with no errors should mean that the program can be released.
"""

import pandas as pd
from common import *
import simulation as s
import interface as ie
import os


matching_param = [
        ["gen", None, "Genome file dir, use ',' to append multiple files", "-gf", "entry"],
        ["primer_pairs", None, "Primer pairs file dir. A particular header must be used in the file", "-pf", "entry"],
        ["output_file", os.path.join(os.getcwd(),"Test_data/output"), "Location of the output files, no extension", "-o", "entry"],
        ["forward missmatches", 5, "Maximum number of missmatches allowed on forward primer", "-fm", "param"],
        ["reverse missmatches", 5, "Maximum number of missmatches allowed on reverse primer", "-rm", "param"], 
        ["Nend miss.", 0, "Missmatches in the last N positions on forward and in the first N pos. on reverse ", "-nend", "info"],
        ["hanging primers", False, "Primers allowed to match between [0-mf,len(genome)+mr] instead of just between genome's length", "--hanging", "param"],
        ["check_integrity", False, "Checks integrity of gen files, integrity of primer file is always checked", "--checki", "param"],
        ["check_uppercase", False, "Checks that all gens are in upper case, lower case gens will trigger an integrity file", "--checku", "param"],
        ["csv_template", None, "Precomputed missmatching template", "-i", "entry"],
        ["verbose", False, "Outputs extra information", "--v", "cmd"]]                      
matching_param  = pd.DataFrame([x[1:] for x in matching_param ], index = [x[0] for x in matching_param ], columns=["value", "description", "flag", "type"])

sim_param = [
        ["template", None, "Precomputed missmatching template", "-i", "entry"],
        ["sample size", 10, "Nº genome samples per simulation step", "-s", "int"],
        ["Beta", 4, "<No description currently>", "-b", "int"],
        ["k", 0.5, "Geometric proportion parameter. Range(0,1)", "-k", "float"],
        ["N", 100, "Nº simulation steps", "-n", "int"],
        ["output_file", os.path.join(os.getcwd(),"simout"), "Location of the output files, no extension", "-o", "entry"],
        ["Confidence Interval", 0.95, "CI used to compute the statistics", "-ci", "float"],
        ["verbose", False, "Outputs extra information", "--v", "cmd"]]

sim_param  = pd.DataFrame([x[1:] for x in sim_param ], index = [x[0] for x in sim_param ], columns=["value", "description", "flag", "type"])


def corrupt_input_test():
    pass


def check_matching(template, correct_template):
    return 0

def check_nend(template, correct_template):
    return 0

def matching_test(gen_files, primer_files):
    #fer matching
    #guardar template
    #carregar template
    #comparar template
    #calcular nend_template
    #comparar nend_template
    
    results = pd.DataFrame(columns=["exec", "save", "load", "match_check", "nend", "nend_check"])
    
    logging.debug("MATCHING TEST...")
    header = [i for i in range(len(TEMPLATE_HEADER))]
    
    input_dir = os.path.join(os.getcwd(), "test_input")
    od = os.path.join(os.getcwd(),"test_output")
    vd = os.path.join(os.getcwd(),"test_validated")
    
    for i in range(len(gen_files)):
        d = os.path.dirname(gen_files[i])
        matching_param.loc["gen", "value"] = os.path.join(input_dir, gen_files[i])
        matching_param.loc["primer_pairs", "value"] = os.path.join(input_dir,primer_files[i])
        matching_param.loc["output_file", "value"] = os.path.join(od, gen_files[i]+primer_files[i])
        
        #MATCHING
        logging.debug("Starting to compute "+gen_files[i]+" "+primer_files[i])
        try:
            template, discarded, gen_record, primer_pairs, raw_stats, cooked_stats = ie.compute(matching_param) 
            results.loc[gen_files[i]+primer_files[i], "exec"] = "OK"
        except(Exception) as e:
            logging.error("Match "+gen_files[i]+" "+primer_files[i])
            logging.error(e)
            results.loc[gen_files[i]+primer_files[i], "exec"] = "FAIL"
        
        #SAVING
        logging.debug("Saving "+gen_files[i]+" "+primer_files[i])
        try:
            input_files = "Fasta = "+str(matching_param.loc["gen", "value"])+"     Primer Pairs = "+matching_param.loc["primer_pairs", "value"]
            ie.save_matching_info(input_files, matching_param.loc["output_file", "value"], template, header, discarded, raw_stats, cooked_stats)
            results.loc[gen_files[i]+primer_files[i], "save"] = "OK"
        except:
            results.loc[gen_files[i]+primer_files[i], "exec"] = "FAIL"
         
        #LOADING
        logging.debug("Loading "+gen_files[i]+" "+primer_files[i])
        try:
            matching_param.loc["csv_template", "value"] = matching_param.loc["output_file", "value"]+"_positive.csv"
            template = ie.load_template(matching_param)
            results.loc[gen_files[i]+primer_files[i], "load"] = "OK"
        except:
            logging.error("Load template "+gen_files[i]+" "+primer_files[i])
            results.loc[gen_files[i]+primer_files[i], "load"] = "FAIL"
            
        #CHECKING MATCH
        logging.debug("Checking results "+gen_files[i]+" "+primer_files[i])
        try:
            matching_param.loc["csv_template", "value"] = matching_param.loc["output_file", "value"]+"_positive.csv"
            template = ie.load_template(matching_param)[0]
            matching_param.loc["csv_template", "value"] = os.path.join(vd, gen_files[i]+primer_files[i]+".csv")
            validated_template = ie.load_template(matching_param)[0]
            if(check_matching(template, validated_template)):
                results.loc[gen_files[i]+primer_files[i], "match_check"] = "OK"
            else:
                results.loc[gen_files[i]+primer_files[i], "match_check"] = "FAIL"
        except(Exception) as e:
            logging.error(e)
            logging.error("Template "+gen_files[i]+primer_files[i]+" has no validated template to compare")
            results.loc[gen_files[i]+primer_files[i], "match_check"] = "FAIL"
            
        #NEND
        logging.debug("Computing Nend "+gen_files[i]+" "+primer_files[i])
        try:
            nend = 5
            template, raw_stats, cooked_stats = ie.get_Nend_match(template, nend)
            output = os.path.join(od, "n"+str(nend)+gen_files[i]+primer_files[i])
            ie.save_matching_info(input_files, output, template, header, discarded, raw_stats, cooked_stats)
            results.loc[gen_files[i]+primer_files[i], "nend"] = "OK"
        except(Exception) as e:
            logging.error(e)
            results.loc[gen_files[i]+primer_files[i], "nend"] = "FAIL"
        
        #NEND CHECK
        logging.debug("Checking Nend "+gen_files[i]+" "+primer_files[i])
        try:
            matching_param.loc["csv_template", "value"] = os.path.join(od, "n"+str(nend)+gen_files[i]+primer_files[i]+"_positive.csv")
            template = ie.load_template(matching_param)[0]
            matching_param.loc["csv_template", "value"] = os.path.join(vd, "n"+str(nend)+gen_files[i]+primer_files[i]+".csv")
            validated_template = ie.load_template(matching_param)[0]
            if(check_nend(template, validated_template)):
                results.loc[gen_files[i]+primer_files[i], "nend_check"] = "OK"
            else:
                results.loc[gen_files[i]+primer_files[i], "nend_check"] = "FAIL"
        except:
            logging.error("Template "+gen_files[i]+primer_files[i]+" has no validated nend template to compare")
            results.loc[gen_files[i]+primer_files[i], "nend_check"] = "FAIL"
            
        print(results)
    
    
    return

def simulation_test():
    template_file = "test_output/insects_positive.csv"
    sample_size = 10
    k = 0.3
    B = 4
    N=10
    ci = 0.90
    template = ie.load_template_only(template_file)
    sim = s.Simulation()
    raw, cooked = sim.simulate(template, sample_size, k, B, N, ci)
    sim.store_data("test_output/sim_test", raw, cooked, template_file, sample_size, k, B, N)
    pass

def get_matching_input():
    #input_dir = os.path.join(os.getcwd(), "test_input")
    gen_files = ["sbog_test.fasta"]
    """
    for i in range(len(gen_files)):
        gen_files[i] = os.path.join(input_dir, gen_files[i])
    """   
    primer_files = ["PP.csv"]
    
    """
    for i in range(len(primer_files)):
        primer_files[i] = os.path.join(input_dir, primer_files[i])
    """
    return gen_files, primer_files
        
    

if(__name__=="__main__"):
    init_logger()
    set_verbosity(2) #debug mode
    
    gen_files, primer_files = get_matching_input()
    matching_test(gen_files, primer_files)
    #simulation_test()
    
    
    close_logger