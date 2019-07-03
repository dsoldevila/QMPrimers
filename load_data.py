#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 12:07:28 2018

@author: David Soldevila
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO 
from Bio.Alphabet import IUPAC
from os import path
import csv
from common import *
import pandas as pd
import numpy as np
import ast

def load_csv_file(file, delimiter=";"):
    """
    This function loads a "Primer" file.
    @returns: List of PrimerPair instances
    """
    pos = {"id": 0, "forwardPrimer": 0, "reversePrimer": 0, "fPDNA": 0, "rPDNA": 0,"ampliconMinLength": 0, "ampliconMaxLength": 0}
    header_len = len(pos)
    primer_dict = {}
    with open(file, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=delimiter)
        headers = next(csvreader)
        if(len(headers) != header_len):
            raise ValueError("Wrong header")
        for i in range(len(headers)): 
            if (headers[i] not in pos):
                raise ValueError("Unknown header "+headers[i])
            pos[headers[i]] = i 
            
        i = 1
        for row in csvreader:
            i+=1
            if(len(row) == header_len):
                fprimer = Seq(row[pos["fPDNA"]], IUPAC.IUPACAmbiguousDNA())
                fprimer = SeqRecord(fprimer)
                fprimer.id = row[pos["forwardPrimer"]]
                
                rprimer = Seq(row[pos["rPDNA"]], IUPAC.IUPACAmbiguousDNA())
                rprimer = SeqRecord(rprimer)
                if(True): #TODO
                    rprimer = rprimer.reverse_complement()
                rprimer.id = row[pos["reversePrimer"]]
                
                primer_pair = PrimerPair((row[pos["id"]]), fprimer, rprimer, int(row[pos["ampliconMinLength"]]), 
                                         int(row[pos["ampliconMaxLength"]]))
                if(check_primer_pair_integrity(primer_pair)):
                    primer_dict[row[pos["id"]]] = primer_pair
                else:
                    logging.warning("Skipping primer pair "+primer_pair.id+", bad sequence")
            else:
                logging.warning("Wrong primer pair in line "+str(i))
            
    return primer_dict

def load_bio_files(files, writable=False, check_uppercase=False):
    """
    This function loads any file whose format is supported by Biopython.
    @return: Dictionary with genomic sequences
    """
    #TODO check all files, not only de first one
    if(check_uppercase==True): writable=True #to modify seqrecord, make it writable is needed
    seq_record = {}
    if(type(files)==str):
        files = [files]
    if(path.isfile(files[0])):
        if(writable or len(files)>1): #if they need to be writable, store in memory
            for file in files:
                seq_record.update(SeqIO.to_dict(SeqIO.parse(file, "fasta")))
        else: #create read_only database
            if(isinstance(files, tuple)):
                files = files[0]
            seq_record = SeqIO.index_db(":memory:", files, "fasta") #TODO specify alphabet? It seems it's only used to catch methodology erros
        if(check_uppercase):
            for kseq in seq_record:
                seq_record[kseq] = seq_record[kseq].upper()
    return seq_record

def check_primer_pair_integrity(primer_pair):
    for nuc in primer_pair.f:
        if(nuc not in IUPAC_AMBIGUOUS_DNA):
            return False
    for nuc in primer_pair.r:
        if(nuc not in IUPAC_AMBIGUOUS_DNA):
            return False
    
    return True
    

def remove_bad_gens(gen_record):
    bio_default = {"id":"<unknown id>", "name":"<unknown name>", "description":"<unknown description>"}
    
    genkeys = gen_record.keys()
    for genkey in genkeys:
        gen = gen_record[genkey]
        if(gen.id == None or gen.id == bio_default["id"]): 
             logging.warning(gen.id+"has no id")
        if(gen.name == None or gen.name == bio_default["name"]):
            logging.warning(gen.id+"has no name")
        if(gen.description == None or gen.description == bio_default["description"]):
            logging.warning(gen.id+"has no name")
        for nucleotide in gen.seq:
            if nucleotide not in IUPAC_AMBIGUOUS_DNA:
                logging.warning(gen.id+" HAS BAD SEQUENCE!, trying to remove it...")
                try:
                    gen_record.pop(genkey, None)
                except:
                    logging.error("Couldn't be removed, it will be skipped during match")
                break
            
    return gen_record

def load_template(template_file):
    template = pd.read_csv(template_file, index_col=0)
    
    return template

def restore_template(template, gen_record, primer_pairs, max_misses):
    """
    @brief
    """
    
    if(template.shape[1]==len(TEMPLATE_HEADER)): #the template is complete, no need to redo all the work
        column_is_list = []
        for col in template.columns.values:
            if(type(template.loc[0, col])==str and template.loc[0, col]=="["): #if this column is containing a list
                column_is_list.append(col)
        
        size = template.shape[0]
        for i in range(size):
            for col in column_is_list:
                template.loc[i, col] = ast.literal_eval(template.loc[i,col])
            print("Loading template "+"{0:.2f}".format(i/size*100)+"%")
            
        return template, pd.DataFrame(), pd.DataFrame(), pd.DataFrame()                     
    else:
        alignment = Alignment(max_misses)
        recovered_template = pd.DataFrame(columns=TEMPLATE_HEADER)
        discarded = pd.DataFrame() # TODO
        columns = template.columns.values
        
        #resize table
        for header in TEMPLATE_HEADER:
            if header not in columns:
                template[header] = None
                
        #reorder columns
        template = template[TEMPLATE_HEADER]
        size = template.shape[0]
        for i in range(size):
            tmp = template.loc[i]
            gen = gen_record[tmp.loc["fastaid"]]
            primer_pair = primer_pairs[str(tmp.loc["primerPair"])]
            fpos = tmp.loc["F_pos"]
            rpos = tmp.loc["R_pos"]
            fmisses = tmp.loc["mismFT"]
            rmisses = tmp.loc["mismRT"]
            amplicon = tmp.loc["amplicon"]
            alignment.complete_from_csv(gen, primer_pair, fpos, rpos, fmisses, rmisses, amplicon)
            recovered_template.loc[recovered_template.shape[0]] = alignment.get_csv()
            print("Restoring template "+"{0:.2f}".format(i/size*100)+"%")
        raw_stats, cooked_stats = alignment.get_stats()
    
        return recovered_template, discarded, raw_stats, cooked_stats


if (__name__=="__main__"):
    gen_record = load_bio_files(["Data/mitochondrion.1.1.genomic.fna"], writable=True)
    print(gen_record["ref|NC_012975.1|"].seq)
