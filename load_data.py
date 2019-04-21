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


def load_csv_file(file, delimiter=";"):
    """
    This function loads a "Primer" file.
    @returns: List of PrimerPair instances
    """
    pos = {"id": 0, "forwardPrimer": 0, "reversePrimer": 0, "fPDNA": 0, "rPDNA": 0,"ampliconMinLength": 0, "ampiconMaxLength": 0}
    header_len = len(pos)
    primer_list = []
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
                
                primer_pair = PrimerPair((row[pos["id"]]), fprimer, rprimer, int(row[pos["ampliconMinLength"]]), int(row[pos["ampliconMinLength"]]))
                if(check_primer_pair_integrity(primer_pair)):
                    primer_list.append(primer_pair)
                else:
                    print("Error: Skipping primer pair "+primer_pair.id+", bad sequence")
            else:
                print("Wrong primer pair in line "+str(i))
            
    return primer_list

def load_bio_files(files, file_format=None, writable=False, check_uppercase=False):
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
        if(file_format==None): #if format not specified, use extension
            extension = path.splitext(files[0])[1]
            file_format = extension[1:]
            #TODO change this patch
            if(file_format == "fna"): file_format = "fasta"
            #TODO is the next line acceptable?
        if(file_format in SeqIO._FormatToIterator): #if format supported by biopython
            if(writable or len(files)>1): #if they need to be writable, store in memory
                for file in files:
                    seq_record.update(SeqIO.to_dict(SeqIO.parse(file, file_format)))
            else: #create read_only database
                if(isinstance(files, tuple)):
                    files = files[0]
                seq_record = SeqIO.index_db(":memory:", files, file_format) #TODO specify alphabet? It seems it's only used to catch methodology erros
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
             print("Warning: "+gen.id+"has no id")
        if(gen.name == None or gen.name == bio_default["name"]):
            print("Warning: "+gen.id+"has no name")
        if(gen.description == None or gen.description == bio_default["description"]):
            print("Warning: "+gen.id+"has no name")
        for nucleotide in gen.seq:
            if nucleotide not in IUPAC_AMBIGUOUS_DNA:
                print("Error: "+gen.id+" HAS BAD SEQUENCE!, trying to remove it...")
                try:
                    gen_record.pop(genkey, None)
                except:
                    print("Warning: "+gen.id+" couldn't be removed, it will be skipped during match")
                break
            
    return gen_record

def load_template(template_file):
    template = pd.read_csv(template_file)
    return template

def restore_template(template, gen_record, primer_pairs):
    alignment = Alignment()
    columns = template.columns.values
    
    #resize table
    for header in TEMPLATE_HEADER:
        if header not in columns:
            template[header] = None
            
    #reorder columns
    template = template[TEMPLATE_HEADER]
    
    for i in range(template.shape[0]):
        gen = gen_record[template.loc[i, "fastaid"]]
        primer_pair = primer_pairs[int(template.loc[i, "primerPair"])-1]
        fpos = template.loc[i, "F_pos"]
        rpos = template.loc[i, "R_pos"]
        fmisses = template.loc[i, "mismFT"]
        rmisses = template.loc[i, "mismRT"]
        amplicon = template.loc[i, "amplicon"]
        alignment.complete_from_csv(gen, primer_pair, fpos, rpos, fmisses, rmisses, amplicon)
        template.loc[i] = alignment.get_csv()
        
    return template
    

if (__name__=="__main__"):
    gen_record = load_bio_files(["Data/mitochondrion.1.1.genomic.fna"], writable=True)
    print(gen_record["ref|NC_012975.1|"].seq)
