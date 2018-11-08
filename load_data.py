#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 12:07:28 2018

@author: david
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO 
from Bio.Alphabet import IUPAC
from os import path
import csv



def loadCSVFile(file, delimiter=";"):
    """
    This function loads a "Primer" file.
    @returns: List of Tuples
    Tuple = (PairID, ForwardPrimer_seq, ReversePRimer_seq)
    """
    pos = {"id": 0, "forwardPrimer": 0, "reversePrimer": 0, "fPDNA": 0, "rPDNA": 0,"ampliconMinLength": 0, "ampiconMaxLength": 0}
    primer_list = []
    with open(file, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=delimiter)
        headers = next(csvreader)
        for i in range(len(headers)): pos[headers[i]] = i 

        for row in csvreader:
            fprimer = Seq(row[pos["fPDNA"]])
            fprimer = SeqRecord(fprimer)
            fprimer.id = row[pos["forwardPrimer"]]
            
            rprimer = Seq(row[pos["rPDNA"]])
            rprimer = SeqRecord(rprimer)
            rprimer.id = row[pos["reversePrimer"]]
            
            primerPair = (row[pos["id"]], fprimer, rprimer)
            primer_list.append(primerPair)
            
    return primer_list

def loadBioFile(file, file_format=None):
    """
    This function loads any file whose format is supported by Biopython.
    @return: Dictionary with genomic sequences
    """
    seq_record = {}
    if(path.isfile(file)):
        if(file_format==None): #if format not specified, use extension
            extension = path.splitext(file)[1]
            file_format = extension[1:]
            
        if(file_format in SeqIO._FormatToIterator): #if format supported by biopython
            seq_record = SeqIO.index(file, file_format, alphabet=None)
    
    return seq_record
  
    

if "__main__":
    gen_record_list = loadBioFile("species_bold_own_genbank.fasta")
    """for gen in gen_record_list:
        print(gen_record_list.get(gen).id)"""
    gen = gen_record_list.get("ACEA1016-14_Aphis_spiraecola_BOLD")
    example_feature = SeqFeature(FeatureLocation(5, 18), type="LCO1490", strand=1)
    gen.features.append(example_feature)
    print(gen.features)
    loadCSVFile("P&PP.csv")