#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 12:07:28 2018

@author: david
"""
from Bio.Seq import Seq
from Bio import SeqIO 
from Bio.Alphabet import IUPAC
from os import path


def __loadBioFile(file_name, extension, alphabet=None):
    """
    This function loads any file whose format is supported by Biopython
    """
    seq_records = SeqIO.index(file_name, extension, alphabet)
    return seq_records
    
def __loadCSVFile(file):
    pass

def LoadFile(file, file_format=None):
    if(path.isfile(file)):
        if(file_format==None): #if format not specified, use extension
            extension = path.splitext(file)[1]
            file_format = extension[1:]
            
        if(file_format in SeqIO._FormatToIterator): #if format supported by biopython
            __loadBioFile(file, file_format, alphabet=None)
            
        elif(file_format == "csv"):
            __loadCSVFile(file)
    #if supported by biopython
    
    pass
    

if "__main__":
    gen_record_list = __loadBioFile("species_bold_own_genbank.fasta")
    print(gen_record_list[0])
    