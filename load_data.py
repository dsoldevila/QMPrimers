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


def load_csv_file(file, delimiter=";"):
    """
    This function loads a "Primer" file.
    @returns: List of PrimerPair instances
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
            if(True): #TODO
                rprimer = rprimer.reverse_complement()
            rprimer.id = row[pos["reversePrimer"]]
            
            primer_pair = PrimerPair(int(row[pos["id"]]), fprimer, rprimer, int(row[pos["ampliconMinLength"]]), int(row[pos["ampliconMinLength"]]))

            primer_list.append(primer_pair)
            
    return primer_list

def load_bio_files(files, file_format=None, writable=False):
    """
    This function loads any file whose format is supported by Biopython.
    @return: Dictionary with genomic sequences
    """
    #TODO check all files, not only de first one
    seq_record = {}
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
    
    return seq_record

def store_matching_results(output_file, gen_matching_list, header=None):
    """
    Stores alignment results
    @param gen_matching_list list of GenMatching instances
    @return None
    """
    with open(output_file, 'w', newline='') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
        if(header):
            filewriter.writerow(header)
        for gm in gen_matching_list:
            gm.write2file(filewriter)   
    return
    

if (__name__=="__main__"):
    gen_record = load_bio_files(["Data/mitochondrion.1.1.genomic.fna"], writable=True)
    print(gen_record["ref|NC_012975.1|"].seq)
