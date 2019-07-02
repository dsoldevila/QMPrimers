#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 17:02:08 2018

@author: david
"""
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from common import *

import numpy as np
import pandas as pd
import load_data as ld
import queue
import threading
import time

class matching_thread_wrapper(threading.Thread):
    def __init__(self, mosi_queue, miso_queue):
        """
        @param mosi_queue: Master Out Slave In --> queue.Queue() from GUI to this thread
        @param miso_queue: Master In Salve Out --> queue.Queue() from this thread to GUI
        """
        self.mosiq = mosi_queue
        self.misoq = miso_queue
        threading.Thread.__init__(self)
        self._stop_event = threading.Event()
        
    def run(self):
        """
        State machine of the matching backend. It keeps waiting for parameters in the queue until a exit command is received
        """
        keep_running = True
        command = "none"
        while(keep_running):
            if(self.mosiq.qsize()>0):
                command = self.mosiq.get(0)
                if(command == "exit"):
                    keep_running = False
                elif(command == "match"):
                    self.match()
                elif(command == "save"):
                    self.save()
                elif(command!="none"):
                    with self.mosiq.mutex: self.mosiq.queue.clear()
                    self.misoq.put(False) #Send error to the master
            logging.debug("Matching thread state: "+command)
            time.sleep(1)
        self.stop()
        self.misoq.put(self.stopped())
                    
    def match(self):
        logging.debug("Matching thread: Entered match state")
        try:
            max_miss_f = self.mosiq.get()
            max_miss_r = self.mosiq.get()
            self.pp_file = self.mosiq.get()
            self.gen_file = self.mosiq.get()
            self.output_file = self.mosiq.get()
            self.check_integrity = self.mosiq.get()
            self.check_uppercase = self.mosiq.get()
            hanging_primers = self.mosiq.get()
            self.misoq.put(True)
        except:
            logging.error("Matching could not retrieve input data")
            self.misoq.put(False)
            return
        logging.debug("Matching thread: Matching paramaters read")
        
        logging.debug("Matching thread: Getting data")
        self.gen_record, self.primer_pairs = self.load_data()
            
        logging.debug("Matching thread: Starting matching..")
        try:
            self.template, self.discarded, self.raw_stats, self.cooked_stats = compute_gen_matching(max_miss_f, max_miss_r, 
                                                              self.primer_pairs, self.gen_record, self.output_file, hanging_primers=hanging_primers)
            self.out_template = self.template
            self.out_raw_stats = self.raw_stats
            self.out_cooked_stats = self.cooked_stats
        except:
            logging.error("Matching thread: Matching crashed")
            self.misoq.put(False)
            return
        logging.debug("Matching thread: Matching done")
        
        self.misoq.put(True)
        self.misoq.put(self.template)
        self.misoq.put(self.discarded)
        self.misoq.put(self.raw_stats)
        self.misoq.put(self.cooked_stats)
        self.misoq.put(self.gen_record)
        self.misoq.put(self.primer_pairs)
        logging.debug("Matching thread: Data written into queue")
            
        
        return
        
    def load_data(self):
        """
        Matching should not care about the looad data options, but it is the siplest way to do it in order to not freeze the GUI with large files.
        """
        logging.debug("Matching thread: Entered load_data function")
        gen_record = ld.load_bio_files(self.gen_file, self.check_integrity, self.check_uppercase)
        primer_pairs = ld.load_csv_file(self.pp_file)
        
        return gen_record, primer_pairs
        
    def save(self):
        logging.debug("Matching thread: Entered save function")
        try:
            self.output_file = self.mosiq.get()
            header = self.mosiq.get()
            nend = self.mosiq.get()
        except:
            logging.error("Saving could not retrieve input data")
            return
        
        if(nend):
            if(nend!=self.previous_nend):
                self.match_nend()
            store_matching_results(self.output_file, self.out_template, header)
            store_stats(self.output_file, self.out_raw_stats, self.out_cooked_stats)
        else:
            store_matching_results(self.output_file, self.template, header)
            store_stats(self.output_file, self.raw_stats, self.cooked_stats)
            
        store_discarded(self.output_file, self.discarded)
        return
        
    def match_nend(self):
        try:
            nend = self.mosiq.get()
            max_misses = self.mosiq.get()
        except:
            logging.error("Matching N-end could not retrieve input data")
            return
        
        self.out_template, self.out_raw_stats, self.out_cooked_stats = get_Nend_template(self.template, nend, max_misses)
        #self.misoq.put(self.out_template)
        #self.misoq.put(self.out_raw_stats)
        #self.misoq.put(self.out_cooked_stats)
        return
    
    def get_data(self):
        #put data into the miso queue
        pass
        
    def stop(self):
        self._stop_event.set()

    def stopped(self):
        return self._stop_event.is_set()
        
                
        

#TODO this func is currently unused, ~line 68
def is_valid(score, min_score):
    score = np.asarray(score)
    min_score = np.asarray(min_score)
    return np.where(score >= min_score, score, 0)

def append_zeros(gen_record, max_miss_f, max_miss_r):
    """
    Adds zeros ('Z') at the beginning and end of the genome to allow the primer to pass the genome's limits
    Example:
        genome AGATTCATT, primer TCTAGA scores 2 at pos 1
        genome ZZZAGATTCATTZZZ, primer TCTAGA scores 3 at pos -3
    """
    for gen_key in gen_record:
        gen_record[gen_key].seq = Seq("Z"*max_miss_f+str(gen_record[gen_key].seq)+"Z"*max_miss_r)
    return gen_record


def _compute_primer_matching(max_misses, primer, len_primer, gen):
    """
    Computes the best matches between a genome and a primer.
    @returns: Numpy matrix of arrays (score, start_pos, end_pos).
    """
    result_matrix = MATCH_TABLE.loc[primer, gen] #get match table

    result_max_len = len(gen)-len_primer+1
    result_raw = np.zeros(result_max_len, dtype='uint8') #TODO 0-255 should be enough, but better to not hardcode this

    result_matrix = result_matrix.values #get numpy matrix, raw indexes are faster than labels
    for i in range(len_primer):
        result_raw = np.add(result_raw, result_matrix[i,i:result_max_len+i]) #Speedup try n2, SpeedUp = 168/26 = 6,4x
    
    is_score_valid = (result_raw>=len_primer-max_misses)
    n_results = is_score_valid.sum()
    result_raw = is_score_valid*result_raw #if score valid =score else =0
    
    result = np.zeros(n_results, dtype=[('score', 'uint8'), ('start', '>i4'), ('end', '>i4')]) #TODO integer 32 too much /  integer 8 enough?
    
    j = 0
    for i in range(result_max_len):        
        if result_raw[i]:
            result[j] = (result_raw[i], i, i+len_primer)
            j=j+1

    return result

def compute_primer_pair_best_alignment(max_miss_f, max_miss_r, primer, gen, hanging_primers, template, discarded, alignment_processor):
    """
    Returns the best alignments between a genome and a primer pair
    @returns: PrimerAlignment instance
    """
    
    max_amplicon = primer.max_amplicon
    search_limit = primer.rlen+max_amplicon
    len_gen = len(gen)
    if(primer.flen+search_limit>len_gen): #If primer pair plus max_amplicon is larger than the genomic sequence, check if with min_amplicon the same happens
        if(primer.flen+primer.rlen+primer.min_amplicon>len_gen): #If primer pair plus min_amplicon is larger than the genomic sequence, abort
            logging.warning("Skipping gen "+gen.id+" primer pair "+str(primer.id))
            discarded.loc[discarded.shape[0]] = [primer.id, gen.id]
            alignment_processor.add_negative_2_stats(primer.id)
            return template, discarded
        else: #else modify max_amplicon to keep the primer_pair within the limits
            max_amplicon = len_gen - (primer.flen + primer.rlen)
    
    best_score = 0
    alignments = []

    forward_matchings = _compute_primer_matching(max_miss_f, primer.f.seq, primer.flen, gen.seq[0:-search_limit]) #compute forward primer best matches
    for fm in forward_matchings: #for each match with forward primer, compute reverse matchings
        start = fm[2]+primer.min_amplicon #forward match start + len(forward) + min amplicon
        end = fm[2]+max_amplicon+primer.rlen #f match start + len(f) + max amplicon + len(r)
        reverse_matchings = _compute_primer_matching(max_miss_r, primer.r.seq, primer.rlen, gen.seq[start:end])
        for rm in reverse_matchings: #get the best or bests matche(s) with this primer pair (alingments)
            score = fm[0] + rm[0]
            if (score > best_score): #if the score is better, erase the previous bests results
                alignments = [(fm, rm)]
                best_score = score
            elif (score == best_score): #elif the score is equaly good, get this alignment too
                alignments.append((fm, rm))
                
    if(alignments==[]):
        discarded.loc[discarded.shape[0]] = [primer.id, gen.id]
        alignment_processor.add_negative_2_stats(primer.id)
    for al in alignments:
        fm = al[0]
        rm= al[1]
        amplicon = primer.min_amplicon+rm[1]
        alignment_processor.get(gen, primer, fm[1], fm[1]-max_miss_f*hanging_primers, rm[1]+amplicon+fm[2], rm[1]+amplicon+fm[2]-max_miss_f*hanging_primers,
                                primer.flen-fm[0], primer.rlen-rm[0], amplicon) #TODO, creating a temp class overkill?
        template.loc[template.shape[0]] = alignment_processor.get_csv()
        
    return template, discarded

def compute_gen_matching(max_miss_f, max_miss_r, primer_pairs, gen_record, output_file, hanging_primers=False):
    """
    @brief Computes the best alignments between each genome and each primer
    @Returns template: pandas DataFrame containing the matching results
    @Returns discarded: pandas DataFrame containing pairs of primer pairs and genome sequence that have not matched
    @Returns raw_stats: pandas Dataframe containing a summary of the matching
    @Returns cooked_stats: pandas Dataframe containing statictics the matching
    """
    if(hanging_primers):
        gen_record = append_zeros(gen_record, max_miss_f, max_miss_r)
    
    primerPair_list = [] #used to sort template by primer pairs
    for pkey in primer_pairs:
        pp = primer_pairs[pkey]
        pp.f.seq = np.array(pp.f)
        pp.r.seq = np.array(pp.r)
        primerPair_list.append(pkey)
        
    size = len(gen_record)
    i = 0       
    template_header = TEMPLATE_HEADER

    template = pd.DataFrame(columns=template_header)
    
    discarded = pd.DataFrame(columns=TEMPLATE_HEADER[0:2])

    alignment_processor = Alignment(max_miss_f + max_miss_r)
    
    template.to_csv(output_file+"_positive.csv", index_label="id")
    discarded.to_csv(output_file+"_negative.csv", index_label="id")
    
    t1 = 0
    d1 = 0
    for gen_key in gen_record:
        print(gen_key, "{0:.2f}".format(i/size*100)+"%")
        i +=1
        gen = gen_record[gen_key]
        gen.seq = np.array(gen.seq)
        for pkey in primer_pairs:
            try:
                t2 = t1
                d2 = d1
                template, discarded = compute_primer_pair_best_alignment(max_miss_f, max_miss_r, primer_pairs[pkey], gen, hanging_primers, template, discarded, alignment_processor)
                t1 = template.shape[0]
                d1 = discarded.shape[0]
                template.loc[template.index[t2:t1]].to_csv(output_file+"_positive.csv", mode='a', index_label="id", header=None)
                discarded.loc[discarded.index[d2:d1]].to_csv(output_file+"_negative.csv", mode='a', index_label="id", header=None)
            except:
                raise
                logging.error("Skipping gen "+gen.id+" primer pair "+str(pkey))
   
    raw_stats, cooked_stats = alignment_processor.get_stats()
    store_stats(output_file+"_stats.txt", raw_stats, cooked_stats)
    print("Template, negative and statistics saved")
    
    
    template["primerPair"] = pd.Categorical(template["primerPair"], categories=primerPair_list, ordered=True)
    template.sort_values(['primerPair', 'fastaid'], inplace=True)
    template.reset_index(drop=True, inplace=True)
    
    discarded["primerPair"] = pd.Categorical(discarded["primerPair"], categories=primerPair_list, ordered=True)
    discarded.sort_values(['primerPair', 'fastaid'], inplace=True)
    discarded.reset_index(drop=True, inplace=True)
         
    
    return template, discarded, raw_stats, cooked_stats

"""
OTHER FUNCTIONS
"""

def store_matching_results(output_file, template, header):
    """
    Stores alignment results
    @param gen_matching_list list of GenMatching instances
    @return None
    """
    try:
        columns = template.columns.values
        for i in range(len(header)):
            header[i] = columns[header[i]]
        template.to_csv(output_file, index_label="id", columns=header)
        print("Positive results saved")
    except:
        logging.error("Could not save positive results")
    return

def store_stats(output_file, raw_stats, cooked_stats):
    try:
        with open(output_file,'w') as outfile:
            raw_stats.to_string(outfile)
            outfile.write("\n\n")
            cooked_stats.to_string(outfile)
            print("Statistics saved")
    except:
        logging.error("Could not save statistics")
    return

def store_discarded(output_file, discarded):
    try:
        discarded.to_csv(output_file, index_label="id")
        print("Negative results saved")
    except:
        logging.error("Could not save negative results")
    return


def get_Nend_template(template, nend, max_misses):
    """
    @brief With the computed template, generate a template but with Nend mismatches
    @Return new template, raw stats, cooked_stats
    """
    try:
        alignment = Alignment(max_misses);
        
        header = ["primerPair","fastaid","primerF","primerR","mismFN"+str(nend),"mismRN"+str(nend),"amplicon", "F_pos", 
                  "mismFN"+str(nend)+"_loc", "mismFN"+str(nend)+"_type", "mismFN"+str(nend)+"_base", "R_pos", "mismRN"+str(nend)+"_loc",
                  "mismRN"+str(nend)+"_type", "mismRN"+str(nend)+"_base"]
        
        nend_template = pd.DataFrame(columns=header)
        size = template.shape[0]
        for i in range(template.shape[0]):
            nend_template.loc[i] = alignment.get_Nend(*template.loc[i], nend)
            print("Getting Nend "+"{0:.2f}".format(i/size*100)+"%")
        raw_stats, cooked_stats = alignment.get_stats()
    except:
        logging.error("Crash when trying to recover N-end template")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    
    return nend_template, raw_stats, cooked_stats
    

if(__name__=="__main__"):
    
    gen_record = ld.load_bio_files(["Data/species_bold_own_genbank.fasta"])
    primer_pairs = ld.load_csv_file("Data/P&PP.csv")


    import time 
    time1 = time.time()
    result = compute_gen_matching(5, 5, primer_pairs, gen_record)
    elapsedTime = ((time.time()-time1))
    print(int(elapsedTime))

    