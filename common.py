#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 15:34:23 2018

@author: David Soldevila
"""
import numpy as np
import pandas as pd
import logging
import os
import sys
import datetime
import re

#Some global constant variables
IUPAC_AMBIGUOUS_DNA = tuple("ACGTWSMKRYBDHVNIZ")
TEMPLATE_HEADER = ["primerPair","fastaid","primerF","primerR","mismFT","mismRT","ampliconLen", "F_pos", "mismFT_loc", "mismFT_type", 
                                     "mismFT_base", "R_pos", "mismRT_loc", "mismRT_type", "mismRT_base", "amplicon"]

#This matrix tells the algorithm whether 2 nucleotides match or don't
SCORE_TABLE = np.array([[1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0],
                        [0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0],
                        [0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0],
                        [0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0],
                        [1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0],
                        [0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                        [1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0],
                        [0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0],
                        [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
                        [1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0],
                        [1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0],
                        [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype='uint8')
MATCH_TABLE = pd.DataFrame(SCORE_TABLE, index=list("ACGTWSMKRYBDHVNIZ"), columns=list("ACGTWSMKRYBDHVNIZ"))


class PrimerPair:
    def __init__(self, pair_id, fprimer, rprimer, min_amplicon, max_amplicon):
        self.id = pair_id
        self.f = fprimer
        self.flen = len(fprimer.seq)
        self.r = rprimer
        self.rlen = len(rprimer.seq)
        self.min_amplicon = min_amplicon
        self.max_amplicon = max_amplicon
        self.fcomplement = self.f.seq.complement()
        self.rcomplement = self.r.seq.complement()
        return
    
    def __str__(self):
        pass
    
class Alignment:
    """
    Alignment info between a genomic sequence and a primer pair
    """
    base_type = {"A":"Pur", "C":"Pyr", "G":"Pur", "T":"Pyr", "R":"Pur", "Y":"Pyr", "Other": "Ind."}

    def __init__(self, max_misses):
        columns = list(range(max_misses+1))
        self.miss_range = columns.copy()
        columns.append("No")
        columns.append("Total")
        self.amplicon_stats = {}
        self.pp_stats = pd.DataFrame(columns=columns, dtype='uint32')
        return
    
    def get(self, gen, primer_pair, fpos, real_fpos, rpos, real_rpos, fmisses, rmisses, amplicon):
    
        """
        self.gen = genomic sequence
        self.primer_pair = primer pair used for matching, instance of PrimerPair class
        self.F_pos = starting position of the forward primer in the genomic sequence, starting at 0
        self.real_fpos --> gen:    ZZAGTAC...     real_fpos = -2, the primer is hanging
                           primer: AGAGT          fpos = 0
        self.R_pos = starting position of the reverse primer in the genomic sequence, starting at the end of genomic sequence
        self.real_rpos = reverse's position depends on forward's position
        self.mismF = number of missmatches on the forward primer
        self.mismF_loc_raw = array of the missmatch locations of the forward primer
        self.mismR = number of missmatches on the reverse primer
        self.mismR_loc_raw = array of the missmatch locations of the reverse primer
        self.amplicon_len = amplicon of the matching, number between the primer pair max and min amplicon.
        """
        self.gen = gen
        self.fastaid = gen.id #TODO patch
        self.primer_pair = primer_pair
        self.primerPair = primer_pair.id #TODO patch
        self.primerF = primer_pair.f.id
        self.primerR = primer_pair.r.id
        
        self.F_pos = int(fpos) #it seems Biopython seqrecord does not support numpy.int32
        self.real_fpos = int(real_fpos)
        logging.debug("Forward pos: "+str(self.real_fpos))
        self.R_pos = int(rpos)
        self.real_rpos = int(real_rpos)
        logging.debug("Reverse pos: "+str(self.real_rpos))
        self.mismF = fmisses
        self.mismR = rmisses
        self.amplicon_len = amplicon
        logging.debug("amplicon: "+str(self.amplicon_len))
        
        
        self.mismF_loc_raw, self.mismF_loc, self.mismR_loc_raw, self.mismR_loc = self._get_missmatch_location()
        self.mismF_type, self.mismR_type = self._get_missmatch_type()
        
        self.mismF_base, self.mismR_base = self._get_missmatch_base_type()
        self.amplicon = self.get_amplicon()
        
        self.add_positive_2_stats()
        
        return
    
    def complete_from_csv(self, gen, primer_pair, real_fpos, real_rpos, fmisses, rmisses, amplicon_len):
        #TODO instead of making a complete output file, calculate only the paramaters needed by the user
        self.gen = gen
        self.fastaid = gen.id #TODO patch
        self.primer_pair = primer_pair
        self.primerPair = primer_pair.id #TODO patch
        self.primerF = primer_pair.f.id
        self.primerR = primer_pair.f.id
        
        self.real_fpos = int(real_fpos)-1
        self.F_pos = int(real_fpos)-1
        self.real_rpos = int(real_rpos)-1
        self.R_pos = int(real_rpos)-1
        self.mismF = fmisses
        self.mismR = rmisses
        self.amplicon_len = amplicon_len
        
        if(self.real_fpos == None):
            raise ValueError("Error: Froward's position NULL")
            
        if(self.mismF == None):
            self.mismF = self.get_missmatches("f")
        else:
            self.mismF = fmisses
            
        if(self.mismR == None):
            self.mismR = self.get_missmatches("r")
        else:
            self.mismR = rmisses
            
        if(self.amplicon_len == None):
            if(self.real_rpos == None):
                raise ValueError("Error Couldn't determine Reverse's position")
            else:
                self.amplicon_len = self.real_rpos - self.real_fpos
        else:
            self.amplicon_len = amplicon
        self.amplicon = self.get_amplicon()
        
        self.mismF_loc_raw, self.mismF_loc, self.mismR_loc_raw, self.mismR_loc = self._get_missmatch_location()
        self.mismF_type, self.mismR_type = self._get_missmatch_type()        
        self.mismF_base, self.mismR_base = self._get_missmatch_base_type()
        
        self.add_positive_2_stats()
            
        return
    def get_Nend(self, primerPair,fastaid,primerF,primerR,mismFT,mismRT,amplicon_len,F_pos,mismFT_loc,mismFT_type,mismFT_base,R_pos,mismRT_loc,
                 mismRT_type,mismRT_base, amplicon, nend):
        """
        @brief Gets alignment and computes the alignment in Nend mode. There are more arguments than the required for convenience.
        """
        #TODO, better remove id
        
        self.primerPair = primerPair
        self.fastaid = fastaid
        self.primerF = primerF
        self.primerR = primerR
        self.F_pos = F_pos-1
        self.real_fpos = F_pos-1 #TODO patch
        self.R_pos = R_pos-1
        self.real_rpos = R_pos-1 #TODO patch
        self.amplicon_len = amplicon_len
        self.amplicon = amplicon
        
        
        self.mismF_loc, self.mismR_loc = self._get_nend_loc(mismFT_loc, mismRT_loc, nend)
        self.mismF = len(self.mismF_loc)
        self.mismR = len(self.mismR_loc)
        self.mismF_type = mismFT_type[-(self.mismF+1): -1]
        self.mismR_type = mismRT_type[0: self.mismR]
        self.mismF_base = mismFT_base[-(self.mismF+1): -1]
        self.mismR_base = mismRT_base[0: self.mismR]
        
        self.add_positive_2_stats()
            
        return self.get_csv()
    
    def add_positive_2_stats(self):
        pp_stats = self.pp_stats
        primer = self.primerPair
        m = self.mismF
        try:
            self.pp_stats.loc[self.primerPair, self.mismF+self.mismR] += 1
            self.amplicon_stats[self.primerPair].append(self.amplicon_len)
        except:
            self.pp_stats.loc[self.primerPair] = 0
            self.pp_stats.loc[self.primerPair, self.mismF+self.mismR] = 1
            self.amplicon_stats[self.primerPair] = [self.amplicon_len] 
        finally:
            self.pp_stats.loc[self.primerPair, "Total"] += 1
        return
        
    
    def add_negative_2_stats(self, primerPair):
        try:
            self.pp_stats.loc[primerPair, "No"] += 1
        except:
            self.pp_stats.loc[primerPair] = 0
            self.pp_stats.loc[primerPair, "No"] = 1
            self.amplicon_stats[primerPair] = [] #init amplicon table
        finally:
            self.pp_stats.loc[primerPair, "Total"] += 1
        return
    
    def get_stats(self):
        
        columns = self.miss_range
        pp_stats = self.pp_stats[columns] #pop "No", "Total", "ampliconLen" columns to compute the stats
        
        
        tmp = pp_stats.multiply(columns)
        total = pp_stats.sum(axis=1)
        #total = total.loc[total!=0] #remove primers with 0 matches in cooked stats
        index = pp_stats.loc[total.index].index.values
        cooked_stats = pd.DataFrame(index=index, columns=["min", "max", "mean", "median", "n_samples", 
                                                          "amplicon_min", "amplicon_max", "amplicon_median"])
        #find mean
        cooked_stats["mean"] = tmp.sum(axis=1).div(total)
        
        cooked_stats["n_samples"] = total
        
        #TODO this simple code could be improved
        #find minimum
        for i in index:
            for j in columns:
                if(pp_stats.loc[i, j]):
                    cooked_stats.loc[i, "min"] = j 
                    break;
                    
        #find maximum
        for i in index:
            for j in columns:
                if(pp_stats.loc[i, j]):
                    cooked_stats.loc[i, "max"] = j 
        #find median
        for i in index:
            a = 0
            for j in columns:
                a += pp_stats.loc[i,j]
                if(a>=total.loc[i]/2):
                    cooked_stats.loc[i, "median"] = j
                    break;
                    
        #amplicon stats
        for pp in self.amplicon_stats:
            if(len(self.amplicon_stats[pp])>0):
                self.amplicon_stats[pp].sort()
                cooked_stats.at[pp, "amplicon_min"] = self.amplicon_stats[pp][0]
                cooked_stats.at[pp, "amplicon_max"] = self.amplicon_stats[pp][-1]
                pos = int(len(self.amplicon_stats[pp])/2)
                cooked_stats.at[pp, "amplicon_median"] = self.amplicon_stats[pp][pos]
            
        cooked_stats = cooked_stats.fillna('-') #replace NaN values with - (NaN values appear when computing stats of primers with no mathces)
        
        return self.pp_stats, cooked_stats
    
    def _get_missmatch_location(self):
        """
        @Brief Returns array with the location of missmatches (on the primer)
        """
        fm_loc = []
        fm_loc_output = [] #fm_loc but inversed and starting at one, formated for the output file
        rm_loc = []
        rm_loc_output = [] #rm_loc, but starting at one, formated for the output file
        
        flen = self.primer_pair.flen
        for i in range(flen):
            if(self.F_pos+i<0 or MATCH_TABLE.loc[self.primer_pair.f.seq[i], self.gen.seq[self.F_pos+i]]!=1):
                fm_loc.append(i)
                fm_loc_output.append(flen-i)
        
        rlen = self.primer_pair.rlen
        leng = len(self.gen)
        for i in range(rlen):
            if(self.R_pos+i>=leng or MATCH_TABLE.loc[self.primer_pair.r.seq[i], self.gen.seq[self.R_pos+i]]!=1):
                    rm_loc.append(i)
                    rm_loc_output.append(i+1)
                
        return fm_loc, fm_loc_output, rm_loc, rm_loc_output
    
    def _get_missmatch_type(self):
        fm_type = []
        rm_type = []
        #TODO ask format of primers, in order to know if the gen should be compared against the compelement
        for m in self.mismF_loc_raw:
            if(self.F_pos+m>0):
                fm_type.append(self.gen.seq[self.F_pos+m]+self.primer_pair.fcomplement[m])
            else:
                fm_type.append("Z"+self.primer_pair.fcomplement[m])
          
        leng = len(self.gen)
        for m in self.mismR_loc_raw:
            if(self.R_pos+m<leng):
                rm_type.append(self.gen.seq[self.R_pos+m]+self.primer_pair.rcomplement[m])
            else:
                fm_type.append("Z"+self.primer_pair.rcomplement[m])
            
        return fm_type, rm_type
    
    def _get_missmatch_base_type(self):
        fm_base_type = []
        fprimer_complement = self.primer_pair.fcomplement
        
        for i in range(self.mismF):
            gen_nucleotide = self.mismF_type[i][0]
            f_nucleotide =fprimer_complement[self.mismF_loc_raw[i]]
            
            if(gen_nucleotide in self.base_type and f_nucleotide in self.base_type):
                gen_nucleotide_base_type = self.base_type[gen_nucleotide]
                f_nucleotide_base_type = self.base_type[f_nucleotide]
                fm_base_type.append(gen_nucleotide_base_type+"-"+f_nucleotide_base_type)
            else:
                fm_base_type.append(self.base_type["Other"])
        
        rm_base_type = []
        rprimer_complement = self.primer_pair.rcomplement
        
        for i in range(self.mismR):
            gen_nucleotide = self.mismR_type[i][0]
            r_nucleotide = rprimer_complement[self.mismR_loc_raw[i]]
            
            if(gen_nucleotide in self.base_type and r_nucleotide in self.base_type):
                gen_nucleotide_base_type = self.base_type[gen_nucleotide]
                r_nucleotide_base_type = self.base_type[r_nucleotide]
                rm_base_type.append(gen_nucleotide_base_type+"-"+r_nucleotide_base_type)
            else:
                rm_base_type.append(self.base_type["Other"])
                
        return  fm_base_type, rm_base_type
    
    def _get_nend_loc(self, mismFT_loc, mismRT_loc, nend):
        """
        @brief Returns mismataches locations that are in a Nend position
        @param mismFT_loc location of forward mismatches, reversed and starting at 1
        @param mismRT_loc location of reverse mismatches, starting at 1
        """
        mismF_loc = []
        for i in mismFT_loc:
            if(i<=nend):
                mismF_loc.append(i)
        
        mismR_loc = []
        for i in mismRT_loc:
            if(i>nend):
                break
            mismR_loc.append(i)
        return mismF_loc, mismR_loc    
    
    def get_amplicon(self):
        amplicon = self.gen[self.F_pos:self.F_pos+self.amplicon_len]
        amplicon = "".join(amplicon)
        return amplicon
    
    def get_csv(self):
        info= [self.primerPair, self.fastaid, self.primerF, self.primerR, self.mismF, self.mismR, 
               self.amplicon_len, self.real_fpos+1, self.mismF_loc, self.mismF_type, self.mismF_base, self.real_rpos+1, self.mismR_loc, 
               self.mismR_type, self.mismR_base, self.amplicon]
        """
        if(self.Nend_misses):
            info.extend([self.mismF_Nend, self.mismR_Nend])
        """
        return info

def get_missmatch_column_name(header, primer="f"):
    if(primer=="f"):
        f_reg = re.compile("^mismF((N\d+)|T)$")
        return list(filter(f_reg.match, header))[0]
    elif(primer=="r"):
        f_reg = re.compile("^mismR((N\d+)|T)$")
        return list(filter(f_reg.match, header))[0]
    
#LOGGER
file_handler = None
console_handler = None

def init_logger():
    global console_handler
    global file_handler
    
    #root_handler = logging.getLogger()
    #logging.basicConfig(filename=os.path.join(os.getcwd(),"log.txt"), filemode='w', level=logging.INFO)
    
    log = logging.getLogger()
    log.setLevel(logging.DEBUG)
    
    for hdlr in log.handlers[:]:  # remove all old handlers
        log.removeHandler(hdlr)

    file_handler = logging.FileHandler(os.path.join(os.getcwd(),"log.txt"), 'a')
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))
    
    log.addHandler(console_handler)
    log.addHandler(file_handler)
    return
        
def set_verbosity(verbosity):
    #verbosity = 2
    if verbosity == True:
        file_handler.setLevel(logging.INFO)
        console_handler.setLevel(logging.WARNING)
    elif verbosity == 2: #True!=2, ugly but as long as it works...
        file_handler.setLevel(logging.DEBUG)
        console_handler.setLevel(logging.DEBUG)
    else:
        file_handler.setLevel(logging.WARNING)
        console_handler.setLevel(logging.ERROR)
    return

def close_logger():
    try:
        log = logging.getLogger()
        for hdlr in log.handlers[:]:  # remove all old handlers
            log.removeHandler(hdlr)
    except:
        logging.error("Can't close console log")
        
    logging.shutdown()