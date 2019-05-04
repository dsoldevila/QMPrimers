#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 15:34:23 2018

@author: David Soldevila
"""
import numpy as np
import pandas as pd

IUPAC_AMBIGUOUS_DNA = tuple("ACGTWSMKRYBDHVNIZ")
TEMPLATE_HEADER = ["primerPair","fastaid","primerF","primerR","mismFT","mismRT","amplicon", "F_pos", "mismFT_loc", "mismFT_type", 
                                     "mismFT_base", "R_pos", "mismRT_loc", "mismRT_type", "mismRT_base"]

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
        columns.append("No")
        self.pp_stats = pd.DataFrame(columns=columns, dtype='uint8')
        return
    
    def get(self, gen, primer_pair, fpos, real_fpos, rpos, real_rpos, fmisses, rmisses, amplicon, Nend_misses):
    
        """
        self.gen = genomic sequence
        self.primer_pair = primer pair used for matching, instance of PrimerPair class
        self.fpos = starting position of the forward primer in the genomic sequence, starting at 0
        self.real_fpos --> gen:    ZZAGTAC...     real_fpos = -2, the primer is hanging
                           primer: AGAGT          fpos = 0
        self.rpos = starting position of the reverse primer in the genomic sequence, starting at the end of genomic sequence
        self.real_rpos = reverse's position depends on forward's position
        self.fm = number of missmatches on the forward primer
        self.fm_loc = array of the missmatch locations of the forward primer
        self.rm = number of missmatches on the reverse primer
        self.rm_loc = array of the missmatch locations of the reverse primer
        self.amplicon = amplicon of the matching, number between the primer pair max and min amplicon.
        """
        self.gen = gen
        self.primer_pair = primer_pair
        self.fpos = int(fpos) #it seems Biopython seqrecord does not support numpy.int32
        self.real_fpos = int(real_fpos)
        self.rpos = int(rpos)
        self.real_rpos = int(real_rpos)
        self.fm = fmisses
        self.rm = rmisses
        self.amplicon = amplicon
        
        
        self.fm_loc, self.fm_loc_output, self.rm_loc, self.rm_loc_output = self._get_missmatch_location()
        self.fm_type, self.rm_type = self._get_missmatch_type()
        
        self.fm_base, self.rm_base = self._get_missmatch_base_type()
        
        self.Nend_misses = Nend_misses
        if(Nend_misses):
            self.fm_Nend, self.rm_Nend = self._get_Nend_missmatches(Nend_misses)
            
        try:
            self.pp_stats.loc[primer_pair.id, fmisses+rmisses] += 1
        except:
            self.pp_stats.loc[primer_pair.id] = 0
            self.pp_stats.loc[primer_pair.id, fmisses+rmisses] = 1
        
        return
    
    def complete_from_csv(self, gen, primer_pair, real_fpos, real_rpos, fmisses, rmisses, amplicon, Nend_misses=None):
        #TODO instead of making a complete output file, calculate only the paramaters needed by the user
        self.gen = gen
        self.primer_pair = primer_pair
        self.real_fpos = int(real_fpos)
        self.fpos = int(real_fpos)
        self.real_rpos = int(real_rpos)
        self.rpos = int(real_rpos)
        self.fm = fmisses
        self.rm = rmisses
        self.amplicon = amplicon
        
        if(self.real_fpos == None):
            raise ValueError("Error: Froward's position NULL")
            
        if(self.fm == None):
            self.fm = self.get_missmatches("f")
        else:
            self.fm = fmisses
            
        if(self.rm == None):
            self.rm = self.get_missmatches("r")
        else:
            self.rm = rmisses
            
        if(self.amplicon == None):
            if(self.real_rpos == None):
                raise ValueError("Error Couldn't determine Reverse's position")
            else:
                self.amplicon = self.real_rpos - self.real_fpos
        else:
            self.amplicon = amplicon
        
        self.fm_loc, self.fm_loc_output, self.rm_loc, self.rm_loc_output = self._get_missmatch_location()
        self.fm_type, self.rm_type = self._get_missmatch_type()        
        self.fm_base, self.rm_base = self._get_missmatch_base_type()
        
        self.Nend_misses = Nend_misses
        if(self.Nend_misses):
            self.fm_Nend, self.rm_Nend = self._get_Nend_missmatches(self.Nend_misses)
        
        try:
            self.pp_stats.loc[primer_pair.id, fmisses+rmisses] += 1
        except:
            self.pp_stats.loc[primer_pair.id] = 0
            self.pp_stats.loc[primer_pair.id, fmisses+rmisses] = 1
            
        return
    
    def add_negative_2_stats(self, primer_pair_id):
        try:
            self.pp_stats.loc[primer_pair_id, "No"] += 1
        except:
            self.pp_stats.loc[primer_pair_id] = 0
            self.pp_stats.loc[primer_pair_id, "No"] = 1
        return
    
    def get_stats(self):
        
        columns = self.pp_stats.columns.values[:-1]
        pp_stats = self.pp_stats[columns] #pop "No" columns to compute the stats
        
        
        tmp = pp_stats.multiply(columns)
        total = pp_stats.sum(axis=1)
        total = total.loc[total!=0]
        index = pp_stats.loc[total.index].index.values
        cooked_stats = pd.DataFrame(index=index, columns=["min", "max", "mean", "median", "n_samples"])
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
            if(MATCH_TABLE.loc[self.primer_pair.f.seq[i], self.gen.seq[self.fpos+i]]!=1):
                fm_loc.append(i)
                fm_loc_output.append(flen-i+1)
        
        rlen = self.primer_pair.rlen
        for i in range(rlen):
            if(MATCH_TABLE.loc[self.primer_pair.r.seq[i], self.gen.seq[self.rpos+i]]!=1):
                    rm_loc.append(i)
                    rm_loc_output.append(i+1)
                
        return fm_loc, fm_loc_output, rm_loc, rm_loc_output
    
    def _get_missmatch_type(self):
        fm_type = []
        rm_type = []
        #TODO ask format of primers, in order to know if the gen should be compared against the compelement
        for m in self.fm_loc:
            fm_type.append(self.gen.seq[self.fpos+m]+self.primer_pair.fcomplement[m])
          
        for m in self.rm_loc:
            rm_type.append(self.gen.seq[self.rpos+m]+self.primer_pair.rcomplement[m])
            
        return fm_type, rm_type
    
    def _get_missmatch_base_type(self):
        fm_base_type = []
        fprimer_complement = self.primer_pair.fcomplement
        
        for i in range(self.fm):
            gen_nucleotide = self.fm_type[i][0]
            f_nucleotide =fprimer_complement[self.fm_loc[i]]
            
            if(gen_nucleotide in self.base_type and f_nucleotide in self.base_type):
                gen_nucleotide_base_type = self.base_type[gen_nucleotide]
                f_nucleotide_base_type = self.base_type[f_nucleotide]
                fm_base_type.append(gen_nucleotide_base_type+"-"+f_nucleotide_base_type)
            else:
                fm_base_type.append(self.base_type["Other"])
        
        rm_base_type = []
        rprimer_complement = self.primer_pair.rcomplement
        
        for i in range(self.rm):
            gen_nucleotide = self.rm_type[i]
            r_nucleotide = rprimer_complement[self.rm_loc[i]]
            
            if(gen_nucleotide in self.base_type and r_nucleotide in self.base_type):
                gen_nucleotide_base_type = self.base_type[gen_nucleotide]
                r_nucleotide_base_type = self.base_type[r_nucleotide]
                rm_base_type.append(gen_nucleotide_base_type+"-"+r_nucleotide_base_type)
            else:
                rm_base_type.append(self.base_type["Other"])
                
        return  fm_base_type, rm_base_type
            
    def _get_Nend_missmatches(self, Nend_misses):
        rm_Nend = 0
        for i in self.rm_loc:
            if(i >= Nend_misses):
                break;
            rm_Nend +=1
            
        fm_Nend = 0
        Nend_misses = self.primer_pair.flen - Nend_misses
        for i in range(1, len(self.fm_loc)+1):
            if(self.fm_loc[-i] < Nend_misses):
                break;
            fm_Nend +=1
        return fm_Nend, rm_Nend
    
    
    def __str__(self):        
        info = ("PRIME PAIR "+str(self.primer_pair.id)+"\n"+
              "Forward's at: "+str(self.real_fpos)+" with "+str(self.fm)+" misses "+ str(self.fm_loc)+" "+str(self.fm_type)+"\n"+
              "Reverse's at: "+str(self.real_rpos)+" with "+str(self.rm)+" misses "+ str(self.rm_loc)+" "+str(self.rm_type)+"\n"+
              "Amplicon: "+str(self.amplicon)+"\n")
        return info
    
    def get_csv(self):
        info= [self.primer_pair.id, self.gen.id, self.primer_pair.f.id, self.primer_pair.r.id, self.fm, self.rm, 
               self.amplicon, self.real_fpos, self.fm_loc_output, self.fm_type, self.fm_base, self.real_rpos, self.rm_loc_output, 
               self.rm_type, self.rm_base]
        if(self.Nend_misses):
            info.extend([self.fm_Nend, self.rm_Nend])
        return info

def get_Nend_missmatches(Nend_misses, rm_loc_output, flen, fm_loc_output):
    """
    @brief: Modified version of Alignment._get_Nend_missmatches. It's needded because missmatch locations in the output_file are formated 
    differently from the ones used by the program internally. Check fm_loc and fm_loc_output.
    """
    rm_Nend = 0
    for i in rm_loc_output:
        if(i > Nend_misses):
            break;
        rm_Nend +=1
        
    fm_Nend = 0
    for i in fm_loc_output:
        if(i > Nend_misses):
            break;
        fm_Nend +=1
            
    return fm_Nend, rm_Nend