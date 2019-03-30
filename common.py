#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 15:34:23 2018

@author: David Soldevila
"""
import numpy as np
import pandas as pd

IUPAC_AMBIGUOUS_DNA = tuple("ACGTWSMKRYBDHVNIZ")

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

    def __init__(self, gen, primer_pair, fpos, real_fpos, rpos, real_rpos, fmisses, rmisses, amplicon, MATCH_TABLE):
    
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
        
        self.fm_loc, self.rm_loc = self._get_missmatch_location(MATCH_TABLE)
        self.fm_type, self.rm_type = self._get_missmatch_type()
        
        self.fm_base, self.rm_base = self._get_missmatch_base_type()
        
        return
    
    def _get_missmatch_location(self, MATCH_TABLE):
        """
        @Brief Returns array with the location of missmatches (on the primer)
        """
        fm_loc = []
        rm_loc = []
        
        for i in range(self.primer_pair.flen):
            if(MATCH_TABLE.loc[self.primer_pair.f.seq[i], self.gen[self.fpos+i]]!=1):
                fm_loc.append(i)
                
        for i in range(self.primer_pair.rlen):
            if(MATCH_TABLE.loc[self.primer_pair.r.seq[i], self.gen.seq[self.rpos+i]]!=1):
                    rm_loc.append(i)
                
        return fm_loc, rm_loc
    
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
    
    def __str__(self):        
        info = ("PRIME PAIR "+str(self.primer_pair.id)+"\n"+
              "Forward's at: "+str(self.real_fpos)+" with "+str(self.fm)+" misses "+ str(self.fm_loc)+" "+str(self.fm_type)+"\n"+
              "Reverse's at: "+str(self.real_rpos)+" with "+str(self.rm)+" misses "+ str(self.rm_loc)+" "+str(self.rm_type)+"\n"+
              "Amplicon: "+str(self.amplicon)+"\n")
        return info
    
    def get_csv(self):
        info= [self.primer_pair.id, self.gen.id, self.primer_pair.f.id, self.primer_pair.r.id, self.fm, self.rm, 
               self.amplicon, self.real_fpos, str(self.fm_loc), str(self.fm_type), str(self.fm_base), self.real_rpos, str(self.rm_loc), str(self.rm_type), str(self.rm_base)]
        return info