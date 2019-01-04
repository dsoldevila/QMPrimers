#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 15:34:23 2018

@author: David Soldevila
"""
import numpy as np
import pandas as pd

class PrimerPair:
    def __init__(self, pair_id, fprimer, rprimer, min_amplicon, max_amplicon):
        self.id = pair_id
        self.f = fprimer
        self.flen = len(fprimer.seq)
        self.r = rprimer
        self.rlen = len(rprimer.seq)
        self.min_amplicon = min_amplicon
        self.max_amplicon = max_amplicon
        return
    
    def __str__(self):
        pass
    
class Alignment:
    """
    Alignment info between a genomic sequence and a primer pair
    """

    def __init__(self, gen, primer_pair, fpos, rpos, fmisses, rmisses, amplicon, MATCH_TABLE):
    
        """
        self.gen = genomic sequence
        self.primer_pair = primer pair used for matching, instance of PrimerPair class
        self.fpos = starting position of the forward primer in the genomic sequence, starting at 0
        self.rpos = starting position of the reverse primer in the genomic sequence, starting at the end of genomic sequence
        self.fm = number of missmatches on the forward primer
        self.fm_loc = array of the missmatch locations of the forward primer
        self.rm = number of missmatches on the reverse primer
        self.rm_loc = array of the missmatch locations of the reverse primer
        self.amplicon = amplicon of the matching, number between the primer pair max and min amplicon.
        """
        self.gen = gen
        self.primer_pair = primer_pair
        self.fpos = int(fpos) #it seems Biopython seqrecord does not support numpy.int32
        self.rpos = int(rpos)
        self.fm = fmisses
        self.rm = rmisses
        self.amplicon = amplicon
        
        self.fm_loc, self.rm_loc = self._get_missmatch_location(MATCH_TABLE)
        self.fm_type, self.rm_type = self._get_missmatch_type()
        
        return
    
    def _get_missmatch_location(self, MATCH_TABLE):
        fm_loc = []
        rm_loc = []
        
        for i in range(self.primer_pair.flen):
            if(MATCH_TABLE.loc[self.gen[self.fpos+i],self.primer_pair.f.seq[i]]!=1):
                fm_loc.append(i)
                
        for i in range(self.primer_pair.rlen):
            if(MATCH_TABLE.loc[self.gen.seq[self.rpos+i], self.primer_pair.r.seq[i]]!=1):
                    rm_loc.append(i)
                
        return fm_loc, rm_loc
    
    def _get_missmatch_type(self):
        fm_type = []
        rm_type = []
    
        for m in self.fm_loc:
            fm_type.append(self.gen.seq[self.fpos+m]+self.primer_pair.f.seq[m])
          
        for m in self.rm_loc:
            rm_type.append(self.gen.seq[self.rpos+m]+self.primer_pair.r.seq[m])
            
            
        return fm_type, rm_type
    
    def __str__(self):
        info = ("PRIME PAIR "+str(self.primer_pair.id)+"\n"+
              "Forward's missmatches: "+str(self.fm)+"\n"+
              "Reverse's missmatches: "+str(self.rm)+"\n"+
              "Pair's amplicon: "+str(self.amplicon)+"\n"+
              "Forward's match at "+ str(self.fpos)+"\n"+
              " "+str(self.gen.seq[self.fpos:self.fpos+self.primer_pair.flen])+"\n"+
              " "+str(self.primer_pair.f.seq) +"\n"+
              "Forward's missmatch Loc: "+ str(self.fm_loc)+"\n"+
              "Forward's missmatch Type: "+ str(self.fm_type)+"\n"+
              "Reverse's match at "+ str(self.rpos)+"\n"+
              " "+str(self.gen.seq[self.rpos:self.rpos+self.primer_pair.rlen])+"\n"+
              " "+str(self.primer_pair.r.seq) +"\n"+
              "Reverse's missmatch Loc: "+ str(self.rm_loc)+"\n"+
              "Reverse's missmatch Type: "+ str(self.rm_type)+"\n")
        return info
    
class AlignmentList:
    """
    It's plausible that a primer pair aligns in more than one place in the sequence, a list is needed
    """
    def __init__(self, primer_pair, gen):
        self.primer_pair = primer_pair
        self.gen = gen
        self._alignment_list = []
        return
    
    def append(self, match):
        if(type(match)==list):
            self._alignment_list.extend(match)
        else:
            self._alignment_list.append(match)
        return
    
    def get_list(self):
        return self._alignment_list
    
    def get_next_alignment(self, pair_id):
        pass
    
    def __str__(self):
        info = ""
        for al in self._alignment_list:
            info += str(al)
        return info
    
class GenMatching:
    def __init__(self, gen):
        self.gen = gen
        self._matching_list = []
        return
    
    def append(self, alignment_list):
        self._matching_list.append(alignment_list)
    
    def __str__(self):
        info = info = "------------\nFOR: "+self.gen.id+"\n-------------\n"
        for m in self._matching_list:
            info += str(m)
        return info
