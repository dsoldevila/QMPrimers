#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 15:34:23 2018

@author: david
"""

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

class Matching:
    def __init__(self, gen, primer_pair, range_f, range_r, misses_f, misses_f_loc, misses_r, misses_r_loc, amplicon):
        self.gen = gen
        self.primer_pair = primer_pair
        self.range_f = range_f
        self.range_r = range_r
        self.misses_f = misses_f
        self.misses_f_loc = self._misses2str(misses_f_loc)
        self.misses_r = misses_r
        self.misses_r_loc = self._misses2str(misses_r_loc)
        self.amplicon = amplicon
        
    def _misses2str(self, misses_loc):
        array = " " if misses_loc[0] else "|"
        for i in range(1, len(misses_loc)):
            char = " " if misses_loc[i] else "|"
            array +=char
        return array
        
    def __str__(self):
        info = ("PRIME PAIR "+str(self.primer_pair.id)+"\n"+
              "Forward's missmatches: "+str(self.misses_f)+"\n"+
              "Backward's missmatches: "+str(self.misses_r)+"\n"+
              "Pair's amplicon: "+str(self.amplicon)+"\n"+
              "Forward's match "+ str(self.range_f)+"\n"+
              " "+str(self.gen.seq[self.range_f[0]:self.range_f[1]])+"\n"+
              " "+self.misses_f_loc+"\n"+
              " "+str(self.primer_pair.f.seq) +"\n"+
              "Backwards's match "+ str(self.range_r)+"\n"+
              " "+str(self.gen.seq[self.range_r[0]:self.range_r[1]])+"\n"+
              " "+self.misses_r_loc+"\n"+
              " "+str(self.primer_pair.r.seq) +"\n")
        return info