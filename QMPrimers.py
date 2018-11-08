#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 15:30:29 2018

@author: david
"""

import load_data as ld
import matching as m

if "__main__":
    gen_seqs = ld.loadBioFile("species_bold_own_genbank.fasta")
    """for gen in gen_seqs:
        print(gen_seqs.get(gen).id)"""
    primer_pairs = ld.loadCSVFile("P&PP.csv")
    template_primer_missmatches = m.compute_matching(5, 5, primer_pairs, gen_seqs, 1)