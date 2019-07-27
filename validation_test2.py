#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 12:58:15 2019

@author: david
"""

import QMPrimers as qm
import os

def matching_test():
    pass

if(__name__=="__main__"):
    rootd = os.getcwd()
    ind = os.path.join(rootd, "test_input")
    outd = os.path.join(rootd, "test_output")
    vald = os.path.join(rootd, "test_validated")
    
    argv = ["--sim", "--help"]
    qm.main(argv)