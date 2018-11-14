#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 11:29:49 2018

@author: david
"""

import _thread

def child(tid):
    print("Hello from thread", tid)
    return

def parent():
    i = 0
    while True:
        i+=1
        _thread.start_new_thread(child, (i,))
        if input() == 'q': break
    return

if (__name__=="__main__"):
    parent()