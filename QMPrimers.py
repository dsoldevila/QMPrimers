#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 15:30:29 2018

@author: david
"""

import load_data as ld
import matching as m
from common import *
from tkinter import filedialog
from tkinter import *
import os

class GUI(Frame):
    def __init__(self, parent=Frame):
        Frame.__init__(self, parent)
        
        self.label_g = Label(text="No Genome")
        self.label_g.pack()
        self.button_g = Button(text="Open", command=self.open_bio_file)
        self.button_g.pack()
        
        self.label_p = Label(text="No Primer Pair(s)")
        self.label_p.pack()
        self.button_p = Button(text="Open", command=self.open_csv_file)
        self.button_p.pack()
        
        self.button_c = Button(text="Compute", command=self.compute)
        self.button_c.pack()
        return
    
    def _open_file(self):
        return filedialog.askopenfilename(initialdir = os.getcwd(),title = "Select file",filetypes = (("all files","*.*"),))
    
    def open_bio_file(self):
        input_file = self._open_file()
        self.label_g.config(text=input_file)
        print(input_file)
        self.gen_record = ld.load_bio_file(input_file, file_format=None)
        print(self.gen_record.get("ACEA1016-14_Aphis_spiraecola_BOLD"))
        return
    
    def open_csv_file(self):
        input_file = self._open_file()
        self.label_p.config(text=input_file)
        print(input_file)
        self.primer_list = ld.load_csv_file(input_file)
        print(self.primer_list[0].f.seq)
        return
    
    def compute(self):
        pass
    

if "__main__":
    root = Tk()
    root.title("QMPrimers")
    main_window = GUI(root)
    root.mainloop()
    