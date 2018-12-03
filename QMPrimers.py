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
        self.main_frame = Frame.__init__(self, parent)
        
        """Menu"""
        self.main_menu = Menu(parent)
        parent.config(menu=self.main_menu)
        file = Menu(self.main_menu, tearoff=False)
        file.add_command(label="New")
        self.main_menu.add_cascade(label="File", menu=file)
        
        """Select Files Frame"""
        self.file_frame = Frame(self.main_frame)
        self.file_frame.pack()
        
        self.entry_g = Entry(self.file_frame)
        self.entry_g.pack(side=LEFT, expand=YES, fill=X)
        self.entry_g.insert(0, "<No Genome>")
        self.entry_g.bind("<Return>", (lambda event: self.open_bio_file(self.entry_g.get())))
        
        self.button_g = Button(self.file_frame, text="Open", command=(lambda: self.open_bio_file(self._open_file())))
        self.button_g.pack(side=LEFT, expand=YES, fill=X)
        
        self.entry_p = Entry(self.file_frame)
        self.entry_p.pack(side=LEFT, expand=YES, fill=X)
        self.entry_p.insert(0, "<No Primer Pairs>")
        self.entry_p.bind("<Return>", (lambda event: self.open_csv_file(self.entry_p.get())))
        
        self.button_p = Button(self.file_frame, text="Open", command=(lambda: self.open_csv_file(self._open_file())))
        self.button_p.pack(side=LEFT, expand=YES, fill=X)
        
        """Check List Frame"""
        self.check_frame = Frame(self.main_frame)
        self.check_frame.pack(expand=YES, fill=X)
        self.check_list = []
        self.check_list_names = ["R is reverse complement", "option2", "option3"]
        for i in range(len(self.check_list_names)):
            check = Checkbutton(self.check_frame,text=self.check_list_names[i])
            check.pack(side=LEFT, expand=YES, fill=X)
            self.check_list.append(check)
        
        """Compute"""
        self.button_c = Button(text="Compute", command=self.compute)
        self.button_c.pack(expand=YES, fill=BOTH)
        
        """Other"""
        self.current_directory = os.getcwd()
        return
        
    def _open_file(self):
        file_name = filedialog.askopenfilename(initialdir=self.current_directory, title = "Select file")
        if(file_name): self.current_directory = os.path.dirname(file_name)
        return file_name
    
    def open_bio_file(self, input_file):
        if(input_file): #is not None
            self.entry_g.delete(0, END)
            self.entry_g.insert(0, input_file)
            print(input_file)
            self.gen_record = ld.load_bio_file(input_file, file_format=None)
            print(self.gen_record.get("ACEA1016-14_Aphis_spiraecola_BOLD")+"\n")
        return
    
    def open_csv_file(self, input_file):
        if(input_file): #is not None
            self.entry_p.delete(0, END)
            self.entry_p.insert(0, input_file)
            print(input_file)
            self.primer_pairs = ld.load_csv_file(input_file)
            print(self.primer_pairs[0].f.seq+"\n")
        return
    
    def compute(self):
        result = m.compute_template_missmatches(10, 10, self.primer_pairs, self.gen_record)
        print(result[0])
        return
    

if (__name__=="__main__"):
    root = Tk()
    root.title("QMPrimers")
    main_window = GUI(root)
    root.mainloop()
    