#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 10:06:42 2019

@author: david
"""
from interface import *
from common import *
from simulation import *
from tkinter import filedialog
from tkinter import *
from tkinter.font import Font
import os
import sys
import _thread
import pandas as pd

class GUI_compute():
    def __init__(self, parent, parameters, output_info, gui_simulate):
        
        self.main_frame = Frame(parent)

        #myFont = Font(family="Noto Sans", size=12, weight="bold")
        
        """Other"""
        self.parameters = parameters
        self.current_directory = os.getcwd()
        self.gui_simulate = gui_simulate
        
        """Frame containing Files and Other Parameters frames"""
        self.first_row_frame = Frame(self.main_frame)
        self.first_row_frame.pack(expand=YES, fill=BOTH)
        
        """Select Files Frame"""
        self.file_frame = Frame(self.first_row_frame)
        self.file_frame.pack(side=LEFT, expand=YES, fill=X)
        Label(self.file_frame, text="Files").pack()
        
        
        self.entries = {}
        self.button = {}
        for name in self.parameters.index:
            if(self.parameters.loc[name, "type"] == "entry"):
                block = Frame(self.file_frame)
                block.pack(expand=YES, fill=BOTH)
                self.entries[name] = Entry(block)
                self.entries[name].pack(side=LEFT, expand=YES, fill=BOTH)
                self.button[name] = Button(block, text="Open")
                self.button[name].pack(side=LEFT, expand=NO, fill=X)
                
        self.entries["gen"].insert(0, "<No Genome>")
        self.entries["gen"].bind("<Return>", (lambda event: self.update_bio_files(self.entries["gen"].get())))
        self.button["gen"].config(command=(lambda: self.update_bio_files(self._open_files())))
        
        self.entries["primer_pairs"].insert(0, "<No Primer Pairs>")
        self.entries["primer_pairs"].bind("<Return>", (lambda event: self.update_csv_file(self.entries["primer_pairs"].get())))
        self.button["primer_pairs"].config(command=(lambda: self.update_primer_file(self._open_files())))
        
        self.entries["output_file"].insert(0, self.parameters.loc["output_file", "value"])
        self.entries["output_file"].bind("<Return>", (lambda event: self.set_output_file(self.entries["output_file"].get())))
        self.button["output_file"].config(text="Set", command=(lambda: self.set_output_file(self.entries["output_file"].get())))
        
        self.entries["csv_template"].insert(0, "<No Precomputed Template>")
        self.entries["csv_template"].bind("<Return>", (lambda event: self.set_template_file(self.entries["csv_template"].get())))
        self.button["csv_template"].config(command=(lambda: self.set_template_file(self._open_files())))
            
            
        """Other parameters"""
        self.other_frame = Frame(self.first_row_frame)
        self.other_frame.pack(side=LEFT, expand=YES, fill=X)
        Label(self.other_frame, text="Parameters").pack()   
        self.other_param = {}
        for name in self.parameters.index:
            if(self.parameters.loc[name, "type"] == "param"):
                if(isinstance(self.parameters.loc[name, "value"], bool)):
                    other = BooleanVar()
                    Checkbutton(self.other_frame, variable=other, text=name).pack(expand=YES)
                    other.set(self.parameters.loc[name, "value"])
                    self.other_param[name] = other
                elif(isinstance(self.parameters.loc[name, "value"], int)):
                    block = Frame(self.other_frame)
                    block.pack(expand=YES)
                    Label(block, text=name).pack(side=RIGHT)
                    other = StringVar()
                    Entry(block, textvariable=other, width=2).pack(side=LEFT)
                    other.set(self.parameters.loc[name, "value"])
                    self.other_param[name] = other
                    
        
        """Output info"""
        self.output_frame = Frame(self.main_frame)
        self.output_frame.pack(expand=YES, fill=X)
        Label(self.output_frame, text="Output Info").grid(row=0, column=0)
        
        self.output_info = output_info
        c=0
        max_row_len=len(self.output_info)/2
        for key in self.output_info:
            var = BooleanVar()
            Checkbutton(self.output_frame, variable=var, text=key).grid(row=int(1+c/max_row_len), column=int(c%max_row_len), sticky=W)
            var.set(self.output_info[key])
            self.output_info[key] = var
            c = c+1
        for name in parameters.loc[parameters["type"]=="info"].index.values:
            block = Frame(self.output_frame)
            other = IntVar()
            Entry(block, textvariable=other, width=2).pack(side=LEFT)
            other.set(self.parameters.loc[name, "value"])
            self.other_param[name] = other 
            Label(block, text=name).pack(side=RIGHT)
            block.grid(row=int(1+c/max_row_len), column=int(c%max_row_len), sticky=W)
            c= c+1
            
                    
        #elf.entries["Nend miss."].bind("<Return>", (lambda event: self.set_Nend(self.entries["Nend miss."].get())))
        #self.entries["Nend miss."].insert(0, self.parameters.loc[name, "value"])
        #TODO limit Nend
        
        self.buttons_frame = Frame(self.main_frame)
        self.buttons_frame.pack(expand=YES, fill=X)
        
        """Compute"""
        self.button_c = Button(self.buttons_frame, text="Compute", command=self.compute_in_thread)
        self.button_c.pack(side=TOP, expand=YES, fill=X)
        
        """Load Template"""
        self.button_lt = Button(self.buttons_frame, text="Load template", command=self.load_template_in_thread)
        #self.button_c.pack(expand=YES, fill=X)
        
        """Save"""
        self.button_s = Button(self.buttons_frame, text="Save", command=self.store_results_in_thread)
        self.button_s.pack(side=BOTTOM, expand=YES, fill=X)
        
        """Other"""
        self.previous_Nend = self.parameters.loc["Nend miss.", "value"]
        
        return
    def pack(self):
        self.main_frame.pack(expand=YES, fill=BOTH)
        return
    
    def pack_forget(self):
        self.main_frame.pack_forget()
        return
        
    def _open_files(self):
        """
        Open file from the filedialog
        @return Tuple of strings
        """
        file_names = filedialog.askopenfilenames(initialdir=self.current_directory, title = "Select file")
        if(file_names): self.current_directory = os.path.dirname(file_names[0])
        if(len(file_names) == 1):
            file_names = file_names[0]
        return file_names
    
    def update_bio_files(self, input_files):
        """
        Loads input genomes 
        """
        if(input_files): #is not None
            self.entries["gen"].delete(0, END)
            self.entries["gen"].insert(0, str(input_files))
            self.parameters.loc["gen", "value"] = input_files
        return
    
    def update_primer_file(self, input_file):
        """
        Loads input primer pairs, stored in csv format
        """
        if(input_file): #is not None
            self.entries["primer_pairs"].delete(0, END)
            self.entries["primer_pairs"].insert(0, str(input_file))
            self.parameters.loc["primer_pairs", "value"] = input_file
        return
    
    def set_output_file(self, output_file):
        if(os.path.isabs(output_file)):
            self.parameters.loc["output_file", "value"] = output_file
        else:
            self.parameters.loc["output_file", "value"] = os.path.join(self.current_directory,output_file)
            self.entries["output_file"].delete(0, END)
            self.entries["output_file"].insert(0, self.parameters.loc["output_file", "value"])
        print("Output file path updated")
        return
    
    def set_Nend(self, Nend):
        self.parameters.loc["Nend miss.", "value"] = int(Nend)
        print("Nend setted to ", Nend)
        return
    
    def set_template_file(self, template_file):
        if(template_file):
            if(os.path.isabs(template_file)):
                self.parameters.loc["csv_template", "value"] = template_file
            else:
                self.parameters.loc["csv_template", "value"] = os.path.join(self.current_directory, template_file)
                
            self.entries["csv_template"].delete(0, END)
            self.entries["csv_template"].insert(0, str(template_file))
                
            self.other_frame.pack_forget()
            self.button["csv_template"].config(text="Unset", command= self.unset_template_file)
            self.button_c.pack_forget()
            self.button_lt.pack(side=TOP, expand=YES, fill=X)
        return
    
    def unset_template_file(self):
        self.button["csv_template"].config(text="Open", command=(lambda: self.set_template_file(self._open_files())))
        self.entries["csv_template"].delete(0, END)
        self.entries["csv_template"].insert(0, "<No Precomputed Template>")
        self.other_frame.pack(side=LEFT, expand=YES, fill=X)
        self.button_lt.pack_forget()
        self.button_c.pack(side=TOP, expand=YES, fill=X)
        return
    
    def compute_in_thread(self):
        _thread.start_new_thread(self.compute, ())
        return
    
    def compute(self):
        for pkey in self.other_param:
            self.parameters.loc[pkey, "value"] = self.other_param[pkey].get()
       
        self.template, self.discarded, self.gen_record, self.primer_pairs, self.raw_stats, self.cooked_stats = compute(self.parameters)
        
        print("Finished!")
        
        self.store_results()
        self.gui_simulate.set_template(self.template)

        
        return
    
    def load_template_in_thread(self):
        _thread.start_new_thread(self.load_template, ())
        return
    
    def load_template(self):
        self.template, self.discarded, self.gen_record, self.primer_pairs, self.raw_stats, self.cooked_stats = load_template(self.parameters)
        self.gui_simulate.set_template(self.template)
        return
    
    def store_results_in_thread(self):
        self.store_results()
        #_thread.start_new_thread(self.store_results, ())
        return
        
    
    def store_results(self):
        header = []
        i = 0 #TODO convert output_info to list
        for key in self.output_info:
            if(self.output_info[key].get()):
                header.append(i)
            i+=1
                   
        try:
            Nend = self.parameters.loc["Nend miss.", "value"] = self.other_param["Nend miss."].get()
        except: #Crash expected if tkinter variable is empty
            Nend = 0
        if(Nend):
            if(self.previous_Nend!=Nend):
                max_misses = int(self.parameters.loc["forward missmatches", "value"])+int(self.parameters.loc["reverse missmatches", "value"])
                self.out_template, self.out_raw_stats, self.out_cooked_stats = get_Nend_match(self.template, Nend, max_misses)
                self.previous_Nend = Nend
        else:
            self.out_template = self.template
            self.out_raw_stats = self.raw_stats
            self.out_cooked_stats = self.cooked_stats
            
        save_matching_info(self.parameters.loc["output_file", "value"], self.out_template, header, self.discarded, self.out_raw_stats, self.out_cooked_stats)
        return
    
class GUI_simulate():
    def __init__(self, parent, parameters, output_info):
        
        self.main_frame = Frame(parent)
        
        """Other"""
        self.parameters = parameters
        self.current_directory = os.getcwd()
        self.entries = {}
        self.buttons = {}
        self.template = None
        
        """Frame containing Files and Other Parameters frames"""
        self.first_row_frame = Frame(self.main_frame)
        self.first_row_frame.pack(expand=YES, fill=BOTH)
        
        """Select Files Frame"""
        self.file_frame = Frame(self.first_row_frame)
        self.file_frame.pack(side=LEFT, expand=YES, fill=X)
        Label(self.file_frame, text="Files").pack()
        
        for name in self.parameters.loc[parameters["type"]=="entry"].index.values:
            block = Frame(self.file_frame)
            block.pack(expand=YES, fill=BOTH)
            self.entries[name] = Entry(block)
            self.entries[name].pack(side=LEFT, expand=YES, fill=BOTH)
            self.buttons[name] = Button(block, text="Open")
            self.buttons[name].pack(side=LEFT, expand=NO, fill=X)
                
        self.entries["template"].insert(0, "<No Template>")
        self.entries["template"].bind("<Return>", (lambda event: self.update_template_file(self.entries["template"].get())))
        self.buttons["template"].config(command=(lambda: self.update_template_file(self._open_file())))
        
        self.entries["output_file"].insert(0, self.parameters.loc["output_file", "value"])
        self.entries["output_file"].bind("<Return>", (lambda event: self.set_output_file(self.entries["output_file"].get())))
        self.buttons["output_file"].config(text="Set", command=(lambda: self.set_output_file(self.entries["output_file"].get())))

        
        """Parameters"""
        self.param_frame = Frame(self.main_frame)
        self.param_frame.pack(expand=YES, fill=X)
        Label(self.param_frame, text="Parameters").pack()
        
        for name in self.parameters.loc[parameters["type"]!="entry"].index.values:
            block = Frame(self.param_frame)
            self.entries[name] = Entry(block, width=4)
            self.entries[name].pack(side=LEFT)
            Label(block, text=name).pack(side=LEFT)
            if(self.parameters.loc[name, "type"] == "int"):
                var = IntVar()
                var.set(int(self.parameters.loc[name, "value"]))
            elif(self.parameters.loc[name, "type"] == "float"):
                var = DoubleVar()
                var.set(self.parameters.loc[name, "value"])
            self.entries[name].config(textvariable=var)
            self.parameters.loc[name, "value"] = var
            block.pack(side=LEFT, expand=YES, fill=BOTH)
            
        """Buttons"""
        self.buttons_frame = Frame(self.main_frame)
        self.buttons_frame.pack(expand=YES, fill=X)
        
        """Compute"""
        self.button_c = Button(self.buttons_frame, text="Simulate", command=self.simulate)
        self.button_c.pack(side=TOP, expand=YES, fill=X)
        
        """Compute"""
        self.button_c = Button(self.buttons_frame, text="Save raw", command=self.save)
        self.button_c.pack(side=TOP, expand=YES, fill=X)
        
        return
    
    def _open_file(self):
        """
        Open file from the filedialog
        @return Tuple of strings
        """
        file_name = filedialog.askopenfilename(initialdir=self.current_directory, title = "Select Template file")
        if(file_name): self.current_directory = os.path.dirname(file_name)
        return file_name
    
    def update_template_file(self, input_file):
        """
        Loads input genomes 
        """
        if(input_file): #is not None
            self.entries["template"].delete(0, END)
            self.entries["template"].insert(0, str(input_file))
            self.parameters.loc["template", "value"] = input_file
            self.template = load_template_only(self.parameters.loc["template", "value"])
        return
    
    def set_output_file(self, output_file):
        if(os.path.isabs(output_file)):
            self.parameters.loc["output_file", "value"] = output_file
        else:
            self.parameters.loc["output_file", "value"] = os.path.join(self.current_directory,output_file)
            self.entries["output_file"].delete(0, END)
            self.entries["output_file"].insert(0, self.parameters.loc["output_file", "value"])
        print("Output file path updated")
        return
    
    def set_template(self, template):
        #TODO use this func to pass template from matching to simulation
        self.template = template
        self.entries["template"].delete(0, END)
        self.entries["template"].insert(0, str("TEMPLATE FROM LAST MATCHING"))
        return
    
    def simulate(self):
        sim = Simulation(self.template, self.parameters.loc["sample size", "value"].get())
        self.raw_stats, self.cooked_stats = sim.simulate(self.parameters.loc["k", "value"].get(), self.parameters.loc["Beta", "value"].get(), self.parameters.loc["N", "value"].get())
        return
    
    def save(self):
        Simulation.store_raw_data(self.parameters.loc["output_file", "value"], self.raw_stats, self.cooked_stats)
        return
    
    def pack(self):
        self.main_frame.pack(expand=YES, fill=BOTH)
        return