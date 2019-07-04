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
import _thread
import pandas as pd

output_info = {}
for key in TEMPLATE_HEADER:
    output_info[key] = True

parameters = [
        ["gen", None, "Genome file dir, no support for multiple files in cl", "-gf", "entry"],
        ["primer_pairs", None, "Primer pairs file dir. A particular header must be used in the file", "-pf", "entry"],
        ["output_file", os.path.join(os.getcwd(),"output"), "Location of the output files, no extension", "-o", "entry"],
        ["forward missmatches", 5, "Maximum number of missmatches allowed on forward primer", "-fm", "param"],
        ["reverse missmatches", 5, "Maximum number of missmatches allowed on reverse primer", "-rm", "param"], 
        ["Nend miss.", 0, "Missmatches in the last N positions on forward and in the first N pos. on reverse ", "-nend", "info"],
        ["hanging primers", False, "Primers allowed to match between [0-mf,len(genome)+mr] instead of just between genome's length", "--hanging", "param"],
        ["check_integrity", False, "Checks integrity of gen files, integrity of primer file is always checked", "--checki", "param"],
        ["check_uppercase", False, "Checks that all gens are in upper case, lower case gens will trigger an integrity file", "--checku", "param"],
        ["csv_template", None, "Precomputed missmatching template", "-i", "entry"],
        ["verbose", False, "Outputs extra information", "--v", "cmd"]]                      
parameters = pd.DataFrame([x[1:] for x in parameters], index = [x[0] for x in parameters], columns=["value", "description", "flag", "type"])



class GUI_matching():
    def __init__(self, parent, gui_simulate):
        
        self.main_frame = Frame(parent)

        #myFont = Font(family="Noto Sans", size=12, weight="bold")
        
        """Other"""
        self.parameters = parameters
        self.current_directory = os.getcwd()
        self.gui_simulate = gui_simulate #pointer to GUI simulation to automatically pass the newest template
        self.op_in_progress = False #Kind of mutex to not allow multiple operations at once
        
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
                    other = IntVar()
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
            
            self.parameters.loc["gen", "value"] = None
            self.parameters.loc["primer_pairs", "value"] = None
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
        
        if(self.op_in_progress==False):
            for pkey in self.other_param:
                self.parameters.loc[pkey, "value"] = self.other_param[pkey].get()
           
            _thread.start_new_thread(self.compute, ())
        else:
            print("An operation is already running")
            
        return
    
    def compute(self):
        self.op_in_progress=True
        self.template, self.discarded, self.gen_record, self.primer_pairs, self.raw_stats, self.cooked_stats = compute(self.parameters)        
        self.gui_simulate.set_template(self.template)        
        self.op_in_progress=False
        return
    
    def load_template_in_thread(self):
        if(self.op_in_progress==False):
            _thread.start_new_thread(self.load_template, ())
        else:
            print("An operation is already running")
        return
    
    def load_template(self):
        self.op_in_progress=True
        self.template, self.discarded, self.gen_record, self.primer_pairs, self.raw_stats, self.cooked_stats = load_template(self.parameters)
        self.gui_simulate.set_template(self.template)
        self.op_in_progress=False
        return
    
    def store_results_in_thread(self):
        if(self.op_in_progress==False):
            _thread.start_new_thread(self.store_results, ())
        else:
            print("An operation is already running")
        return
        
    
    def store_results(self):
        self.op_in_progress=True
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
        
        if(self.parameters.loc["gen", "value"]!=None): 
            input_files = "Fasta = "+self.parameters.loc["gen", "value"]+"     Primer Pairs = "+self.parameters.loc["primer_pairs", "value"]
        else:
            input_files = "From Template: "+self.parameters.loc["csv_template", "value"]
            
        save_matching_info(input_files, self.parameters.loc["output_file", "value"], self.out_template, header, self.discarded, self.out_raw_stats, self.out_cooked_stats)
        self.op_in_progress=False
        return

def get_help(paramaters):
    output_info_keys = [key for key in output_info.keys()]
    other_info = {"-info <multiple strings>": "\tSelect which info to output, all info by default\n\t"+str(output_info_keys),
                  "PRIMER PAIRS HEADER": "tid;forwardPrimer;fPDNA;reversePrimer;rPDNA;ampliconMinLength;ampiconMaxLength"}

    print("QMPRIMERS HELP PAGE")
    pd.options.display.max_colwidth = 100
    print(paramaters.to_string())

    for param in other_info:
        print(param, ": ", other_info[param])
    return

def matching_cl(args):
    parameters.loc["help"] = [False, "Display this list", "--help", ""]
        
    flags = parameters["flag"].values
   
    i = 0
    nargs = len(args)
    
    last_option= None
    
    while i < nargs:
        argv = args[i]
        if(last_option==None):
            if(argv not in flags):
                if(argv=="-info"):
                    last_option = "info"
                    for key in output_info: output_info[key] = False
                else:
                    print("Parameter "+str(sys.argv[i])+" unknown")
                    print("Use --help to display the manual")
                    exit();
            elif(argv[:2]=="--"):
                index = parameters[parameters["flag"]==argv].index
                parameters.loc[index, "value"] = True
            elif(argv[0]=="-"):
                last_option = parameters[parameters["flag"]==argv].index.values
            i+=1
        else:
            if(last_option == "info" and argv in output_info):
                output_info[argv] = True
                i+=1
            elif(last_option in parameters.index.values):
                parameters.loc[last_option, "value"] = argv
                last_option = None
                i+=1
            else:
                last_option = None
                
    set_verbosity(parameters.loc["verbose", "value"])
    parameters.loc["gen", "value"] = (parameters.loc["gen", "value"]) 
    parameters.loc["Nend miss.", "value"] = int(parameters.loc["Nend miss.", "value"]) 

    if(parameters.loc["help", "value"]):
        get_help(parameters)
    else:
        
        if(parameters.loc["csv_template", "value"]):
            template, discarded, gen_record, primer_pairs, raw_stats, cooked_stats = load_template(parameters)

        elif(parameters.loc["gen", "value"]):
            parameters.loc["gen", "value"] = (parameters.loc["gen","value"],) #TODO multiple files not implemented in cl
            template, discarded, gen_record, primer_pairs, raw_stats, cooked_stats = compute(parameters)
            
        i = 0 #TODO convert output_info to list
        header = []
        for key in output_info:
            if(output_info[key]):
                header.append(i)
            i+=1

        Nend = parameters.loc["Nend miss.", "value"]
        if(Nend):
            if(parameters.loc["csv_template", "value"]):
                template = recalculate_Nend(template, primer_pairs, Nend, 0)
            header.extend(["mismFN"+str(Nend), "mismRN"+str(Nend)])

            
        save_matching_info(parameters.loc["output_file", "value"], template, header, discarded, raw_stats, cooked_stats)