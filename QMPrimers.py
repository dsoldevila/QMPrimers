
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 15:30:29 2018

@author: David Soldevila
"""

import load_data as ld
import matching as m
from common import *
from tkinter import filedialog
from tkinter import *
import os
import sys
import _thread

parameters = {"gen": "", "primer_pairs": "", "output_file": os.path.join(os.getcwd(),"out.csv"), "forward missmatches": 5, "reverse missmatches": 5, "hanging primers": False} #TODO Matching Type
output_info = {"primerPair": True,"fastaid": True,"primerF": True,"primerR": True,"mismFT": True,"mismRT": True,"amplicon": True,
               "F_pos": True, "mismFT_loc": True, "mismFT_type": True, "R_pos": True, "mismRT_loc": True, "mismRT_type": True}

class GUI(Frame):
    def __init__(self, parent=Frame):
        
        self.main_frame = Frame.__init__(self, parent)
        
        """Other"""
        self.parameters = parameters
        self.current_directory = os.getcwd()
        
        """Menu"""
        self.main_menu = Menu(parent)
        parent.config(menu=self.main_menu)
        file = Menu(self.main_menu, tearoff=False)
        file.add_command(label="New")
        self.main_menu.add_cascade(label="File", menu=file)
        
        """Frame containing Files and Other Parameters frames"""
        self.first_row_frame = Frame(self.main_frame)
        self.first_row_frame.pack(expand=YES, fill=BOTH)
        
        """Select Files Frame"""
        self.file_frame = Frame(self.first_row_frame)
        self.file_frame.pack(side=LEFT, expand=YES, fill=X)
        Label(self.file_frame, text="Files").pack()
        
        self.entries = {}
        self.button = {}
        for pkey in self.parameters:
            if(isinstance(self.parameters[pkey], str)):
                block = Frame(self.file_frame)
                block.pack(expand=YES, fill=BOTH)
                self.entries[pkey] = Entry(block)
                self.entries[pkey].pack(side=LEFT, expand=YES, fill=X)
                self.button[pkey] = Button(block, text="Open")
                self.button[pkey].pack(side=LEFT, expand=NO, fill=X)
                
        self.entries["gen"].insert(0, "<No Genome>")
        self.entries["gen"].bind("<Return>", (lambda event: self.open_bio_files(self.entry_g.get())))
        self.button["gen"].config(command=(lambda: self.update_bio_files(self._open_file())))
        
        self.entries["primer_pairs"].insert(0, "<No Primer Pairs>")
        self.entries["primer_pairs"].bind("<Return>", (lambda event: self.open_csv_file(self.entry_p.get())))
        self.button["primer_pairs"].config(command=(lambda: self.update_primer_file(self._open_file())))
        
        self.entries["output_file"].insert(0, os.path.join(self.current_directory,"output.csv"))
        self.entries["output_file"].bind("<Return>", (lambda event: self.set_output_file(self.entry_out.get())))
        self.button["output_file"].config(text="Set", command=(lambda: self.set_output_file(self.entry_out.get())))
            
            
        """Other parameters"""
        self.other_frame = Frame(self.first_row_frame)
        self.other_frame.pack(side=LEFT, expand=YES, fill=X)
        Label(self.other_frame, text="Parameters").pack()
        self.other_param = {}
        for pkey in self.parameters:
            if(isinstance(self.parameters[pkey], bool)):
                other = BooleanVar()
                Checkbutton(self.other_frame, variable=other, text=pkey).pack(expand=YES)
                other.set(self.parameters[pkey])
                self.other_param[pkey] = other
            elif(isinstance(self.parameters[pkey], int)):
                block = Frame(self.other_frame)
                block.pack(expand=YES)
                Label(block, text=pkey).pack(side=RIGHT)
                other = StringVar()
                Entry(block, textvariable=other).pack(side=LEFT)
                other.set(self.parameters[pkey])
                self.other_param[pkey] = other
        
        """Output info"""
        self.output_frame = Frame(self.main_frame)
        self.output_frame.pack(expand=YES, fill=X)
        Label(self.output_frame, text="Output Info").grid(row=0, column=0)
        Label(self.output_frame, text="(NO EFFECT)").grid(row=0, column=1)
        
        self.output_info = output_info
        c=0
        max_row_len=len(self.output_info)/2
        for key in self.output_info:
            var = BooleanVar()
            Checkbutton(self.output_frame, variable=var, text=key).grid(row=int(1+c/max_row_len), column=int(c%max_row_len), sticky=W)
            var.set(self.output_info[key])
            self.output_info[key] = var
            c = c+1
        
        self.buttons_frame = Frame(self.main_frame)
        self.buttons_frame.pack(expand=YES, fill=X)
        
        """Compute"""
        self.button_c = Button(self.buttons_frame, text="Compute", command=self.compute)
        self.button_c.pack(expand=YES, fill=X)
        
        return
        
    def _open_file(self):
        """
        Open file from the filedialog
        @return Tuple of strings
        """
        file_names = filedialog.askopenfilenames(initialdir=self.current_directory, title = "Select file")
        if(file_names): self.current_directory = os.path.dirname(file_names[0])
        return file_names
    
    def update_bio_files(self, input_files):
        """
        Loads input genomes 
        """
        if(input_files): #is not None
            self.entries["gen"].delete(0, END)
            self.entries["gen"].insert(0, str(input_files))
            self.parameters["gen"] = input_files
        return
    
    def update_primer_file(self, input_files):
        """
        Loads input primer pairs, stored in csv format
        """
        if(input_files): #is not None
            self.entries["primer_pairs"].delete(0, END)
            self.entries["primer_pairs"].insert(0, str(input_files))
            self.parameters["primer_pairs"] = input_files[0]
        return
    
    def set_output_file(self, output_file):
        if(os.path.isabs(output_file)):
            self.parameters["output_file"] = output_file
        else:
            self.parameters["output_file"] = output_file = os.path.join(self.current_directory,output_file)
        return
    
    def compute(self):
        for pkey in self.other_param:
            self.parameters[pkey] = self.other_param[pkey].get()
       
        result = compute(self.parameters)
        
        header = []
        for key in self.output_info:
            if(self.output_info[key].get()):
                header.append(key)
                
        ld.store_matching_results(self.parameters["output_file"], result, header=header)
        print("Finished!")
        return

def get_help():
    parameters_help = {"--help": "Display this list", 
                       "-mf <number>": "Maximum number of missmatches allowed in the forward primer", 
                       "-mr <number>": "Maximum number of missmatches allowed in the reverse primer", 
                       "-gf <path/to/file>": "Location of the genome file, currently no support for multiple files in command line",
                       "-gformat <string>": "(NOT IMPLEMENTED, Optional) Format of the genome file (fasta, etc)", 
                       "-pf </path/to/file>": "Location of the primer pairs, the following header must be included in the file. The order does not matter: \
                           id;forwardPrimer;fPDNA;reversePrimer;rPDNA;ampliconMinLength;ampiconMaxLength", 
                       "--nogui": "GUI is not loaded", 
                       "--hanging-primers": "Primer pairs are allowed to match between [0-mf,len(genome)+mr] instead of just \
                           between the length of the genome",
                       "-info <multiple strings>": "Select which info to output, all info by default"}

    print("QMPRIMERS HELP PAGE")

    for param in parameters_help:
        print(param, ": ", parameters_help[param])
    return

def compute(parameters):
    """
    Manages the program in terminal mode
    """
    gen_record = ld.load_bio_files(parameters["gen"], writable=parameters["hanging primers"])
    primer_pairs = ld.load_csv_file(parameters["primer_pairs"])
    result = m.compute_gen_matching(int(parameters["forward missmatches"]), int(parameters["reverse missmatches"]), primer_pairs, gen_record, hanging_primers=parameters["hanging primers"])

    return result

if (__name__=="__main__"):
    only_cl_param = {"help": False, "command_line": False}
    all_parameters = {**parameters, **only_cl_param}
    flags = {"--help": "help", "--nogui": "command_line", "-mf": "forward missmatches", "-mr": "reverse missmatches", "-gf": "gen", 
             "-gformat": None, "-pf": "primer_pairs", "--hanging-primers": "hanging primers", "-o": "output_file"}
   
    i = 1
    nargs = len(sys.argv)
    
    while i < nargs:
        argv = sys.argv[i]
        if(argv not in flags):
            print("Parameter "+str(sys.argv[i])+" unknown")
            print("Use --help to display the manual")
            exit();
            if(argv=="-info"): #Set custom output
                for key in output_info: output_info[key] = False
                i+=1
                argv = sys.argv[i]
                while(argv[0]!="-"):
                    output_info[argv] = True
        if(argv[:2]=="--"):
            all_parameters[flags[sys.argv[i]]] = True
        elif(argv[0]=="-"):
             all_parameters[flags[sys.argv[i]]] = sys.argv[i+1]
             i+=1
        i+=1  
        
    if(all_parameters[flags["--help"]]):
        get_help()
    elif(all_parameters[flags["--nogui"]]):
        all_parameters["gen"] = (all_parameters["gen"],) #TODO multiple files not implemented in cl
        compute(all_parameters)
    else:
        root = Tk()
        root.title("QMPrimers")
        root.geometry('800x400')
        main_window = GUI(root)
        root.mainloop()
    