
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


class GUI(Frame):
    def __init__(self, parent=Frame):
        
        self.main_frame = Frame.__init__(self, parent)
        
        self.parameters = {"gen": None, "primer_pairs": None, "output_file": "out.csv", "mf": 5, "mr": 5, "hanging-primers": False, "matching_type(NOT IMPLEMENTED": 0}
        
        """Menu"""
        self.main_menu = Menu(parent)
        parent.config(menu=self.main_menu)
        file = Menu(self.main_menu, tearoff=False)
        file.add_command(label="New")
        self.main_menu.add_cascade(label="File", menu=file)
        
        """Select Files Frame"""
        self.file_frame = Frame(self.main_frame)
        self.file_frame.pack(expand=YES, fill=X)
        
        self.entry_g = Entry(self.file_frame)
        self.entry_g.pack(side=LEFT, expand=YES, fill=X)
        self.entry_g.insert(0, "<No Genome>")
        self.entry_g.bind("<Return>", (lambda event: self.open_bio_files(self.entry_g.get())))
        
        self.button_g = Button(self.file_frame, text="Open", command=(lambda: self.update_bio_files(self._open_file())))
        self.button_g.pack(side=LEFT, expand=NO, fill=X)
        
        self.entry_p = Entry(self.file_frame)
        self.entry_p.pack(side=LEFT, expand=YES, fill=X)
        self.entry_p.insert(0, "<No Primer Pairs>")
        self.entry_p.bind("<Return>", (lambda event: self.open_csv_file(self.entry_p.get())))
        
        self.button_p = Button(self.file_frame, text="Open", command=(lambda: self.update_primer_file(self._open_file())))
        self.button_p.pack(side=LEFT, expand=NO, fill=X)
        
        self.entry_out = Entry(self.file_frame)
        self.entry_out.pack(side=LEFT, expand=YES, fill=X)
        self.entry_out.insert(0, "output.csv")
        self.entry_out.bind("<Return>", (lambda event: self.set_output_file(self.entry_out.get())))
        
        self.button_out = Button(self.file_frame, text="Set", command=(lambda: self.set_output_file(self.entry_out.get())))
        self.button_out.pack(side=LEFT, expand=NO, fill=X)
        
            
        """Other parameters"""
        self.other_frame = Frame(self.main_frame)
        self.other_frame.pack(expand=YES, fill=X)
        self.other_param_dict = {}
        for pkey in self.parameters:
            if(isinstance(self.parameters[pkey], bool)):
                other = BooleanVar()
                Checkbutton(self.other_frame, variable=other, text=pkey).pack(side=LEFT, expand=YES)
                other.set(self.parameters[pkey])
                self.other_param_dict[pkey] = other
            elif(isinstance(self.parameters[pkey], int)):
                block = Frame(self.other_frame)
                block.pack(side=LEFT, expand=YES)
                Label(block, text=pkey).pack(side=LEFT)
                other = StringVar()
                Entry(block, textvariable=other).pack(side=LEFT, expand=YES)
                other.set(self.parameters[pkey])
                self.other_param_dict[pkey] = other
        
        
        
        """Compute"""
        self.button_c = Button(text="Compute", command=self.compute)
        self.button_c.pack(side=BOTTOM, expand=YES, fill=BOTH)
        
        """Other"""
        self.current_directory = os.getcwd()
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
            self.entry_g.delete(0, END)
            self.entry_g.insert(0, str(input_files))
            self.parameters["gen"] = input_files
        return
    
    def update_primer_file(self, input_files):
        """
        Loads input primer pairs, stored in csv format
        """
        if(input_files): #is not None
            self.entry_p.delete(0, END)
            self.entry_p.insert(0, str(input_files))
            self.parameters["primer_pairs"] = input_files[0]
        return
    
    def set_output_file(self, output_file):
        if(os.path.isabs(output_file)):
            self.parameters["output_file"] = output_file
        else:
            self.parameters["output_file"] = output_file = os.path.join(self.current_directory,output_file)
        return
    
    def compute(self):
        for pkey in self.other_param_dict:
            self.parameters[pkey] = self.other_param_dict[pkey].get()
        print(self.parameters)
        self.gen_record = ld.load_bio_files(self.parameters["gen"], writable=self.parameters["hanging-primers"])
        self.primer_pairs = ld.load_csv_file(self.parameters["primer_pairs"])
        result = m.compute_gen_matching(int(self.parameters["mf"]), int(self.parameters["mr"]), self.primer_pairs, self.gen_record,
                                        hanging_primers=self.parameters["hanging-primers"])
        header = ["primerPair","fastaid","primerF","primerR","mismFT","mismRT","amplicon", "F_pos", "mismFT_loc", "mismFT_type", "R_pos", "mismRT_loc", "mismRT_type"]
        ld.store_matching_results(self.parameters["output_file"], result, header=header)
        print("Finished!")
        return

def get_help():
    parameters_help = {"--help": "Display this list", 
                       "-mf <number>": "Maximum number of missmatches allowed in the forward primer", 
                       "-mr <number>": "Maximum number of missmatches allowed in the reverse primer", "-gf <path/to/file>": "Location of the genome file",
                       "-gformat <string>": "(Optional) Format of the genome file (fasta, etc)", 
                       "-pf </path/to/file>": "Location of the primer pairs, the following header must be included in the file. The order does not matter: \
                           id;forwardPrimer;fPDNA;reversePrimer;rPDNA;ampliconMinLength;ampiconMaxLength", 
                       "--nogui": "GUI is not loaded", 
                       "--hanging-primers": "Primer pairs are allowed to match between [0-mf,len(genome)+mr] instead of just \
                           between the length of the genome"}

    print("QMPRIMERS HELP PAGE")

    for param in parameters_help:
        print(param, ": ", parameters_help[param])
    return

def compute_from_cl(parameters):
    """
    Manages the program in terminal mode
    """
    gen_record = ld.load_bio_files([parameters["-gf"]],parameters["-gformat"], writable=parameters["--hanging-primers"])
    primer_pairs = ld.load_csv_file(parameters["-pf"])
    result = m.compute_gen_matching(int(parameters["-mf"]), int(parameters["-mr"]), primer_pairs, gen_record, hanging_primers=parameters["--hanging-primers"])

    return result

if (__name__=="__main__"):
    parameters = {"--help": False, "--nogui": False, "-mf": 5, "-mr": 5, "-gf": None, "-gformat": None, "-pf": None, "--hanging-primers": False}
    i = 1
    nargs = len(sys.argv)
    
    while i < nargs:
        if(sys.argv[i] not in parameters):
            print("Parameter "+str(sys.argv[i])+" unknown")
            print("Use --help to display the manual")
            exit();
        if(sys.argv[i][:2]=="--"):
            parameters[sys.argv[i]] = True
        elif(sys.argv[i][0]=="-"):
             parameters[sys.argv[i]] = sys.argv[i+1]
             i+=1
        i+=1  
        
    if(parameters["--help"]):
        get_help()
    elif(parameters["--nogui"]):
        compute_from_cl(parameters)
    else:
        root = Tk()
        root.title("QMPrimers")
        main_window = GUI(root)
        root.mainloop()
    