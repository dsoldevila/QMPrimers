
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 15:30:29 2018

@author: David Soldevila
"""

from interface import *
from common import *
from tkinter import filedialog
from tkinter import *
import os
import sys
import _thread

parameters = {"gen": "", "primer_pairs": "", "output_file": os.path.join(os.getcwd(),"output.csv"), "forward missmatches": 5, "reverse missmatches": 5, "hanging primers": False} #TODO Matching Type
output_info = {"primerPair": True,"fastaid": True,"primerF": True,"primerR": True,"mismFT": True,"mismRT": True,"amplicon": True, "F_pos": True,
               "mismFT_loc": True, "mismFT_type": True, "mismFT_base": True, "R_pos": True, "mismRT_loc": True, "mismRT_type": True, "mismRT_base": True}

class TextRedirector(object):
    def __init__(self, widget, tag="stdout"):
        self.widget = widget
        self.tag = tag

    def write(self, str):
        self.widget.configure(state="normal")
        self.widget.insert("end", str, (self.tag,))
        self.widget.yview(END)
        
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
        self.entries["gen"].bind("<Return>", (lambda event: self.update_bio_files(self.entries["gen"].get())))
        self.button["gen"].config(command=(lambda: self.update_bio_files(self._open_file())))
        
        self.entries["primer_pairs"].insert(0, "<No Primer Pairs>")
        self.entries["primer_pairs"].bind("<Return>", (lambda event: self.update_csv_file(self.entries["primer_pairs"].get())))
        self.button["primer_pairs"].config(command=(lambda: self.update_primer_file(self._open_file())))
        
        self.entries["output_file"].insert(0, self.parameters["output_file"])
        self.entries["output_file"].bind("<Return>", (lambda event: self.set_output_file(self.entries["output_file"].get())))
        self.button["output_file"].config(text="Set", command=(lambda: self.set_output_file(self.entries["output_file"].get())))
            
            
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
        self.button_c = Button(self.buttons_frame, text="Compute", command=self.compute_in_thread)
        self.button_c.pack(expand=YES, fill=X)
        
        """Save"""
        self.button_s = Button(self.buttons_frame, text="Save", command=self.store_results)
        self.button_s.pack(expand=YES, fill=X)
        
        
        """Terminal"""
        self.terminal_frame = Frame(self.main_frame)
        self.terminal_frame.pack(expand=YES, fill=X)
        self.text = Text(self.terminal_frame, wrap="word")
        self.text.pack(side="top", fill="both", expand=True)
        self.text.tag_configure("stderr", foreground="#b22222")
                            
        sys.stdout = TextRedirector(self.text, "stdout")
        #sys.stderr = TextRedirector(self.text, "stderr")
        
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
    
    def compute_in_thread(self):
        _thread.start_new_thread(self.compute, ())
        return
    
    def compute(self):
        for pkey in self.other_param:
            self.parameters[pkey] = self.other_param[pkey].get()
       
        self.template = compute(self.parameters)
        
        print("Finished!")
        return
    
    def store_results(self):
        header = []
        for key in self.output_info:
            if(self.output_info[key].get()):
                header.append(key)
        save_template_primer_missmatches(self.parameters["output_file"], self.template, header=header)
        print("Saved")
        return

def get_help():
    parameters_help = {"--help": "Display this list", 
                       "-mf <number>": "\t\tMaximum number of missmatches allowed in the forward primer, default 5", 
                       "-mr <number>": "\t\tMaximum number of missmatches allowed in the reverse primer, default 5", 
                       "(*) -gf <path/to/file>": "\tLocation of the genome file, currently no support for multiple files in command line",
                       "-gformat <string>": "\t\t(NOT IMPLEMENTED, Optional) Format of the genome file (fasta, etc)", 
                       "(*) -pf </path/to/file>": "\tLocation of the primer pairs, the following header must be included in the file. The order does not matter:",
                       "\tFORMAT": "\t\tid;forwardPrimer;fPDNA;reversePrimer;rPDNA;ampliconMinLength;ampiconMaxLength", 
                       "--nogui": "\t\t\tGUI is not loaded", 
                       "--hanging-primers": "\t\tPrimer pairs are allowed to match between [0-mf,len(genome)+mr] instead of just between the length of the genome",
                       "-info <multiple strings>": "\tSelect which info to output, all info by default\n\t"+str(output_info.keys()),
                       "-o </path/to/output": "\t\t./out.csv is the default path",
                       "(*)": "Parameters marked with this are mandatory, the other either have a default value or are optional"}

    print("QMPRIMERS HELP PAGE")

    for param in parameters_help:
        print(param, ": ", parameters_help[param])
    return

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
            if(argv=="-info"): #Set custom output
                for key in output_info: output_info[key] = False
                i+=1
                argv = sys.argv[i]
                while(argv[0]!="-"):
                    output_info[argv] = True
            else:
                print("Parameter "+str(sys.argv[i])+" unknown")
                print("Use --help to display the manual")
                exit();
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
        template = compute(all_parameters)
        header = []
        for key in output_info:
            if(output_info[key]):
                header.append(key)
        save_template_primer_missmatches(parameters["output_file"], template, header=header)
    else:
        saved_sys_stdout = sys.stdout
        saved_sys_stderr = sys.stderr
        root = Tk()
        root.title("QMPrimers")
        root.geometry('800x400')
        main_window = GUI(root)
        root.mainloop()
        sys.stdout = saved_sys_stdout
        sys.stderr = saved_sys_stderr