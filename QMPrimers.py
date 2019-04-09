
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
import pandas as pd

output_info = {}
for key in TEMPLATE_HEADER:
    output_info[key] = True

parameters = [
        ["gen", "", "Genome file dir, no support for multiple files in cl", "-gf"],
        ["primer_pairs", "", "Primer pairs file dir. A particular header must be used in the file", "-pf"],
        ["output_file", os.path.join(os.getcwd(),"output.csv"), "Location of the output file", "-o"],
        ["forward missmatches", 5, "Maximum number of missmatches allowed in the forward primer", "-fm"],
        ["reverse missmatches", 5, "Maximum number of missmatches allowed in the reverse primer", "-rm"], 
        ["hanging primers", False, "Primers allowed to match between [0-mf,len(genome)+mr] instead of just between genome's length", "--hanging"],
        ["check_integrity", False, "Checks integrity of gen files, integrity of primer file is always checked", "--checki"],
        ["check_uppercase", False, "Checks that all gens are in upper case, lower case gens will trigger an integrity file", "--checku"]]
                        
parameters = pd.DataFrame([x[1:] for x in parameters], index = [x[0] for x in parameters], columns=["value", "description", "flag"])

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
        for name in self.parameters.index:
            if(isinstance(self.parameters.loc[name, "value"], str)):
                block = Frame(self.file_frame)
                block.pack(expand=YES, fill=BOTH)
                self.entries[name] = Entry(block)
                self.entries[name].pack(side=LEFT, expand=YES, fill=X)
                self.button[name] = Button(block, text="Open")
                self.button[name].pack(side=LEFT, expand=NO, fill=X)
                
        self.entries["gen"].insert(0, "<No Genome>")
        self.entries["gen"].bind("<Return>", (lambda event: self.update_bio_files(self.entries["gen"].get())))
        self.button["gen"].config(command=(lambda: self.update_bio_files(self._open_file())))
        
        self.entries["primer_pairs"].insert(0, "<No Primer Pairs>")
        self.entries["primer_pairs"].bind("<Return>", (lambda event: self.update_csv_file(self.entries["primer_pairs"].get())))
        self.button["primer_pairs"].config(command=(lambda: self.update_primer_file(self._open_file())))
        
        self.entries["output_file"].insert(0, self.parameters.loc["output_file", "value"])
        self.entries["output_file"].bind("<Return>", (lambda event: self.set_output_file(self.entries["output_file"].get())))
        self.button["output_file"].config(text="Set", command=(lambda: self.set_output_file(self.entries["output_file"].get())))
            
            
        """Other parameters"""
        self.other_frame = Frame(self.first_row_frame)
        self.other_frame.pack(side=LEFT, expand=YES, fill=X)
        Label(self.other_frame, text="Parameters").pack()
        self.other_param = {}
        for name in self.parameters.index:
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
                Entry(block, textvariable=other).pack(side=LEFT)
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
            self.parameters.loc["gen", "value"] = input_files
        return
    
    def update_primer_file(self, input_files):
        """
        Loads input primer pairs, stored in csv format
        """
        if(input_files): #is not None
            self.entries["primer_pairs"].delete(0, END)
            self.entries["primer_pairs"].insert(0, str(input_files))
            self.parameters.loc["primer_pairs", "value"] = input_files[0]
        return
    
    def set_output_file(self, output_file):
        if(os.path.isabs(output_file)):
            self.parameters.loc["output_file", "value"] = output_file
        else:
            self.parameters.loc["output_file", "value"] = output_file = os.path.join(self.current_directory,output_file)
        return
    
    def compute_in_thread(self):
        _thread.start_new_thread(self.compute, ())
        return
    
    def compute(self):
        for pkey in self.other_param:
            self.parameters.loc[pkey, "value"] = self.other_param[pkey].get()
       
        self.template = compute(self.parameters)
        
        print("Finished!")
        return
    
    def store_results(self):
        header = []
        for key in self.output_info:
            if(self.output_info[key].get()):
                header.append(key)
        save_template_primer_missmatches(self.parameters.loc["output_file", "value"], self.template, header=header)
        print("Saved")
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

if (__name__=="__main__"):
    only_cl_parameters = [
        ["help", False, "Display this list", "--help"],
        ["command_line", False, "Triggers the command line mode", "--nogui"]]
    cl_parameters = parameters.copy()
    for param in only_cl_parameters:
        cl_parameters.loc[param[0]] = param[1:]
        
    flags = cl_parameters["flag"].values
   
    i = 1
    nargs = len(sys.argv)
    
    last_option= None
    
    while i < nargs:
        argv = sys.argv[i]
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
                index = cl_parameters[cl_parameters["flag"]==sys.argv[i]].index
                cl_parameters.loc[index, "value"] = True
            elif(argv[0]=="-"):
                last_option = cl_parameters[cl_parameters["flag"]==sys.argv[i]].index.values
                #index = cl_parameters[cl_parameters["flag"]==sys.argv[i]].index
                #cl_parameters.loc[index, "value"] = sys.argv[i+1]
            i+=1
        else:
            if(last_option == "info" and argv in output_info):
                output_info[argv] = True
                i+=1
            elif(last_option in cl_parameters.index.values):
                cl_parameters.loc[last_option, "value"] = argv
                last_option = None
                i+=1
            else:
                last_option = None
    
    cl_parameters.loc["gen", "value"] = (cl_parameters.loc["gen", "value"]) 
    
    """
    while i < nargs:
        argv = sys.argv[i]
        if(argv not in flags):
            
            if(argv=="-info"): #Set custom output
                for key in output_info: output_info[key] = False
                i+=1
                argv = sys.argv[i]
                while(sys.argv[i+1]!="-"):
                    output_info[argv] = True
                    i+=1
                    argv = sys.argv[i]
                output_info[argv] = True
                print(output_info)
            
            else:
                print("Parameter "+str(sys.argv[i])+" unknown")
                print("Use --help to display the manual")
                exit();
        if(argv[:2]=="--"):
            index = cl_parameters[cl_parameters["flag"]==sys.argv[i]].index
            cl_parameters.loc[index, "value"] = True
        elif(argv[0]=="-"):
            index = cl_parameters[cl_parameters["flag"]==sys.argv[i]].index
            cl_parameters.loc[index, "value"] = sys.argv[i+1]
            i+=1
        i+=1  
    """ 
    if(cl_parameters.loc["help", "value"]):
        get_help(cl_parameters)
    elif(cl_parameters.loc["command_line","value"]):
        cl_parameters.loc["gen", "value"] = (cl_parameters.loc["gen","value"],) #TODO multiple files not implemented in cl
        template = compute(cl_parameters)
        header = []
        for key in output_info:
            if(output_info[key]):
                header.append(key)
        save_template_primer_missmatches(cl_parameters.loc["output_file", "value"], template, header=header)
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