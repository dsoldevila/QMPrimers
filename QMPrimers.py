
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
from GUI import *

output_info = {}
for key in TEMPLATE_HEADER:
    output_info[key] = True

parameters = [
        ["gen", None, "Genome file dir, no support for multiple files in cl", "-gf", "entry"],
        ["primer_pairs", None, "Primer pairs file dir. A particular header must be used in the file", "-pf", "entry"],
        ["output_file", os.path.join(os.getcwd(),"output.csv"), "Location of the output file", "-o", "entry"],
        ["forward missmatches", 5, "Maximum number of missmatches allowed on forward primer", "-fm", "param"],
        ["reverse missmatches", 5, "Maximum number of missmatches allowed on reverse primer", "-rm", "param"], 
        ["Nend miss.", 0, "Missmatches in the last N positions on forward and in the first N pos. on reverse ", "-nend", "info"],
        ["hanging primers", False, "Primers allowed to match between [0-mf,len(genome)+mr] instead of just between genome's length", "--hanging", "param"],
        ["check_integrity", False, "Checks integrity of gen files, integrity of primer file is always checked", "--checki", "param"],
        ["check_uppercase", False, "Checks that all gens are in upper case, lower case gens will trigger an integrity file", "--checku", "param"],
        ["csv_template", None, "Precomputed missmatching template", "-i", "entry"]]
                        
parameters = pd.DataFrame([x[1:] for x in parameters], index = [x[0] for x in parameters], columns=["value", "description", "flag", "type"])

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
        file.add_command(label="Option 1")
        file.add_command(label="Option 2...")
        self.main_menu.add_cascade(label="File", menu=file)
        
        
        self.gui_compute = GUI_compute(self.main_frame, self.parameters, output_info)
        self.gui_compute.pack()
        
        """Terminal"""
        self.terminal_frame = Frame(self.main_frame)
        self.terminal_frame.pack(side=BOTTOM, expand=YES, fill=X)
        self.text = Text(self.terminal_frame, wrap="word")
        self.text.pack(fill=BOTH, expand=True)
        self.text.tag_configure("stderr", foreground="#b22222")
                            
        sys.stdout = TextRedirector(self.text, "stdout")
        #sys.stderr = TextRedirector(self.text, "stderr")
        
    

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
        ["help", False, "Display this list", "--help", ""],
        ["command_line", False, "Triggers the command line mode", "--nogui", ""]]
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
    cl_parameters.loc["Nend miss.", "value"] = int(cl_parameters.loc["Nend miss.", "value"]) 

    if(cl_parameters.loc["help", "value"]):
        get_help(cl_parameters)
    elif(cl_parameters.loc["command_line","value"]):
        
        if(cl_parameters.loc["gen", "value"]):
            cl_parameters.loc["gen", "value"] = (cl_parameters.loc["gen","value"],) #TODO multiple files not implemented in cl
            template, gen_record, primer_pair = compute(cl_parameters)
        elif(cl_parameters.loc["csv_template", "value"]):
            template, gen_record, primer_pair = load_template(cl_parameters)
            
        header = []
        for key in output_info:
            if(output_info[key]):
                header.append(key)
        Nend = cl_parameters.loc["Nend miss.", "value"]
        if(Nend):
            header.extend(["mismFN"+str(Nend), "mismRN"+str(Nend)])
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