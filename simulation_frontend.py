#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 15:45:06 2019

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

parameters = [
        ["template", None, "Precomputed missmatching template", "-i", "entry"],
        ["sample size", 10, "Nº genome samples per simulation step", "-s", "int"],
        ["Beta", 4, "<No description currently>", "-b", "int"],
        ["k", 0.5, "Geometric proportion parameter. Range(0,1)", "-k", "float"],
        ["N", 100, "Nº simulation steps", "-n", "int"],
        ["output_file", os.path.join(os.getcwd(),"simout"), "Location of the output files, no extension", "-o", "entry"],
        ["Confidence Interval", 0.95, "CI used to compute the statistics", "-ci", "float"],
        ["verbose", False, "Outputs extra information", "--v", "cmd"]]

parameters = pd.DataFrame([x[1:] for x in parameters], index = [x[0] for x in parameters], columns=["value", "description", "flag", "type"])

class GUI_simulate():
    def __init__(self, parent):
        
        self.main_frame = Frame(parent)
        
        """Other"""
        self.parameters = parameters[parameters["type"]!="cmd"]
        self.current_directory = os.getcwd()
        self.entries = {}
        self.buttons = {}
        self.template = None
        self.is_simulating = False #Kind of mutex to not allow multiple simulations at once
        
        """Frame containing Files and Other Parameters frames"""
        self.first_row_frame = Frame(self.main_frame)
        self.first_row_frame.pack(expand=YES, fill=BOTH)
        
        """Select Files Frame"""
        self.file_frame = Frame(self.first_row_frame)
        self.file_frame.pack(side=LEFT, expand=YES, fill=X)
        Label(self.file_frame, text="Files").pack()
        
        for name in self.parameters.loc[self.parameters["type"]=="entry"].index.values:
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
        
        for name in self.parameters.loc[self.parameters["type"]!="entry"].index.values:
            block = Frame(self.param_frame)
            self.entries[name] = Entry(block, width=4)
            self.entries[name].pack(side=LEFT)
            Label(block, text=name).pack(side=LEFT)
            if(self.parameters.loc[name, "type"] == "int"):
                logging.debug(str(name)+" loaded as int")
                var = IntVar()
                var.set(int(self.parameters.loc[name, "value"]))
            elif(self.parameters.loc[name, "type"] == "float"):
                logging.debug(str(name)+" loaded as float")
                var = DoubleVar()
                var.set(self.parameters.loc[name, "value"])
            self.entries[name].config(textvariable=var)
            self.parameters.at[name, "value"] = var
            block.pack(side=LEFT, expand=YES, fill=BOTH)
            
        """Buttons"""
        self.buttons_frame = Frame(self.main_frame)
        self.buttons_frame.pack(expand=YES, fill=X)
        
        """Compute"""
        self.button_c = Button(self.buttons_frame, text="Simulate", command=self.simulate_in_thread)
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
            self.parameters.at["template", "value"] = input_file
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
    
    def simulate_in_thread(self):
        if(self.is_simulating==False):
            _thread.start_new_thread(self.simulate, ())
        return
    
    def simulate(self):
        self.is_simulating = True
        sim = Simulation()
        try:
            self.raw_stats, self.cooked_stats = sim.simulate(self.template, int(self.parameters.loc["sample size", "value"].get()),
                                                             self.parameters.loc["k", "value"].get(), 
                                                             self.parameters.loc["Beta", "value"].get(), 
                                                             self.parameters.loc["N", "value"].get(), 
                                                             self.parameters.loc["Confidence Interval", "value"].get())
        except(Exception) as e:
            logging.error(e)
        
        self.is_simulating = False
        return
    
    def save(self):
        Simulation.store_data(self.parameters.loc["output_file", "value"], self.raw_stats, self.cooked_stats, 
                                  self.parameters.loc["template", "value"], self.parameters.loc["sample size", "value"].get(),
                                  self.parameters.loc["k", "value"].get(), self.parameters.loc["Beta", "value"].get(), 
                                  self.parameters.loc["N", "value"].get())
        return
    
    def pack(self):
        self.main_frame.pack(expand=YES, fill=BOTH)
        return
    
def get_help(paramaters):
    print("QMPRIMERS SIMULATION HELP PAGE")
    pd.options.display.max_colwidth = 100
    print(paramaters.to_string())

    return

def sim_cl(args):
    
    parameters.loc["help"] = [False, "Display this list", "--help", ""]
        
    flags = parameters["flag"].values
   
    i = 0
    nargs = len(args)
    
    last_option= None
    
    while i < nargs:
        argv = args[i]
        if(last_option==None):
            if(argv not in flags):
                print("Parameter "+str(argv)+" unknown")
                print("Use --help to display the manual")
                exit();
            elif(argv[:2]=="--"):
                index = parameters[parameters["flag"]==argv].index
                parameters.loc[index, "value"] = True
            elif(argv[0]=="-"):
                last_option = parameters[parameters["flag"]==argv].index.values
            i+=1
        else:
            if(last_option in parameters.index.values):
                parameters.loc[last_option, "value"] = argv
                last_option = None
                i+=1
            else:
                last_option = None
                
    set_verbosity(parameters.loc["verbose", "value"])
    
    if(parameters.loc["help", "value"]):
        get_help(parameters)
    else:
        template = load_template_only(parameters.loc["template","value"])
        if(template.empty): return
        
        sim = Simulation()
        
        try: #wrap exception of type int("a")
            raw_stats, cooked_stats = sim.simulate(template, int(parameters.loc["sample size", "value"]), float(parameters.loc["k", "value"]), 
                                                   int(parameters.loc["Beta", "value"]), int(parameters.loc["N", "value"]),  
                                                   float(parameters.loc["Confidence Interval", "value"]))
        except(Exception) as e:
            logging.error(e)
            return
        
        Simulation.store_data(parameters.loc["output_file", "value"], raw_stats, cooked_stats, 
                              parameters.loc["template", "value"], parameters.loc["sample size", "value"],
                              parameters.loc["k", "value"], parameters.loc["Beta", "value"], parameters.loc["N", "value"])
    return