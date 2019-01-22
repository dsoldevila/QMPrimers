
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
        self.entry_g.bind("<Return>", (lambda event: self.open_bio_files(self.entry_g.get())))
        
        self.button_g = Button(self.file_frame, text="Open", command=(lambda: self.open_bio_files(self._open_file())))
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
        file_names = filedialog.askopenfilenames(initialdir=self.current_directory, title = "Select file")
        if(file_names): self.current_directory = os.path.dirname(file_names[0])
        return file_names
    
    def open_bio_files(self, input_files):
        if(input_files): #is not None
            self.entry_g.delete(0, END)
            self.entry_g.insert(0, str(input_files))
            self.gen_record = ld.load_bio_files(list(input_files))
        return
    
    def open_csv_file(self, input_files):
        if(input_files): #is not None
            self.entry_p.delete(0, END)
            self.entry_p.insert(0, str(input_files))
            print(input_files)
            self.primer_pairs = ld.load_csv_file(input_files[0])
            print(self.primer_pairs[0].f.seq+"\n")
        return
    
    def compute(self):
        result = m.compute_gen_matching(10, 10, self.primer_pairs, self.gen_record)
        print(result[0])
        return

def get_help():
    print("QMPrimers help page")
    return

def compute_from_cmd(parameters):
    gen_record = ld.load_bio_files([parameters["-gf"]],parameters["-gformat"])
    primer_pairs = ld.load_csv_file(parameters["-pf"])
    result = m.compute_gen_matching(int(parameters["-mf"]), int(parameters["-mr"]), primer_pairs, gen_record)

    return result

if (__name__=="__main__"):
    
    parameters = {"--help": False, "-mf": 10, "-mr": 10, "-gf": None, "-gformat": None, "-pf": None, "--nogui": False}
    
    i = 1
    nargs = len(sys.argv)
    
    while i < nargs:
        if(sys.argv[i] not in parameters):
            print("Parameter "+str(sys.argv[i])+" unknown")
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
        compute_from_cmd(parameters)
    else:
        root = Tk()
        root.title("QMPrimers")
        main_window = GUI(root)
        root.mainloop()
    