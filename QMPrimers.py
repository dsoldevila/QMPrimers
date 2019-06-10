
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
from tkinter.font import Font
import tkinter.ttk as ttk
import os
import sys
import pandas as pd
from simulation_frontend import *
from matching_frontend import *


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

        
        """Tabs"""
        nb = ttk.Notebook(self.main_frame)
        
        page2 = ttk.Frame(nb)
        self.gui_simulate = GUI_simulate(page2)
        
        page1 = ttk.Frame(nb)
        self.gui_matching = GUI_matching(page1, self.gui_simulate)
        
        self.gui_matching.pack()
        nb.add(page1, text="Matching")
        
        
        self.gui_simulate.pack()
        nb.add(page2, text="Simulation")

        nb.pack(fill=BOTH, expand=YES)
        
        

        
        """Terminal"""
        self.terminal_frame = Frame(self.main_frame)
        self.terminal_frame.pack(side=BOTTOM, expand=YES, fill=X)
        self.text = Text(self.terminal_frame, wrap="word")
        self.text.pack(fill=BOTH, expand=True)
        self.text.tag_configure("stderr", foreground="#b22222")
                            
        sys.stdout = TextRedirector(self.text, "stdout")
        #sys.stderr = TextRedirector(self.text, "stderr")
        
    

def get_help():
    info = {"PARAMETERS" : "\t-------------",
            "--help": "\tDisplay this message",
            "--sim": "\tSimulation mode",
            "--match": "\tMatchinf mode",
            "GUI mode": "To open the graphical mode, do not pass any parameter"}

    print("QMPRIMERS HELP PAGE");

    for param in info:
        print(param, ": ", info[param])
    return

if (__name__=="__main__"):
    if(len(sys.argv)>1 ):
        if(sys.argv[1]=="--sim"):
            sim_cl(sys.argv[2:])
        elif(sys.argv[1]=="--match"):
            matching_cl(sys.argv[2:])
            
        elif(sys.argv[1]=="--help"):
            get_help()
        else:
            print("Unknown command ",sys.argv[2],". Use --help to display the manual.")
            
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