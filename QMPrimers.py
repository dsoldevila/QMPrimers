
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 15:30:29 2018

@author: David Soldevila
"""

from tkinter import filedialog
from tkinter import *
from tkinter.font import Font
import tkinter.ttk as ttk
import os
#import sys
import pandas as pd
from simulation_frontend import *
from matching_frontend import *
import queue

# The new Stream Object which replaces the default stream associated with sys.stdout
# This object just puts data in a queue!
class WriteStream(object):
    def __init__(self,queue):
        self.queue = queue

    def write(self, text):
        self.queue.put(text)

class TextRedirector(object):
    def __init__(self, widget, tag="stdout"):
        self.widget = widget
        self.tag = tag

    def write(self, str):
        self.widget.configure(state="normal")
        self.widget.insert("end", str, (self.tag,))
        self.widget.yview(END)
        
class GUI(Frame):
    def __init__(self, parent, match_mosi_queue, match_miso_queue):
        
        parent.protocol('WM_DELETE_WINDOW', (lambda: self.finnish(parent))) #overwrites the close button function
        Frame.__init__(self, parent)
        self.main_frame = parent
        
        """Other"""
        self.current_directory = os.getcwd()
        self.queue = queue  
        """Menu Bar demo"""
        """
        self.main_menu = Menu(parent)
        parent.config(menu=self.main_menu)
        file = Menu(self.main_menu, tearoff=False)
        file.add_command(label="Option 1")
        file.add_command(label="Option 2...")
        self.main_menu.add_cascade(label="File", menu=file)
        """
        
        """Tabs"""
        nb = ttk.Notebook(self.main_frame)
        
        page2 = ttk.Frame(nb)
        self.gui_simulate = GUI_simulate(page2)
        
        page1 = ttk.Frame(nb)
        self.gui_matching = GUI_matching(page1, self.gui_simulate, match_mosi_queue, match_miso_queue)
        
        self.gui_matching.pack()
        nb.add(page1, text="Matching")
        
        
        self.gui_simulate.pack()
        nb.add(page2, text="Simulation")

        nb.pack(fill=BOTH, expand=NO)
        
        """Terminal"""
        self.terminal_frame = Frame(self.main_frame)
        self.text = Text(self.terminal_frame, wrap="word")
        self.text.pack(fill=BOTH, expand=YES)
        self.text.tag_configure("stderr", foreground="#b22222")
                                
        """STDOUT redirect"""
        #sys.stdout = TextRedirector(self.text, "stdout")
        #sys.stderr = TextRedirector(self.text, "stderr")

        """Verbose"""
        self.extra_frame = Frame(self.main_frame)
        self.is_verbose = BooleanVar()
        Checkbutton(self.extra_frame, variable=self.is_verbose, text="verbose", command=self.update_verbosity).pack(expand=NO)
        self.is_verbose.set(False)
        init_logger()
        self.update_verbosity()
    
    
        #Packing Terminal and Verbose
        self.extra_frame.pack(side=BOTTOM, expand=NO, fill=X)
        self.terminal_frame.pack(side=BOTTOM, expand=YES, fill=BOTH)
        
        self.main_frame.after(1000, self.update_terminal)


    def update_terminal(self):
        s = self.queue.qsize()
        if(s>0):
            items = ""
            for i in range(s):
                items+=self.queue.get()
            self.text.configure(state="normal")
            self.text.insert("end", items, ("stdout",))
            self.text.yview(END)
            
        self.main_frame.after(1000, self.update_terminal)
    
    def update_verbosity(self):
        set_verbosity(self.is_verbose.get())
        
    def finnish(self, parent):
        close_logger()
        parent.destroy()
    

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
        init_logger()
        if(sys.argv[1]=="--sim"):
            sim_cl(sys.argv[2:])
        elif(sys.argv[1]=="--match"):
            matching_cl(sys.argv[2:])
            
        elif(sys.argv[1]=="--help"):
            get_help()
        else:
            print("Unknown command ",sys.argv[2],". Use --help to display the manual.")
        close_logger() 
    else:
        saved_sys_stdout = sys.stdout
        saved_sys_stderr = sys.stderr
        
        # Create Queue and redirect sys.stdout to this queue
        queue = queue.Queue()
        sys.stdout = WriteStream(queue)
        
        #create matching thread
        match_miso_queue = queue.Queue()
        match_mosi_queue = queue.Queue()
        matching = matching_thread_wrapper(match_mosi_queue, match_miso_queue)
        matching.run()
        
        root = Tk()
        root.title("QMPrimers")
        root.geometry('900x400')
        main_window = GUI(root, match_mosi_queue, match_miso_queue)
        root.mainloop()
        sys.stdout = saved_sys_stdout
        sys.stderr = saved_sys_stderr