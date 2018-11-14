#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 12:14:47 2018

@author: david
"""

from tkinter import *

class GUI(Frame):
    def __init__(self, parent=None):
        Frame.__init__(self, parent)
        button = Button(self, text="Clone", command=self.clone)
        button.pack()
        return
    
    def clone(self):
        return GUI(self).pack(side="right")


if (__name__=="__main__"):
    mainwin = Tk()
    mainwin.title("Test")
    Label(mainwin, text="I'm you father").pack(expand=YES, fill=X)
    
    #popup = Toplevel()
    #g = GUI(popup).pack()
    mainwin.mainloop()