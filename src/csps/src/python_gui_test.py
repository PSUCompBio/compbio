#!/usr/bin/python

from Tkinter import *

def onclick():
   pass

root = Tk()
text = Text(root)
text.insert(INSERT, "Hello ")
text.insert(END, "world")
text.pack()

root.mainloop()
