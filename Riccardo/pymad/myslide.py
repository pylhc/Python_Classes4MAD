# File: bind1.py

from Tkinter import *

root = Tk()

def callback(e):
    print "clicked at", e.x, e.y
    print e.widget['text']

frame = Frame(root, width=100, height=100)
frame.bind("<B1-Motion>", callback)
frame.pack()

label =Label(root,text="hello")
label.bind("<B1-Motion>", callback)
label.pack()

root.mainloop()
