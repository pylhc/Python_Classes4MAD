#!/usr/bin/env pythonafs
 
 
from Numeric import *
import sys, os
from string import find, split
from Interpolation import InterpolatingFunction



#########################
class csv:
#########################
    "comma separated read file"

    def __init__(self, filename):
        i=0
        self.time=[]
        self.value=[]
        for line in open(filename):
            i=i+1
            if i>2 and split(line,',')[0] is not "":
                self.time.append(float(split(line,',')[0]))
                self.value.append(float(split(line,',')[1]))

        self.interfunc=InterpolatingFunction((array(self.time)[:],), array(self.value))

