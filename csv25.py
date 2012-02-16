#!/usr/bin/env python2.5
 
from numpy import *
import sys, os
from string import find, split
from scipy import interpolate



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
            if i>2:
                self.time.append(float(split(line,',')[0]))
                self.value.append(float(split(line,',')[1]))

        self.interfunc=interpolate.interp1d(array(self.time), array(self.value))

