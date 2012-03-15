#-------------------------------------------------------------------------------
# This file is part of PyMad.
#
# Copyright (c) 2011, CERN. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# 	http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#-------------------------------------------------------------------------------
# -*- coding: utf-8 -*-
'''
.. module:: tfs

Function to read tfs tables into Python objects

.. moduleauthor:: Yngve Inntjore Levinsen <Yngve.Inntjore.Levinsen@cern.ch>
'''
import numpy
import os
from string import lower


class LookupDict:
    '''
    A dictionary-like structure, which exposes the values of the keys also as attributes with the key names
    '''

    def __init__(self, values):
        ''' Initializes the class with the values.

        Parameters:
        values -- A dictionary with strings as keys and lists as values
        '''
        # we store the values in a new dict internally, to unify the keys to lowercase
        self._values = dict()
        for key, val in values.items():
            self._values[self._unify_key(key)] = val

    def __iter__(self):
        return iter(self._values)

    def _get_val_or_raise_error(self, key, error):
        ukey = self._unify_key(key)

        if (self._values.has_key(ukey)):
            return self._values[key]
        else:
            raise(error)

    def __getattr__(self, name):
        ''' Exposes the variables as attributes. This allows to use code like the following:

        tfs = TfsTable(...)
        print tfs.x

        '''
        return self._get_val_or_raise_error(name, AttributeError())

    def __getitem__(self, key):
        ''' Emulates the [] operator, so that TfsTable can be used just like a dictionary.
            The keys are considered case insensitive.
        '''
        return self._get_val_or_raise_error(key, KeyError())

    def _unify_key(self, key):
        return lower(key)

    def keys(self):
        '''
         Similar to dictionary.keys()...
        '''
        return self._values.keys()


class TfsTable(LookupDict):
    ''' A class to hold the results of a twiss '''
    def __init__(self, values):
        LookupDict.__init__(self, values)
        if self._values.has_key('name'):
            self._names = self._values['name']

    @property
    def names(self):
        ''' Returns the names of the elements in the twiss table '''
        return self._names

class TfsSummary(LookupDict):
    ''' A class to hold the summary table of a twiss with lowercase keys '''
    pass

def tfs(inputfile):
    '''
     Returns table and summary information
     as LookUp dictionaries. These extend on normal
     dictionary syntax. We recommend using this function
     for reading tfs files.

     :param string inputfile: tfs file, full path
    '''
    table,params=tfsDict(inputfile)
    return TfsTable(table), TfsSummary(params)

def tfsDict(inputfile):
    '''
    .. py:function:: tfsDict(inputfile)

    Read a tfs table and returns table/summary info

    The function takes in a tfs file. It will add
    all parameters into one dictionary, and the table
    into another dictionary.

    :param string inputfile: tfs file, full path
    :raises ValueError: In case file path is not found
    :rtype: tuple containing dictionaries (tfs table , summary)

    '''
    params={}
    if not os.path.isfile(inputfile):
        if os.path.isfile(inputfile+'.tfs'):
            inputfile+='.tfs'
        elif os.path.isfile(inputfile+'.TFS'):
            inputfile+='.TFS'
        else:
            raise ValueError("ERROR: "+inputfile+" is not a valid file path")
    f=file(inputfile,'r')
    l=f.readline()
    while(l):
        if l.strip()[0]=='@':
            _addParameter(params,l)
        if l.strip()[0]=='*': # beginning of vector list...
            names=l.split()[1:]
            table=_read_table(f,names)
        l=f.readline()
    return table, params

##
# Add parameter to object
#
# Any line starting with an @ is a parameter.
# If that is found, this function should be called and given the line
#
# @param line The line from the file that should be added
def _addParameter(params,line):
    lname=line.split()[1].lower()
    if line.split()[2]=='%le':
        params[lname]=float(line.split()[3])
    if line.split()[2][-1]=='s':
        params[lname]=line.split('"')[1]
    if line.split()[2]=='%d':
        params[lname]=int(line.split()[3])

##
# Reads in a table in tfs format.
# Input the file stream at the location
# where the names of the columns have just been read.
def _read_table(fstream,names):
    l=fstream.readline()
    types=[]
    table={}
    for n in names:
        table[n.lower()]=[]
    while(l):
        if l.strip()[0]=='$':
            types=l.split()[1:]
        else:
            for n,el in zip(names,l.split()):
                table[n.lower()].append(el)
        l=fstream.readline()
    for n,typ in zip(names,types):
        if typ=='%le':
            table[n.lower()]=numpy.array(table[n.lower()],dtype=float)
        elif typ=='%d':
            table[n.lower()]=numpy.array(table[n.lower()],dtype=int)
        elif typ=='%s':
            for k in xrange(len(table[n.lower()])):
                table[n.lower()][k]=table[n.lower()][k].split('"')[1]
    return table

