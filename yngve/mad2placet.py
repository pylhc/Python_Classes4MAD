import numpy as np

import os,sys,shutil

'''
.. module:: yngve.mad2placet

Holds various functions for converting things in Mad-X into 
Placet input.

.. moduleauthor:: Yngve Inntjore Levinsen <Yngve.Inntjore.Levinsen@cern.ch>

'''

def errors(errtable,order=8,outfile='placet_new_errors.tcl'):
    '''
     Writes a placet multipole file for
     all elements in the errtable
     Typically add this in the centre of the
     main element...

     :param string errtable: Input file, tfs table, full path
     :param int order: The highest multipole order to be added
     :param string outfile: [optional] The output file name
    '''
    fout=file(outfile,'w')
    first=True
    in_cmd='Multipole -name "%(name)s" -length %(l)f -synrad $mult_synrad -type %(o)i -strength [expr %(k)e*$e0]\n'
    in_cmd_s=in_cmd.strip()+' -tilt %(tilt)e\n'
    def _write_multipoles(i,o,skew=False,ismain=False):
        '''
         Input: index, order, skew [bool]
        '''
        _o=str(o)
        o+=1
        l=0
        if skew:
            name='MULTS'+_o+errtable.name[i]
            k=getattr(errtable,'k'+_o+'sl')[i]
            tilt=np.pi/2/o
            fout.write(in_cmd_s % locals())
        else:
            name='MULTN'+_o+errtable.name[i]
            k=getattr(errtable,'k'+_o+'l')[i]
            fout.write(in_cmd % locals())
    for i in xrange(len(errtable.name)):
        fout.write('# Element name in Mad-X: '+errtable.name[i]+'\n')
        for o in xrange(order):
            if getattr(errtable,'k'+str(o)+'l')[i]!=0:
                _write_multipoles(i,o)
            if getattr(errtable,'k'+str(o)+'sl')[i]!=0:
                _write_multipoles(i,o,skew=True)
