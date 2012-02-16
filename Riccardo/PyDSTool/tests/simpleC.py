# simpleC.py
#
# Example of exporting a C definition of a vector field and helper
# functions with the minimum of set up and fuss!
# Example for a Hodgkin-Huxley cell

from PyDSTool import *

vfn_str = '(I-ionic(v,m,h,n))/C'
mfn_str = 'ma(v)*(1-m)-mb(v)*m'
nfn_str = 'na(v)*(1-n)-nb(v)*n'
hfn_str = 'ha(v)*(1-h)-hb(v)*h'

auxdict = {'ionic': (['vv', 'mm', 'hh', 'nn'],
                        'gna*mm*mm*mm*hh*(vv-vna) + gk*nn*nn*nn*nn*(vv-vk) + gl*(vv-vl)'),
           'ma': (['v'], '0.32*(v+54)/(1-exp(-(v+54)/4))'),
           'mb': (['v'], '0.28*(v+27)/(exp((v+27)/5)-1)'),
           'ha': (['v'], '.128*exp(-(50+v)/18)'),
           'hb': (['v'], '4/(1+exp(-(v+27)/5))'),
           'na': (['v'], '.032*(v+52)/(1-exp(-(v+52)/5))'),
           'nb': (['v'], '.5*exp(-(57+v)/40)')}

DSargs = {}
DSargs['name'] = 'my_name'
DSargs['varspecs'] = {'v': vfn_str, 'm': mfn_str,
                         'h': hfn_str, 'n': nfn_str
                         }
#DSargs['pars'] = par_args
DSargs['fnspecs'] = auxdict

# don't need the Cgen object itself, but this call creates the C file
# 'my_name_vf.c' in a local dopri853_temp/ directory.
Cgen = Generator.Dopri_ODEsystem(DSargs)
