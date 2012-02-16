#!/usr/bin/env python

import tpg
from abpyutils import container

elems=container()
currseq=None;

def mkelem(d):
  madname,args=d
  name=madname.replace('.','_')
  if not hasattr(elems,name):
    if args and args[0][0]=='parent':
      parent=args[0][1]
      if hasattr(elems,parent):
        parent=getattr(elems,parent)
      else:
        parent=mkelem([parent,[]])
      class cl(parent): pass
      parent._ref.append(cl)
      args.pop(0)
    else:
      class cl(container): pass
    if parent=='sequence': currseq=cl
  else:
    cl=getattr(elems,name)
  for i,j in args:
    if i=='from': i='From'
    setattr(cl,i,j)
  cl.madname=madname
  cl.__name__=madname
  cl._ref=[]
  if name=='endsequence': currseq=None
  if currseq: cl.sequence=currseq
  setattr(elems,name,cl)
  return cl

vars=container()

def mkvar(d):
  madname,e=d
  name=madname.replace('.','_')
  setattr(vars,name,e)

class madx(tpg.VerboseParser):
    r"""
    set lexer_multiline = True
    separator blank     '[ \t\r\n]+'  ;
    separator comment1  '![^\n]*\n'   ;
    separator comment2  '//[^\n]*\n'  ;
    separator comment3  '/\*.*/\*'    ;
    token reserved  '(real const)|(shared)';
    #token name      '[_a-zA-Z\.][_a-zA-Z0-9\.]+';
    token mathexpr  '[_/a-zA-Z0-9\.\-\+\*\^\(\)]+';
    token eq      '(:=|=)';
    token eqp     ':';
    token comma   ',';
    START    -> linelist;
    linelist -> (line ) (line )* ;
    line    -> ( vardef/l $mkvar(l)$  | cmddef/l $mkelem(l)$ ) ';';
    vardef/v  -> (reserved)* mathexpr/n eq/q expr/e $v=[n,e]$ ;
    expr/m    -> mathexpr/m |
                 ('{' mathexpr/p $m=[p]$ )(comma mathexpr/p $m.append(p)$ )* '}';
    cmddef/c  -> (reserved)*
                 ( mathexpr/n $c=[n,[]]$)
                 ( parent/p   $c[1].append(['parent',p])$ )?
                 ( varlist/v   $c[1].extend(v)$ )?;
    parent/p  -> eqp mathexpr/p ;
    cmdopt/n  -> (vardef/n | mathexpr/m $n=[m,True]$);
    varlist/v -> ( comma cmdopt/n $v=[n]$ ) ( comma cmdopt/n $v.append(n)$ )* ;
    """
    verbose=0

seq = madx()
seq(open('V6.5.seq').read())
