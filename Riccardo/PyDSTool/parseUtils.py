#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      13 September 2006
#$Revision: 1.1.4 $
"""
    Parser utilities.

    Robert Clewley, September 2005.
    
    Includes AST code by Pearu Peterson and Ryan Gutenkunst.
"""

# ------------------------------------------------------------------------
# Version History
#
# Version 1.1.4
#   Added 'mod' function to list of protected names.
#   Fixed obscure bug in addArgToCalls (see code comments).
#   Fixed obscure bug in parserObject handling of numeric literals
#     involving exponent notation (names ending in 'e' prior to a '+' or
#     '-' were being treated as a number of the form <X>e+<EXP>).
#
# Version 1.1.3
#   Added 'getbound' in-built auxiliary function.
#
# Version 1.1.2
#   Added count_sep function.
#   Fixed simplify function to properly process symbolic division.
#   Usage of temp_macro_names function in convertPowers to avoid syntax
#     exception being raised when this function is used on expressions
#     containing the 'for' or 'if' pydstool macros. These use a different
#     syntax to regular python.
#
# Version 1.1.1
#   Conversion of powers in symbolic expressions can now be done between
#     any of the syntactic forms (using the 'pow', '**' or '^' operators).
#   Added function mapNames for mapping names for many types of object --
#     especially useful for mapping hierarchical names -> FuncSpec-
#     compatible names (using '_' separator instead of '.').
#   Updated symbolMapClass so that it can be used more like a dictionary.
#
# Version 1.1
#   Incorporated Pearu Peterson's 'symbolic' legacy code, with further
#     extensions and fixes.
#   Added ensurebare() function to allow dofun() to simplify expressions that
#     are already secured inside a function call and do not need additional braces.
#   Added simplify function to run appropriate 'do' functions on original
#     *arguments* to math function calls, which were previously left alone
#     inside those calls.
#   Added scopes to eval statements to prevent argument names clashing with
#     names appearing in the eval'd string.
#   Changed str() to repr() in trysimple to increase floating point accuracy.
#   Minor bug fixes in original code.
#
# Version 1.0.4
#   Support for builtin aux functions / macros "if" and "for" in eval().
#
# Version 1.0.3
#   Minor bug fixes, including to isCompound to catch expressions involving
#     function calls.
#
# Version 1.0.2
#   Added preserveSpace option to parser.
#
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------

# IMPORTS
from __future__ import division
from errors import *
from common import *
import re
import math, random
from scipy import alltrue, sometrue, any, all
from copy import copy
import parser, symbol, token

# --------------------------------------------------------------------------
# GLOBAL BEHAVIOUR CONTROLS -- leave both as True for use with PyDSTool.
# switch to change output from using 'x**y' (DO_POW=False)
# to 'pow(x,y)' (DO_POW=TRUE)
DO_POW=True
# switch to activate treatment of integer constants as decimals in expressions
# and simplifications of expressions.
DO_DEC=True

# ------------------------------------------------------------------------
### Protected names

protected_scipynames = ['sign', 'mod']
protected_mathnames = filter(lambda s: not s.startswith('__'), \
                                 dir(math))
protected_randomnames = filter(lambda s: not s.startswith('_'), \
                                 dir(random))  # yes, just single _
# We add internal default auxiliary function names for use by
# functional specifications.
builtin_auxnames = ['globalindepvar', 'initcond', 'heav', 'if',
                    'getindex', 'getbound']

protected_macronames = ['for', 'if']

reserved_keywords = ['and', 'del', 'for', 'is', 'raise', 'assert', 'elif',
                'from', 'lambda', 'return', 'break', 'else', 'global',
                'not', 'try', 'class', 'except', 'if', 'or', 'while',
                'continue', 'exec', 'import', 'pass', 'yield', 'def',
                'finally', 'in', 'print', 'as', 'None']

# 'abs' is defined in python core, so doesn't appear in math
protected_allnames = protected_mathnames + protected_scipynames \
                    + protected_randomnames \
                    + builtin_auxnames + protected_macronames \
                    + ['abs', 'pow']

# signature lengths for builtin auxiliary functions and macros, for
# use by ModelSpec in eval() method (to create correct-signatured temporary
# functions).
builtinFnSigInfo = {'globalindepvar': 1, 'initcond': 1, 'heav': 1, 'getindex': 1,
                    'if': 3, 'for': 4, 'getbound': 2}

# ------------------------------------------------------------------------

# EXPORTS
_functions = ['readArgs', 'findEndBrace', 'makeParList', 'joinStrs',
                'parseMatrixStrToDictStr', 'joinAsStrs', 'replaceSep',
                'wrapArgInCall', 'addArgToCalls', 'findNumTailPos',
                'isToken', 'isNameToken', 'isNumericToken', 'count_sep',
                'isHierarchicalName', 'replaceSepInv', 'replaceSepListInv',
                'replaceSepList', 'convertPowers', 'mapNames']

_objects = ['protected_auxnamesDB', 'protected_allnames', 'protected_macronames',
            'protected_mathnames', 'protected_randomnames', 'builtin_auxnames',
            'protected_scipynames', 'builtinFnSigInfo']

_classes = ['symbolMapClass', 'parserObject', 'auxfnDBclass']

_constants = ['name_chars_RE', 'num_chars', 'ZEROS', 'ONES', 'NAMESEP',
              '_indentstr']

_symbfuncs = ['simplify', 'simplify_str', 'ensurebare', 'ensureparen',
             'trysimple', 'ensuredecimalconst', 'doneg', 'dosub', 'doadd',
             'dodiv', 'domul', 'dopower', 'splitastLR', 'ast2string', 'string2ast',
             'sym2name', 'ast2shortlist', 'splitargs', 'mapPowStr',
             'toPowSyntax', 'ensureparen_div']

_symbconsts = ['syms']

__all__ = _functions + _classes + _objects + _constants + _symbfuncs + _symbconsts

#-----------------------------------------------------------------------------

## constants for parsing
name_chars_RE = re.compile('\w')
alphabet_chars_RE = re.compile('[a-zA-Z0-9]')   # without the '_'
num_chars = map(lambda i: str(i), range(10))

if DO_POW:
    POW_STR = 'pow(%s,%s)'
else:
    POW_STR = '%s**%s'

if DO_DEC:
    ZEROS = ['0','0.0','0.','(0)','(0.0)','(0.)']
    ONES = ['1','1.0','1.','(1)','(1.0)','(1.)']
    TENS = ['10', '10.0', '10.','(10)','(10.0)','(10.)']
else:
    ZEROS = ['0']
    ONES = ['1']
    TENS = ['10']

# separator for compound, hierarchical names
NAMESEP = '.'

# for python indentation in automatically generated code, use 4 spaces
_indentstr = "    "


# ----------------------------------------------------------------------------
# This section: code by Pearu Peterson, adapted by Ryan Gutenkunst
#   and Robert Clewley.

syms=token.tok_name
for s in symbol.sym_name.keys():
    syms[s]=symbol.sym_name[s]


def mapPowStr(t, p='**'):
    """Input an expression of the form 'pow(x,y)'. Outputs an expression of the
    form x**y, x^y, or pow(x,y) and where x and y have also been processed to the
    target power syntax.
    Written by R. Clewley"""
    ll = splitargs(ast2string(t[2])[1:-1])
    if p=='**':
        lpart = dopower(ensureparen(ast2string(toDoubleStarSyntax(string2ast(ll[0]))),1),
                       ensureparen(ast2string(toDoubleStarSyntax(string2ast(ll[1]))),1),
                       '%s**%s')
        if len(t) > 3:
            res = ensureparen(ast2string(toDoubleStarSyntax(['power',string2ast('__LPART__')]+t[3:])))
            return res.replace('__LPART__', lpart)
        else:
            return ensureparen(lpart,1)
    elif p=='^':
        lpart = dopower(ensureparen(ast2string(toCircumflexSyntax(string2ast(ll[0]))),1),
                       ensureparen(ast2string(toCircumflexSyntax(string2ast(ll[1]))),1),
                       '%s^%s')
        if len(t) > 3:
            res = ensureparen(ast2string(toCircumflexSyntax(['power',string2ast('__LPART__')]+t[3:])),1)
            return res.replace('__LPART__', lpart)
        else:
            return ensureparen(lpart,1)
    elif p=='pow':
        lpart = dopower(ensurebare(ast2string(toPowSyntax(string2ast(ll[0])))),
                       ensurebare(ast2string(toPowSyntax(string2ast(ll[1])))),
                       'pow(%s,%s)')
        if len(t) > 3:
            res = ensurebare(ast2string(toPowSyntax(['power',string2ast('__LPART__')]+t[3:])))
            return res.replace('__LPART__', lpart)
        else:
            return ensureparen(lpart,1)
    else:
        raise ValueError, "Invalid power operator"
        

def toCircumflexSyntax(t):
    # R. Clewley
    if isinstance(t[0], str):
        if t[0] == 'power':
            if t[2][0] == 'DOUBLESTAR':
                return string2ast(ensureparen(dopower(ast2string(toCircumflexSyntax(t[1])),
                                    ast2string(toCircumflexSyntax(t[3])),
                         '%s^%s'),1))
            if t[1] == ['NAME', 'pow']:
                return string2ast(ensureparen(mapPowStr(t,'^'),1))
    o = []
    for i in t:
        if isinstance(i,list):
            if type(i[0]) == str and i[0].islower():
                o.append(toCircumflexSyntax(i))
            else:
                o.append(i)
        else:
            o.append(i)
    return o


def toDoubleStarSyntax(t):
    # R. Clewley
    if isinstance(t[0], str):
        if t[0] == 'xor_expr' and t[2][0]=='CIRCUMFLEX':
            # ^ syntax has complex binding rules in python parser's AST!
            # trick - easy to convert to ** first. then, using a bit of a hack
            # convert to string and back to AST so that proper AST for **
            # is formed.
            tc = copy(t)
            tc[0] = 'power'
            tc[2] = ['DOUBLESTAR', '**']
            return toDoubleStarSyntax(string2ast(ast2string(tc)))   # yes, i mean this
        if t[0] == 'power' and t[1] == ['NAME', 'pow']:
            return string2ast(ensureparen(mapPowStr(t,'**'),1))
    o = []
    for i in t:
        if isinstance(i,list):
            if type(i[0]) == str and i[0].islower():
                o.append(toDoubleStarSyntax(i))
            else:
                o.append(i)
        else:
            o.append(i)
    return o


def toPowSyntax(t):
    # R. Clewley
    if isinstance(t[0],str):
        if t[0] == 'power':
            try:
                if t[2][0]=='DOUBLESTAR':
                    try:
                        return string2ast(dopower(ensurebare(ast2string(toPowSyntax(t[1]))),
                                            ensurebare(ast2string(toPowSyntax(t[3]))),
                                 'pow(%s,%s)'))
                    except IndexError:
                        # there's a pow statement already here, not a **
                        # so ignore
                        return t
                elif t[1][1] == 'pow':
                    return string2ast(ensureparen(mapPowStr(t,'pow'),1))
                elif len(t)>3 and t[3][0]=='DOUBLESTAR':
                    try:
                        return string2ast(dopower(ensurebare(ast2string(toPowSyntax(t[1:3]))),
                                            ensurebare(ast2string(toPowSyntax(t[4]))),
                                 'pow(%s,%s)'))
                    except IndexError:
                        # there's a pow statement already here, not a **
                        # so ignore
                        return t
            except:
                print t
                print ast2string(t)
                raise
        elif t[0] == 'xor_expr' and t[2][0]=='CIRCUMFLEX':
            # ^ syntax has complex binding rules in python parser's AST!
            # trick - easy to convert to ** first. then, using a bit of a hack
            # convert to string and back to AST so that proper AST for **
            # is formed.
            tc = copy(t)
            tc[0] = 'power'
            tc[2] = ['DOUBLESTAR', '**']
            return toPowSyntax(string2ast(ast2string(tc)))   # yes, i mean this
    o = []
    for i in t:
        if isinstance(i,list):
            if type(i[0]) == str and i[0].islower():
                o.append(toPowSyntax(i))
            else:
                o.append(i)
        else:
            o.append(i)
    return o


def temp_macro_names(s):
    t = s
    for m in reserved_keywords:
        t = t.replace(m, '__'+m+'__')
    return t

def temp_macro_names_inv(s):
    t = s
    for m in reserved_keywords:
        t = t.replace('__'+m+'__', m)
    return t
    

def convertPowers(s, target="pow"):
    """convertPowers takes a string argument and maps all occurrences
    of power sub-expressions into the chosen target syntax. That option
    is one of "**", "^", or "pow"."""
    # temp_macro_names switches python reserved keywords to adapted names
    # that won't make the python parser used in string2ast raise a syntax error
    # ... but only invoke this relatively expensive function in the unlikely
    # event a clash occurs.
    # R. Clewley
    if target=="**":
        try:
            return ast2string(toDoubleStarSyntax(string2ast(s)))
        except SyntaxError:
            s = temp_macro_names(s)
            return temp_macro_names_inv(ast2string(toDoubleStarSyntax(string2ast(s))))
    elif target=="^":
        try:
            return ast2string(ensureints(toCircumflexSyntax(string2ast(s))))
        except SyntaxError:
            s = temp_macro_names(s)
            return temp_macro_names_inv(ast2string(ensureints(toCircumflexSyntax(string2ast(s)))))
    elif target=="pow":
        try:
            return ast2string(toPowSyntax(string2ast(s)))
        except SyntaxError:
            s = temp_macro_names(s)
            return temp_macro_names_inv(ast2string(toPowSyntax(string2ast(s))))
    else:
        raise ValueError, "Invalid target syntax"


def ensureints(t):
    """Ensure that any floating point constants appearing in t that are
    round numbers get converted to integer representations."""
    if type(t)==str:
        return ast2string(ensureints(string2ast(t)))
    o = []
    for i in t:
        if type(i) == list:
            if type(i[0]) == str and i[0].islower():
                o.append(ensureints(i))
            elif i[0]=='NUMBER':
                # CAPS for constants
                o.append(string2ast(trysimple(ast2string(i))))
            else:
                o.append(i)
                
        else:
            o.append(i)
    return o
                

def splitargs(da, lbraces=['('], rbraces=[')']):
    """Function to split string-delimited arguments in a string without
    being fooled by those that occur in function calls.
    Written by Pearu Peterson. Adapted by Rob Clewley to accept different
    braces."""
    if all([da.find(lbrace)<0 for lbrace in lbraces]):
        return da.split(',')
    ll=[];o='';ii=0
    for i in da:
        if i==',' and ii==0:
            ll.append(o)
            o=''
        else:
            if i in lbraces: ii=ii+1
            if i in rbraces: ii=ii-1
            o=o+i
    ll.append(o)
    return ll

def ast2shortlist(t):
    if type(t) is parser.ASTType: return ast2shortlist(t.tolist())
    if not isinstance(t, list): return t
    if t[1] == '': return None
    if not isinstance(t[1], list): return t
    if len(t) == 2 and isinstance(t[1], list):
        return ast2shortlist(t[1])
    o=[]
    for tt in map(ast2shortlist, t[1:]):
        if tt is not None:
            o.append(tt)
    if len(o)==1: return o[0]
    return [t[0]]+o

def sym2name(t):
    if type(t) is parser.ASTType: return sym2name(t.tolist())
    if not isinstance(t, list): return t
    return [syms[t[0]]]+map(sym2name,t[1:])

def string2ast(t):
    return sym2name(ast2shortlist(parser.expr(t)))

def ast2string(t):
    #if isinstance(t, str): return t
    if type(t) is parser.ASTType: return ast2string(t.tolist())
    if not isinstance(t, list): return None
    if not isinstance(t[1], list): return t[1]
    o=''
    for tt in map(ast2string,t):
        if isinstance(tt, str):
            o=o+tt
    return o

def splitastLR(t):
    lft=t[1]
    rt=t[3:]
    if len(rt)>1:
        rt=[t[0]]+rt
    else:
        rt=rt[0]
    return lft,rt

def dopower(l,r,pow_str=POW_STR):
    if r in ZEROS: return '1'
    if l in ZEROS: return '0'
    if l in ONES: return '1'
    if r in ONES: return (l)
    if pow_str=='%s**%s':
        return (trysimple('%s**%s'%(ensureparen(l),ensureparen(r))))
    elif pow_str == '%s^%s':
        return (trysimple('%s^%s'%(ensuredecimalconst(ensureparen(l)),ensuredecimalconst(ensureparen(r)))))
    elif pow_str == 'pow(%s,%s)':
        return (trysimple('pow(%s,%s)'%(ensurebare(l),ensurebare(r))))
    else:
        raise ValueError, "Invalid target power syntax"

def domul(l,r):
    if l in ZEROS or r in ZEROS: return '0'
    if l in ONES: return (r)
    if r in ONES: return (l)
    return (trysimple('%s*%s'%(ensureparen(l),ensureparen(r))))

def dodiv(l,r):
    if r in ZEROS: raise ValueError, "Division by zero in expression"
    if l in ZEROS: return '0'
    if r in ONES: return (l)
    if r in ['-'+o for o in ONES]: return (doneg(l))
    if r==l: return '1'
    return (trysimple('%s/%s'%(ensureparen(ensuredecimalconst(l)),
                               ensureparen(ensuredecimalconst(r),1))))

def doadd(l,r):
    if l in ZEROS and r in ZEROS: return '0'
    if l in ZEROS: return r
    if r in ZEROS: return l
    if l==r: return trysimple(domul('2',l))
    if r[0]=='-': return trysimple('%s%s'%(l,r))
    return ensureparen(trysimple('%s+%s'%(l,r)))

def dosub(l,r):
    if l in ZEROS and r in ZEROS: return '0'
    if l in ZEROS: return doneg(r)
    if r in ZEROS: return l
    if l==r: return '0'
    if r[0]=='-': return ensureparen(trysimple('%s-%s'%(l,r)))
    return ensureparen(trysimple('%s-%s'%(l,r)))

def doneg(l):
    if l in ZEROS: return '0'
    if l[0]=='-': return l[1:]
    return trysimple('-%s'%l)

def ensuredecimalconst(t):
    return trysimple(t,do_decimal=True)

def trysimple(t,do_decimal=False):
    try:
        t_e = eval(t, {}, {})
        if do_decimal and type(t_e) == int and t_e != 0 and DO_DEC:
            t = repr(t_e)+"."
        elif type(t_e) == float and int(t_e)==t_e:
            t = repr(int(t_e))
        else:
            t = repr(t_e)
    except:
        pass
    return t

def ensureparen_div(tt):
    if tt[0] == 'term' and tt[2][0] == 'SLASH' and len(tt[3:])>1:
        return [tt[0], string2ast(ensureparen(ast2string(tt[1])+'/'+ast2string(tt[3]),1))] + tt[4:]
    else:
        return tt

def ensureparen(t,flag=0):
    t=trysimple(t)
    if t[0]=='-':
##        print "0: ", t, "->", '(%s)'%t
        return '(%s)'%t
    tt=string2ast(t)
    if tt[0]=='arith_expr':
##        print "1: ", t, "->", '(%s)'%t
        return '(%s)'%t
    if flag>0:
        if tt[0] == 'term':
            return '(%s)'%t
        elif tt[0] == 'power':
            if tt[1] == ['NAME', 'pow']:
                return t
            else:
                # ** case
                for x in tt[1:]:
                    if x[0] == 'arith_expr':
                        return '(%s)'%t
        elif tt[0] == 'xor_expr':    # added xor_expr for ^ powers
            if len(tt)>3:
                for x in tt[1:]:
                    if x[0] == 'arith_expr':
                        return '(%s)'%t
            else:
                return t
##        if t[0]=='(' and t[-1]==')' and t[1:-1].find('(') < t[1:-1].find(')'):
##            # e.g. (1+x) doesn't need another set of braces
##            print "2: ", t, "->", t
##            return t
##        else:
##            # e.g. (1+x)-(3-y) does need braces
##            print "3: ", t, "->", '(%s)'%t
##            return '(%s)'%t
    return t

def ensurebare(t):
    """Ensure no braces in string expression (where possible).
    Written by Robert Clewley"""
    t=trysimple(t)
    try:
        if t[0]=='(' and t[-1]==')':
            if t[1:-1].find('(') < t[1:-1].find(')'):
                return t[1:-1]
            else:
                # false positive, e.g. (1+x)-(3-y)
                return t
        else:
            return t
    except IndexError:
        return t

def simplify_str(s):
    """String output version of simplify"""
    t=string2ast(s)
    if isinstance(t, list) and len(t) == 1:
        return ast2string(t[0])
    if t[0] in ['NUMBER', 'NAME']: return s
    if t[0]=='factor':
        return doneg(simplify_str(ast2string(t)[1:]))
    if t[0]=='arith_expr':
        o='0'
        for i in range(1,len(t),2):
            if t[i-1][1]=='-':
                o=dosub(o,ast2string(simplify(t[i])))
            else:
                o=doadd(o,ast2string(simplify(t[i])))
        return o
    if t[0]=='term':
        ttemp = ensureparen_div(t)
        lft,rt=map(ast2string,splitastLR(ttemp))
        if ttemp[2][0]=='STAR':
            return domul(ast2string(simplify(string2ast(lft))),ast2string(simplify(string2ast(rt))))
        if ttemp[2][0]=='SLASH':
            return dodiv(ast2string(simplify(string2ast(lft))),ast2string(simplify(string2ast(rt))))
    if t[0]=='power': # covers math functions like sin and cos as well
        if t[2][0]=='trailer':
            if len(t)>3:
                return ast2string([t[0],string2ast('(%s)'%ast2string(simplify(t[:3]))),simplify(t[3:])])
            elif len(t)==3 and t[1][1]=='pow':
                # 'pow' syntax case
                return ast2string(toPowSyntax(simplify(string2ast(mapPowStr(t)))))
            # else ignore and carry on
        ts=[];o=[];
        for i in t[1:]:
            if i[0]=='DOUBLESTAR':
                if len(o)==1: o=o[0]
                ts.append(o);
                o=[]
            else: o.append(simplify(i))
        if len(o)==1: o=o[0]
        ts.append(o)
        if t[2][0]=='DOUBLESTAR':
            st,lft,rt=map(ast2string,[t,ts[0],ts[1]])
            return dopower(ast2string(simplify(string2ast(lft))),ast2string(simplify(string2ast(rt))))
        if t[2][0]=='trailer':
            return ast2string(simplify(ts[0]))
    if t[0] in ['arglist','testlist']:
        o=[]
        for i in ast2string(t).split(','):
            o.append(ast2string(simplify(string2ast(i))))
        return ','.join(o)
    if t[0]=='atom':
        return ensureparen(ast2string(simplify(t[2:-1])))
    if t[1][0]=='trailer': # t=[[NAME,f],[trailer,[(],[ll],[)]]]
        # just simplify arguments to functions
        return ast2string([['NAME',t[0][1]], ['trailer',['LPAR','('],simplify(t[1][2]),['RPAR',')']]])
    return s


def simplify(t):
    """Attempt to simplify symbolic expression string.
    Adapted by R. Clewley from original DiffStr() code by Pearu Peterson and Ryan Gutenkunst.
    
    Essentially the same format as DiffStr() except we jus moves down the
    syntax tree calling appropriate 'do' functions on the parts to
    simplify them."""
    if isinstance(t, list) and len(t) == 1:
        return t[0]
    if t[0] in ['NUMBER', 'NAME']: return t
    if t[0]=='factor':
        return string2ast(doneg(simplify_str(ast2string(t)[1:])))
    if t[0]=='arith_expr':
        o='0'
        for i in range(1,len(t),2):
            if t[i-1][1]=='-':
                o=dosub(o,ast2string(simplify(t[i])))
            else:
                o=doadd(o,ast2string(simplify(t[i])))
        return string2ast(o)
    if t[0]=='term':
        ttemp = ensureparen_div(t)
        lft,rt=map(ast2string,splitastLR(ttemp))
        if ttemp[2][0]=='STAR':
            return string2ast(domul(ast2string(simplify(string2ast(lft))),ast2string(simplify(string2ast(rt)))))
        if ttemp[2][0]=='SLASH':
            return string2ast(dodiv(ast2string(simplify(string2ast(lft))),ast2string(simplify(string2ast(rt)))))
    if t[0]=='xor_expr' and t[2]=='CIRCUMFLEX': # covers alternative power syntax
        alt = copy(t)
        alt[0] = 'power'; alt[2]=='DOUBLESTAR'
        return toCircumflexSyntax(simplify(alt))
    if t[0]=='power': # covers math functions like sin and cos as well
        if t[2][0]=='trailer':
            if len(t)>3:
                return [t[0],string2ast('(%s)'%ast2string(simplify(t[:3]))),simplify(t[3:])]
            elif len(t)==3 and t[1][1]=='pow':
                # 'pow' syntax case
                return toPowSyntax(simplify(string2ast(mapPowStr(t))))
            # else ignore and carry on
        ts=[];o=[];
        for i in t[1:]:
            if i[0]=='DOUBLESTAR':
                if len(o)==1: o=o[0]
                ts.append(o);
                o=[]
            else: o.append(simplify(i))
        if len(o)==1: o=o[0]
        ts.append(o)
        if t[2][0]=='DOUBLESTAR':
            st,lft,rt=map(ast2string,[t,ts[0],ts[1]])
            return string2ast(dopower(ast2string(simplify(string2ast(lft))),ast2string(simplify(string2ast(rt)))))
        if t[2][0]=='trailer':
            return simplify(ts[0])
    if t[0] in ['arglist','testlist']:
        o=[]
        for i in ast2string(t).split(','):
            o.append(ast2string(simplify(string2ast(i))))
        return string2ast(','.join(o))
    if t[0]=='atom':
        return string2ast(ensureparen(ast2string(simplify(t[2:-1]))))
    if t[1][0]=='trailer': # t=[[NAME,f],[trailer,[(],[ll],[)]]]
        # just simplify arguments to functions
        return [['NAME',t[0][1]], ['trailer',['LPAR','('],simplify(t[1][2]),['RPAR',')']]]
    return t


# ----------------------------------------------------------------------------

class symbolMapClass(object):
    """Abstract class for hassle-free symbol re-mappings."""
    def __init__(self, symbolMap={}):
        if isinstance(symbolMap, symbolMapClass):
            self.lookupDict = copy(symbolMap.lookupDict)
        else:
            self.lookupDict = copy(symbolMap)

    def __call__(self, arg):
        if isinstance(arg, list):
            return [self.__getitem__(symbol) for symbol in arg]
        else:
            return self.__getitem__(arg)

    def __eq__(self, other):
        try:
            return self.lookupDict == other.lookupDict
        except AttributeError:
            return False

    def __neq__(self, other):
        return not self.__eq__(other)

    def __setitem__(self, symbol, mappedsymbol):
        self.lookupDict[symbol] = mappedsymbol

    def __getitem__(self, symbol):
        try:
            return self.lookupDict[symbol]
        except KeyError:
            return symbol

    def __contains__(self, symbol):
        return self.lookupDict.__contains__(symbol)

    def keys(self):
        return self.lookupDict.keys()

    def values(self):
        return self.lookupDict.values()

    def items(self):
        return self.lookupDict.items()

    def iterkeys(self):
        return self.lookupDict.iterkeys()

    def itervalues(self):
        return self.lookupDict.itervalues()

    def iteritems(self):
        return self.lookupDict.iteritems()

    def update(self, amap):
        try:
            # perhaps passed a symbolMapClass object
            self.lookupDict.update(amap.lookupDict)
        except AttributeError:
            # was passed a dict
            self.lookupDict.update(amap)

    def copy(self):
        return symbolMapClass(self)


    def __repr__(self):
        return "Symbol mapping"

    __str__ = __repr__

    def has_key(self, k):
        return k in self.lookupDict


class auxfnDBclass(object):
    """Auxiliary function database, for use by parsers."""
    def __init__(self):
        self.auxnames = {}

    def addAuxFn(self, auxfnName, parserObj):
        if parserObj not in self.auxnames:
            self.auxnames[parserObj] = auxfnName
        else:
            raise ValueError("Parser object " + parserObj.name + " already "
                             "exists in auxiliary function database")

    def __repr__(self):
        return "ModelSpec internal helper class: auxfnDBclass object"

    __str__ = __repr__

    def __call__(self, parserObj=None):
        if parserObj is None:
            # return all auxiliary functions known
            return self.auxnames.values()
        else:
            try:
                return [self.auxnames[parserObj]]
            except KeyError:
                return []

    def removeAuxFn(self, auxfnName):
        flagdelete = None
        for k, v in self.auxnames.iteritems():
            if v == auxfnName:
                flagdelete = k
                break
        if flagdelete is not None:
            del self.auxnames[k]

    def clear(self, parserObj):
        if parserObj in self.auxnames:
            del self.auxnames[parserObj]


    def clearall(self):
        self.auxnames = {}

# only need one of these per session
global protected_auxnamesDB
protected_auxnamesDB = auxfnDBclass()


class parserObject(object):
    """Alphanumeric symbol (pseudo-)parser for mathematical expressions.

    An AST is not properly implemented -- rather, we tokenize,
    identify free symbols, and apply a small number of syntactic rule checks.
    The target language parser is relied upon for full syntax checking.
    """

    def __init__(self, specStr, includeProtected=True,
                 treatMultiRefs=False, ignoreTokens=[],
                 preserveSpace=False):
        # freeSymbols does not include protected names, and so forth.
        # freeSymbols only contains unrecognized alphanumeric symbols.
        self.usedSymbols = []
        self.freeSymbols = []
        self.preserveSpace = preserveSpace
        # token by token list of the specStr
        self.tokenized = []
        if type(specStr) is str:
            self.specStr = specStr
        else:
            print "Found type", type(specStr), ": ", specStr
            raise TypeError, "specStr must be a string"
        self.treatMultiRefs = treatMultiRefs
        # record init options in case want to reset
        self.ignoreTokens = copy(ignoreTokens)
        self.includeProtected = includeProtected
        # process specStr with empty symbol map to create
        # self.freeSymbols and self.usedSymbols
        self.parse(ignoreTokens, None, includeProtected, reset=True)


    def isCompound(self, ops=['+', '-', '*', '/']):
        """Function to verify whether an expression is 'compound',
        in the sense that it has an operator at the root of its syntax
        parse tree (i.e. not inside braces)."""
        result = False
        nested = 0
        if len(self.tokenized) > 2:
            stage = 0
            for s in self.tokenized:
                if stage == 0:
                    if s in self.usedSymbols and nested == 0:
                        stage = 1
                    elif s == ')':
                        stage = 1
                        nested = max([0, nested-1])
                    elif s == '(':
                        nested += 1
                elif stage == 1:
                    if s in ops and nested == 0:
                        stage = 2
                    elif s == '(':
                        nested += 1
                    elif s == ')':
                        nested = max([0, nested-1])
                    elif nested == 0:
                        stage = 0
                elif stage == 2:
                    if s in self.usedSymbols and nested == 0:
                        stage = 3
                    elif s == '(':
                        stage = 3
                        nested += 1
                    elif s == ')':
                        nested = max([0, nested-1])
                    elif nested == 0:
                        stage = 0
                if stage == 3:
                    result = True
                    break
        return result

    def __call__(self, specialtoks=None, symbolMap=None, includeProtected=True):
        if specialtoks is None:
            if self.ignoreTokens is not None:
                specialtoks = self.ignoreTokens
            else:
                specialtoks = []
        if self.tokenized == []:
            return self.parse(specialtoks, symbolMap, includeProtected)
        else:
            if symbolMap is None:
                return "".join(self.tokenized)
            else:
                return "".join(symbolMap(self.tokenized))


    def parse(self, specialtoks, symbolMap=None, includeProtected=True,
              reset=False):
        if reset:
            self.usedSymbols = []
            self.freeSymbols = []
        if symbolMap is None:
            symbolMap = symbolMapClass()
        specialtokens = specialtoks + ['('] + self.usedSymbols
        if includeProtected:
            specialtokens.extend(protected_allnames)
            protected_auxnames = protected_auxnamesDB(self)
        else:
            protected_auxnames = []
        if self.treatMultiRefs:
            specialtokens.append('[')
        allnames = specialtokens + protected_auxnames
        specstr = self.specStr   # eases notation
        returnstr = ""
        if specstr == "":
            # Hack for empty specs to pass without losing their last characters
            specstr = " "
        elif specstr[-1] != ')':
            # temporary hack because strings not ending in ) lose their last
            # character!
            # Problem could be with line "if scount < speclen - 1:" below
            # ... should it be <= ?
            specstr += " "
        scount = 0
        speclen = len(specstr)
        # temp holders for used and free symbols
        used = copy(self.usedSymbols)
        free = copy(self.freeSymbols)
        tokenized = []
        # current token being built is: s
        # current character being processed is: stemp
        s = ''
        foundtoken = False
        while scount < speclen:
            stemp = specstr[scount]
            ## DEBUGGING PRINT STATEMENTS
##            print "\n*********************************************"
##            print specstr[:scount]
##            print stemp
##            print s
##            print tokenized
            scount += 1
            if name_chars_RE.match(stemp) is None:
                # then stemp is a non-alphanumeric char
                if s not in ['', ' ', '\n', '\t']:
                    # checking allnames catches var names etc. that are valid
                    # in auxiliary functions but are not special tokens
                    # and must be left alone
                    if s in allnames:
                        snew = symbolMap(s)
                        tokenized.append(snew)
                        if snew not in used:
                            used.append(snew)
                        returnstr += snew
                    else:
                        # isdecimal cases:
                        #  s = number followed by a '.'
                        #  s = number with 'e' followed by +/-
                        #      (just a number dealt with elsewhere)
                        isnumtok = isNumericToken(s)
                        issimpledec = stemp == '.' and isnumtok
                        isexpdec = s[-1] in ['e','E'] and s[0] not in ['e','E']\
                             and stemp in ['-','+'] and isnumtok
                        isdecimal = issimpledec or isexpdec
                        ishierarchicalname = stemp == '.' and isNameToken(s)
                        if isdecimal or ishierarchicalname:
                            # continue building token
                            s += stemp
                            continue
                        else:
                            # We have found a complete token
                            snew = symbolMap(s)
                            if s[0] not in num_chars + ['+','-'] \
                                    and snew not in free:
                                free.append(snew)
                            tokenized.append(snew)
                            if snew not in used:
                                used.append(snew)
                            returnstr += snew
                if stemp in ['+', '-']:
                    # may be start of a unary number
                    try:
                        next_stemp = specstr[scount]
                    except IndexError:
                        tokenized.append(stemp)
                        returnstr += stemp
                        s = ''
                        continue
                    if (tokenized==[] or tokenized[-1]=='(') and \
                       next_stemp in num_chars + ['.']:
                        # continue to build token
                        s += stemp
                        continue
                    elif len(tokenized)>0 and tokenized[-1] == '+':
                        # process double sign
                        if stemp == '-':
                            tokenized[-1] = '-'
                            returnstr = returnstr[:-1] + '-'
                        # else do nothing -- keep just one '+'
                        s = ''
                        continue
                    elif len(tokenized)>0 and tokenized[-1] == '-':
                        # process double sign
                        if stemp == '-':
                            tokenized[-1] = '+'
                            returnstr = returnstr[:-1] + '+'
                        # else do nothing -- keep just one '-'
                        s = ''
                        continue
                    else:
                        tokenized.append(stemp)
                        returnstr += stemp
                        s = ''
                        continue
                elif stemp in ['`', '!', '@', '#', '$', '{',
                               '}', "\\"]:
                    if stemp in specialtokens:
                        tokenized.append(stemp)
                        returnstr += stemp
                        s = ''
                        continue
                    else:
##                        if stemp == '^':
##                            raise ValueError('Symbol ^ is not allowed. '
##                                             'Please use the pow() call')
##                        else:
                        print "Problem with string '%s'"%specstr
                        raise ValueError('Symbol %s is illegal. '%stemp)
                elif stemp == '[':
                    # self.treatMultiRefs == False and '[' in specialtokens
                    # means it was probably in the ignoreToken list in __init__
                    if self.treatMultiRefs and len(tokenized)>0 \
                            and tokenized[-1].isalnum():
                        # then this is probably an actual multiRef
                        s = '['
                    elif '[' in specialtokens:
                        returnstr += '['
                        # s already to tokenized in this case
                        tokenized.append('[')   # was just '['
                        s = ''
                        continue
                    else:
                        raise ValueError("Syntax error: Square braces not to "
                                         "be used outside of multiple Quantity"
                                         " definitions and references")
                # only use next clause if want to process function call
                # arguments specially
#                elif stemp == '(':
#                    returnstr += s
#                    s = stemp
                else:
                    if stemp == "*":
                        if len(returnstr)>1 and returnstr[-1] == "*":
                            # check for ** case
                            if tokenized[-1] == '*':
                                tokenized[-1] = '**'
                            else:
                                tokenized.append('**')
                            s = ''
                            returnstr += stemp
                            continue   # avoids returnstr += stemp below
##                            if "**" in specialtokens:
##                                if tokenized[-1] == '*':
##                                    tokenized[-1] = '**'
##                                else:
##                                    tokenized.append('**')
##                                s = ''
##                                returnstr += "**"
##                                continue   # avoids returnstr += stemp below
##                            else:
##                                raise ValueError('Operator ** is not allowed. '
##                                       'Please use the pow() call')
                        else:
                            # just a single *
                            tokenized.append('*')
                    elif stemp in [" ","\t","\n"]:
                        if self.preserveSpace: tokenized.append(stemp)
                    else:
                        tokenized.append(stemp)
                    s = ''
                    returnstr += stemp
                    continue
            else:
                s += stemp
            if s in specialtokens:
                # only use next clause if want to process function call
                # arguments specially
#                if s == '(':
#                    foundtoken = False  # so that don't enter next if statement
#                    tokenized.append(s)
#                    returnstr += s
#                    s = ''
                if s == '[' and self.treatMultiRefs and len(tokenized)>0 \
                        and tokenized[-1].isalnum():
                    # then this is probably an actual multiRef ...
                    # only treat as multiRef if there's an alphanumeric
                    # token directly preceding '[', e.g. z[i,0,1]
                    # otherwise will catch array usage, e.g. [[x,y],[a,b]]
                    foundtoken = False    # don't go into next clause
                    # copy the square-bracketed clause to the output
                    # verbatim, and add whole thing as a token.
                    try:
                        rbpos = specstr[scount:].index(']')
                    except ValueError:
                        raise ValueError, "Mismatch [ and ] in spec"
                    # expr includes both brackets
                    expr = specstr[scount-1:scount+rbpos+1]
                    # find index name in this expression. this should
                    # be the only free symbol, otherwise syntax error.
                    temp = parserObject(expr[1:-1],
                                        includeProtected=False)
                    if len(temp.freeSymbols) == 1:
                        free.extend(temp.freeSymbols)
                        # not sure whether should add to usedSymbols
                        # for consistency with freeSymbols, or leave
                        # it out for consistency with tokenized
#                                used.append(temp.freeSymbols)
                    else:
                        raise ValueError("Invalid index clause in "
                                         "multiple quantity reference -- "
                                         "multiple index names used in [...]")
                    # start next parsing iteration after the ']'
                    scount += rbpos+1
                    # treat whole expression as one symbol
                    returnstr += expr
                    tokenized.append(expr)
                    used.append(expr)
                    s = ''
                else:
                    if scount < speclen - 1:
                        if name_chars_RE.match(specstr[scount]) is None:
                            foundtoken = True
                        else:
                            if s[-1] in ['e','E'] and s[0] not in ['e','E'] and \
                               name_chars_RE.match(specstr[scount]).group() \
                                       in num_chars+['-','+']:
                                # not expecting an arithmetic symbol or space
                                # ... we *are* expecting a numeric
                                foundtoken = True
                    else:
                        foundtoken = True
                if foundtoken:
                    if includeProtected:
                        if s == 'for':
                            # check next char is '('
                            if specstr[scount] != '(':
                                print "Next char found:", specstr[scount]
                                raise ValueError("Invalid 'for' macro syntax")
                            # find next ')' (for statement should contain no
                            # braces itself)
                            try:
                                rbpos = specstr[scount:].index(')')
                            except ValueError:
                                raise ValueError("Mismatch ( and ) in 'for' "
                                                 "macro")
                            # expr includes both brackets
                            expr = specstr[scount:scount+rbpos+1]
                            # find index name in this expression. this should
                            # be the only free symbol, otherwise syntax error.
                            temp = parserObject(expr, includeProtected=False)
                            macrotests = [len(temp.tokenized) == 7,
                               temp.tokenized[2] == temp.tokenized[4] == ',',
                               temp.tokenized[5] in ['+','*']]
                            if not all(macrotests):
                                print "specstr was: ", specstr
                                print "tokens: ", temp.tokenized
                                print "test results: ", macrotests
                                raise ValueError("Invalid sub-clause in "
                                                 "'for' macro")
                            # start next parsing iteration after the ')'
                            scount += rbpos+1
                            # keep contents of braces as one symbol
                            returnstr += s+expr
                            tokenized.extend([s,expr])
                            if s not in used:
                                used.append(s)
                            if expr not in used:
                                used.append(expr)
                        elif s == 'abs':
                            snew = symbolMap(s)
                            returnstr += snew
                            tokenized.append(snew)
                            if snew not in used:
                                used.append(snew)
                        elif s in protected_scipynames:
                            snew = symbolMap(s)
                            returnstr += snew
                            tokenized.append(snew)
                            if snew not in used:
                                used.append(snew)
                        elif s in protected_mathnames:
                            if s in ['e','E']:
                                # special case where e is either = exp(0)
                                # as a constant or it's an exponent in 1.0e-4
                                if len(returnstr)>0:
                                    if returnstr[-1] not in num_chars+['.']:
                                        snew = symbolMap(s.lower())
                                        returnstr += snew
                                        tokenized.append(snew)
                                        if snew not in used:
                                            used.append(snew)
                                    else:
                                        returnstr += s
                                else:
                                    snew = symbolMap(s.lower())
                                    returnstr += snew
                                    tokenized.append(snew)
                                    if snew not in used:
                                        used.append(snew)
                            else:
                                snew = symbolMap(s)
                                returnstr += snew
                                tokenized.append(snew)
                                if snew not in used:
                                    used.append(snew)
                        elif s in protected_randomnames:
                            snew = symbolMap(s)
                            returnstr += snew
                            tokenized.append(snew)
                            if snew not in used:
                                used.append(snew)
                        elif s in protected_auxnames:
                            snew = symbolMap(s)
                            returnstr += snew
                            tokenized.append(snew)
                        else:
                            # s is e.g. a declared argument to an aux fn but
                            # only want to ensure it is present. no action to
                            # take.
                            snew = symbolMap(s)
                            tokenized.append(snew)
                            returnstr += snew
                    else:
                        # not includeProtected names case, so just map
                        # symbol
                        snew = symbolMap(s)
                        tokenized.append(snew)
                        returnstr += snew
                    # reset for next iteration
                    s = ''
                    foundtoken = False
        # end of scount while loop
        # strip extraneous whitespace
        if reset:
            self.usedSymbols = used
            self.freeSymbols = free
            self.specStr = returnstr
            self.tokenized = tokenized
        return returnstr.strip()


#-----------------------------------------------------------------------------


def isToken(s, treatMultiRefs=False):
    # token must be alphanumeric string (with no punctuation, operators, etc.)
    if type(s) is not str:
        return False
    try:
        temp = parserObject(s, includeProtected=False,
                            treatMultiRefs=treatMultiRefs)
    except ValueError:
        return False
    if treatMultiRefs and s.find('[')>0:
        lenval = 2
    else:
        lenval = 1
    return not temp.isCompound() and len(temp.usedSymbols) == lenval \
            and len(temp.freeSymbols) == lenval \
            and len(temp.tokenized) == lenval


def isNameToken(s, treatMultiRefs=False):
    return isToken(s, treatMultiRefs=treatMultiRefs) \
        and s[0] not in num_chars and not (len(s)==1 and s=="_")


def isNumericToken(arg):
    # supports unary + / - at front, and checks for usage of exponentials
    # (using 'E' or 'e')
    s = arg.lower()
    try:
        if s[0] in ['+','-']:
            s_rest = s[1:]
        else:
            s_rest = s
    except IndexError:
        return False
    pts = s.count('.')
    exps = s.count('e')
    pm = s_rest.count('+') + s_rest.count('-')
    if pts > 1 or exps > 1 or pm > 1:
        return False
    if exps == 1:
        exp_pos = s.find('e')
        pre_exp = s[:exp_pos]
        # must be numbers before and after the 'e'
        if not any([n in num_chars for n in pre_exp]):
            return False
        if s[-1]=='e':
            # no chars after 'e'!
            return False
        if not any([n in num_chars for n in s[exp_pos:]]):
            return False
        # check that any additional +/- occurs directly after 'e'
        if pm == 1:
            pm_pos = max([s_rest.find('+'), s_rest.find('-')])
            if s_rest[pm_pos-1] != 'e':
                return False
            e_rest = s_rest[pm_pos+1:]   # safe due to previous check
        else:
            e_rest = s[exp_pos+1:]
        # only remaining chars in s after e and possible +/- are numbers
        if '.' in e_rest:
            return False
    # cannot use additional +/- if not using exponent
    if pm == 1 and exps == 0:
        return False
    return bool(all([n in num_chars + ['.', 'e', '+', '-'] for n in s_rest]))


def findNumTailPos(s):
    """Find position of numeric tail in alphanumeric string.

    e.g. findNumTailPos('abc678') = 3"""
    try:
        l = len(s)
        if l > 1:
            if s[-1] not in num_chars or s[0] in num_chars:
                raise ValueError("Argument must be an alphanumeric string "
                              "starting with a letter and ending in a number")
            for i in range(1, l+1):
                if s[-i] not in num_chars:
                    return l-i+1
        else:
            raise ValueError("Argument must be alphanumeric string starting "
                             "with a letter and ending in a number")
    except TypeError:
        raise ValueError("Argument must be alphanumeric string starting "
                             "with a letter and ending in a number")


def isHierarchicalName(s, sep=NAMESEP, treatMultiRefs=False):
    s_split = s.split(sep)
    return len(s_split) > 1 and bool(all([isNameToken(t, treatMultiRefs) \
                                         for t in s_split]))


def replaceSepList(speclist):
    return [replaceSep(spec) for spec in speclist]

def replaceSepListInv(speclist):
    return [replaceSepInv(spec) for spec in speclist]

def replaceSepInv(spec):
    """Invert default name separator replacement."""
    return replaceSep(spec, "_", NAMESEP)


def replaceSep(spec, sourcesep=NAMESEP, targetsep="_"):
    """Replace hierarchy separator character with another and return spec
    string. e.g. "." -> "_"
    Only replaces the character between name tokens, not between numbers."""

    for char in sourcesep:
        if alphabet_chars_RE.match(char) is not None:
            raise ValueError("Source separator must be non-alphanumeric")
    for char in targetsep:
        if alphabet_chars_RE.match(char) is not None:
            raise ValueError("Target separator must be non-alphanumeric")    
    if isinstance(spec, str):
        return replaceSepStr(spec, sourcesep, targetsep)
    else:
        # safe way to get string definition from either Variable or QuantSpec
        return replaceSepStr(str(spec()), sourcesep, targetsep)


def mapNames(themap, target):
    """Map names in <target> argument using the symbolMapClass
    object <themap>, returning a renamed version of the target."""
    try:
        themap.lookupDict
    except AttributeError:
        t = repr(type(themap))
        raise TypeError, \
          "Map argument must be of type symbolMapClass, not %s"%t
    if hasattr(target, 'mapNames'):
        ct = copy(target)
        ct.mapNames(themap)  # these methods work in place
        return ct
    elif isinstance(target, list):
        return themap(target)
    elif isinstance(target, tuple):
        return tuple(themap(target))
    elif isinstance(target, dict):
        o = {}
        for k, v in target.iteritems():
            o[themap(k)] = v
        return o
    elif isinstance(target, str):
        return themap(target)
    elif target is None:
        return None
    else:
        raise TypeError, "Invalid target type %s"%repr(type(target))


## internal functions for replaceSep
def replaceSepQSpec(spec, sourcesep, targetsep):
    # currently unused because many QuantSpec's with hierarchical names
    # are not properly tokenized! They may have the tokenized form
    # ['a_parent','.','a_child', <etc.>] rather than
    # ['a_parent.a_child', <etc.>]
    outstr = ""
    # search for pattern <nameToken><sourcechar><nameToken>
    try:
        for t in spec[:]:
            if isHierarchicalName(t, sourcesep):
                outstr += t.replace(sourcesep, targetsep)
            else:
                # just return original token
                outstr += t
    except AttributeError:
        raise ValueError("Invalid QuantSpec passed to replaceSep()")
    return outstr


def replaceSepStr(spec, sourcesep, targetsep):
    # spec is a string (e.g. 'leaf.v'), and p.tokenized contains the
    # spec split by the separator (e.g. ['leaf', '.', 'v']) if the separator
    # is not "_", otherwise it will be retained as part of the name
    outstr = ""
    treatMultiRefs = '[' in spec and ']' in spec
    p = parserObject(spec, treatMultiRefs=treatMultiRefs,
                     includeProtected=False)
    # search for pattern <nameToken><sourcechar><nameToken>
    # state 0: not in pattern
    # state 1: found nameToken
    # state 2: found nameToken then sourcechar
    # state 3: found complete pattern
    state = 0
    for t in p.tokenized:
        if isNameToken(t):
            if sourcesep in t:
                # in case sourcesep == '_'
                tsplit = t.split(sourcesep)
                if all([isNameToken(ts) for ts in tsplit]):
                    outstr += targetsep.join(tsplit)
                else:
                    outstr += t
            else:
                if state == 0:
                    state = 1
                    outstr += t
                elif state == 1:
                    state = 0
                    outstr += t
                else:
                    # state == 2
                    state = 1   # in case another separator follows
                    outstr += targetsep + t
            continue
        elif t == sourcesep:
            if state == 1:
                state = 2
                # delay output until checked next token
            else:
                # not part of pattern, so reset
                state = 0
                outstr += t
        else:
            state = 0
            outstr += t
    return outstr


def joinStrs(strlist):
    """Join a list of strings into a single string (in order)."""

    return ''.join(strlist)


def joinAsStrs(objlist,sep=""):
    """Join a list of objects in their string representational form."""

    retstr = ''
    for o in objlist:
        if type(o) is str:
            retstr += o + sep
        else:
            retstr += str(o) + sep
    avoidend = len(sep)
    if avoidend > 0:
        return retstr[:-avoidend]
    else:
        return retstr


def count_sep(specstr, sep=','):
    """Count number of specified separators (default = ',') in given string,
    avoiding occurrences of the separator inside nested braces"""
    
    num_seps = 0
    brace_depth = 0
    for s in specstr:
        if s == sep and brace_depth == 0:
            num_seps += 1
        elif s == '(':
            brace_depth += 1
        elif s == ')':
            brace_depth -= 1
    return num_seps


def parseMatrixStrToDictStr(specstr, specvars, m=0):
    """Convert string representation of m-by-n matrix into a single
    string, assuming a nested comma-delimited list representation in
    the input, and outputting a dictionary of the sub-lists, indexed
    by the ordered list of names specvars (specified as an
    argument)."""

    specdict = {}
    # matrix is n by m
    n = len(specvars)
    if n == 0:
        raise ValueError("parseMatrixStrToDictStr: specvars was empty")
    if m == 0:
        # assume square matrix
        m = len(specvars)
    # strip leading and trailing whitespace
    spectemp1 = specstr.strip()
    assert spectemp1[0] == '[' and spectemp1[-1] == ']', \
        ("Matrix must be supplied as a Python matrix, using [ and ] syntax")
    # strip first [ and last ] and then all whitespace and \n
    spectemp2 = spectemp1[1:-1].replace(' ','').replace('\n','')
    splitdone = False
    entrycount = 0
    startpos = 0
    try:
        while not splitdone:
            nextrbrace = findEndBrace(spectemp2[startpos:],
                                       '[', ']') + startpos
            if nextrbrace is None:
                raise ValueError, \
                      "Mismatched braces before end of string"
            specdict[specvars[entrycount]] = \
                                spectemp2[startpos:nextrbrace+1]
            entrycount += 1
            if entrycount < n:
                nextcomma = spectemp2.find(',', nextrbrace)
                if nextcomma > 0:
                    nextlbrace = spectemp2.find('[', nextcomma)
                    if nextlbrace > 0:
                        startpos = nextlbrace
                    else:
                        raise ValueError, \
                          "Not enough comma-delimited entries"
                else:
                    raise ValueError, \
                          "Not enough comma-delimited entries"
            else:
                splitdone = True
    except:
        print "Error in matrix specification"
        raise
    return specdict


def readArgs(argstr, lbchar='(', rbchar=')'):
    """Parse arguments out of string beginning and ending with round braces."""
    bracetest = argstr[0] == lbchar and argstr[-1] == rbchar
    return [bracetest, argstr[1:-1].replace(" ","").split(","), len(argstr)]


def findEndBrace(s, lbchar='(', rbchar=')'):
    """Find position in string (or list of strings) s at which final matching
    brace occurs (if at all)."""
    pos = 0
    assert s[0] == lbchar, 'string argument must begin with left brace'
    stemp = s
    leftbrace_count = 0
    notDone = True
    while len(stemp) > 0 and notDone:
        # for compatibility with s being a list, use index method
        try:
            left_pos = stemp.index(lbchar)
        except ValueError:
            left_pos = -1
        try:
            right_pos = stemp.index(rbchar)
        except ValueError:
            right_pos = -1
        if left_pos >= 0:
            if left_pos < right_pos:
                if left_pos >= 0:
                    leftbrace_count += 1
                    pos += left_pos+1
                    stemp = s[pos:]
                else:
                    # no left braces found. next brace is right.
                    if leftbrace_count > 0:
                        leftbrace_count -= 1
                        pos += right_pos+1
                        stemp = s[pos:]
            else:
                # right brace found first
                leftbrace_count -= 1
                pos += right_pos+1
                stemp = s[pos:]
        else:
            if right_pos >= 0:
                # right brace found first
                leftbrace_count -= 1
                pos += right_pos+1
                stemp = s[pos:]
            else:
                # neither were found (both == -1)
                raise ValueError, 'End of string found before closing brace'
        if leftbrace_count == 0:
            notDone = False
            # adjust for 
            pos -= 1
    if leftbrace_count == 0:
        return pos
    else:
        return None


# wrap objlist into a comma separated string of str(objects)
def makeParList(objlist, prefix=''):
    parlist = ', '.join(map(lambda i: prefix+str(i), objlist))
    return parlist


def wrapArgInCall(source, callfn, wrapL, wrapR=None, argnums=[0],
                  notFirst=False):
    """Add delimiters to single argument in function call."""
    done = False
    output = ""
    currpos = 0
    first_occurrence = True
    if wrapR is None:
        # if no specific wrapR is specified, just use wrapL
        # e.g. for symmetric delimiters such as quotes
        wrapR = wrapL
    assert isinstance(wrapL, str) and isinstance(wrapR, str), \
           "Supplied delimiters must be strings"
    while not done:
        # find callfn in source
        findposlist = [source[currpos:].find(callfn+'(')]
        try:
            findpos = min(filter(lambda x:x>=0,findposlist))+currpos
        except ValueError:
            done = True
        if not done:
            # find start and end braces
            startbrace = source[findpos:].find('(')+findpos
            endbrace = findEndBrace(source[startbrace:])+startbrace
            output += source[currpos:startbrace+1]
            # if more than one argument present, apply wrapping to specified
            # arguments
            currpos = startbrace+1
            numargs = source[startbrace+1:endbrace].count(',')+1
            if max(argnums) >= numargs:
                raise ValueError, "Specified argument number out of range"
            if numargs > 1:
                for argix in range(numargs-1):
                    nextcomma = source[currpos:endbrace].find(',')
                    argstr = source[currpos:currpos + nextcomma]
                    # get rid of leading or tailing whitespace
                    argstr = argstr.strip()
                    if argix in argnums:
                        if first_occurrence and notFirst:
                            output += source[currpos:currpos + nextcomma + 1]
                        else:
                            output += wrapL + argstr + wrapR + ','
                        first_occurrence = False
                    else:
                        output += source[currpos:currpos + nextcomma + 1]
                    currpos += nextcomma + 1
                if numargs-1 in argnums:
                    if first_occurrence and notFirst:
                        # just include last argument as it was if this is
                        # the first occurrence
                        output += source[currpos:endbrace+1]
                    else:
                        # last argument needs to wrapped too
                        argstr = source[currpos:endbrace]
                        # get rid of leading or tailing whitespace
                        argstr = argstr.strip()
                        output += wrapL + argstr + wrapR + ')'
                    first_occurrence = False
                else:
                    # just include last argument as it was
                    output += source[currpos:endbrace+1]
            else:
                if argnums[0] != 0:
                    raise ValueError, "Specified argument number out of range"
                if first_occurrence and notFirst:
                    output += source[currpos:endbrace+1]
                else:
                    argstr = source[currpos:endbrace]
                    # get rid of leading or tailing whitespace
                    argstr = argstr.strip()
                    output += wrapL + argstr + wrapR + ')'
                first_occurrence = False
            currpos = endbrace+1
        else:
            output += source[currpos:]
    return output


def addArgToCalls(source, callfns, arg, notFirst=''):
    """Add an argument to calls in source, to the functions listed in callfns."""
    # This function used to work on lists of callfns directly, but I can't
    # see why it stopped working. So I just added this recursing part at
    # the front to reduce the problem to a singleton function name each time.
    if isinstance(callfns, list):
        if len(callfns) > 1:
            res = source
            for f in callfns:
                res = addArgToCalls(res, [f], arg, notFirst)
            return res
    else:
        raise TypeError, "Invalid list of function names"
    done = False
    output = ""
    currpos = 0
    while not done:
        # find any callfns in source (now just always a singleton)
        findposlist_candidates = [source[currpos:].find(fname+'(') for fname in callfns]
        findposlist = []
        for candidate_pos in findposlist_candidates:
            # remove any that are actually only full funcname matches to longer strings
            # that happen to have that funcname at the end: e.g. functions ['if','f']
            # and find a match 'f(' in the string 'if('
            if currpos+candidate_pos-1 >= 0:
                if not isNameToken(source[currpos+candidate_pos-1]):
                    findposlist.append(candidate_pos)
            else:
                # no earlier character in source, so must be OK
                findposlist.append(candidate_pos)
        try:
            findpos = min(filter(lambda x:x>=0,findposlist))+currpos
        except ValueError:
            done = True
        if not done:
            # find start and end braces
            startbrace = source[findpos:].find('(')+findpos
            endbrace = findEndBrace(source[startbrace:])+startbrace
            sub_source = source[startbrace+1:endbrace]
            embedded_calls = [sub_source.find(fname+'(') for fname in callfns]
            try:
                subpositions = filter(lambda x:x>0,embedded_calls)
                if subpositions == []:
                    filtered_sub_source = sub_source
                else:
                    filtered_sub_source = addArgToCalls(sub_source,callfns,arg,notFirst)
            except ValueError:
                pass
            # add unchanged part to output
            # insert arg before end brace
            if currpos==0 and callfns == [notFirst]:
                notFirst = ''
                addStr = ''
            else:
                addStr = ', ' + arg
            output += source[currpos:startbrace+1] + filtered_sub_source \
                    + addStr + ')'
            currpos = endbrace+1
        else:
            output += source[currpos:]
    return output


if __name__ == '__main__':
    s1 = 'x+3^((t+1)^5)*(4-5)'
    s2 = 'x+3^((t+1)^5+x)+(4-5)'
    i1 = string2ast(s1)
    i2 = string2ast(s2)
    assert ast2string(toCircumflexSyntax(toPowSyntax(i1)))=='x+3.^(t+1)^5.*(4-5)'
    assert ast2string(toCircumflexSyntax(toPowSyntax(i2)))=='x+3.^((t+1)^5.+x)+(4-5)'
    assert ast2string(toCircumflexSyntax(toDoubleStarSyntax(toPowSyntax(i1))))=='x+3.^(t+1)^5.*(4-5)'
    assert ast2string(toPowSyntax(string2ast('pow(x,4)**4')))=='pow(pow(x,4),4)'
    assert ast2string(toCircumflexSyntax(string2ast('x**5**2')))=='x^5.^2.'
    j=string2ast('x+pow(t^5+1,2)**5*(4-5)')
    assert ast2string(toPowSyntax(j))=='x+pow(pow(pow(t,5)+1,2),5)*(4-5)'
    assert ast2string(toCircumflexSyntax(toDoubleStarSyntax(toPowSyntax(j))))=='x+(t^5.+1)^2.^5.*(4-5)'
    assert ast2string(simplify(toCircumflexSyntax(toDoubleStarSyntax(toPowSyntax(j)))))
    assert ast2string(toPowSyntax(string2ast('x+pow(t^5+1,2)^5*(4-5)')))=='x+pow(pow(pow(t,5)+1,2),5)*(4-5)'
    assert convertPowers('5*p^3*u+2*q^3')=='5*pow(p,3)*u+2*pow(q,3)'
    assert ast2string(ensureints(string2ast('x+(3.*pow(3.,y)^4.+1.)')))=='x+(3*pow(3,y)^4+1)'
    assert convertPowers('x+3**((t+1)**5+x)+(4-5)', '^')=='x+3^((t+1)^5+x)+(4-5)'
    assert ast2string(ensureparen_div(string2ast('(a+b)/x*y'))) == '((a+b)/x)*y'
    assert ast2string(ensureparen_div(string2ast('x/2*m'))) == '(x/2)*m'
    # -----------------------------------------------------------------
    s1='K.na(V)*(1-K.n)-K.nb(V)*K.n'
    s2='3-(cell1.K.na)*a.b'
    print replaceSep(s1)
    print replaceSep(s2)
    a21='a21'
    ab21='ab21'
    assert findNumTailPos(a21)==1
    assert findNumTailPos(ab21)==2
    assert isNameToken('a.b12')
    assert isNameToken('x300.y.z0')
    from Symbolic import QuantSpec, Var
    q1='v' - QuantSpec('k', '-2')
    assert q1.isCompound()
    q2 = QuantSpec('d', 'log(2)/mdt_216')
    assert q2.isCompound()
    q3 = QuantSpec('d', '1.026/(__mu_39)-32')
    assert q3.isCompound()
    print "----------------"
    q4=QuantSpec('q', 'exp((V + 27.0) / 5) - 1')
    assert q4.isCompound()
    q5=QuantSpec('q', 'exp(V+20)')
    assert not q5.isCompound()
    q6=QuantSpec('q', '(1-1/(a+bnb)*2)')
    assert not q6.isCompound()
    q7=QuantSpec('q', '(2*(t**-3))')
    assert str(q7.eval())=='2*Pow(t,-3)'
    q8=QuantSpec('q', '[[1,2],[3,var0]]')
    assert '[1' not in q8
    q9=Var('1/2*m','mf')
    assert str(q9.eval(m='0')) == '0'
    v = Var('10+4', 'v')
    assert v.spec[:]==['10','+','4']
    vfor = Var('1-I+for(atype,I,+)', 'fortest')
    vfor.mapNames({'I': 'dummy.I'})
    print vfor.spec[:]
    assert len(vfor.spec[:])==6 and vfor.spec[:][5]=='(atype,I,+)'
    print "------------------"
    v=Var('((k-(-10)))/(3+4)', 'v')
    print "v=Var('((k-(-10)))/(3+4)', 'v')"
    print "v.eval() -->", v.eval()
    print "(note simplification of braces around 10 and in the denominator)"
    print "v.eval(k=3) -->", v.eval(k=3)
    print "------------------"
    r=QuantSpec('dummy', 'kd1c1_140+kd2c1_144*(ec1n3_88*CLN3_25+ec1k2_86*BCK2_1+ec1n2_87*CLN2_24+ec1b5_85*CLB5_22+ec1b2_84*CLB2_20)/(Jd2c1_111+SIC1_50+C2_4+C5_6+SIC1P_51+C2P_5+C5P_7)')
    print r[:]
    print "------------------"
    sa='1 - call(a, b) + 4/6 - call(c, d)/2'
    print addArgToCalls(sa, ['call'], '_p')
    print addArgToCalls(sa, ['call'], '_p', notFirst='call')
    print "------------------"
    sb='def fun(x):\n\tcall(a, b) + 4/6 - call(c, d)/2'
    print addArgToCalls(sb, ['call','fun'], '_p')
    assert addArgToCalls(sb, ['call','fun'], '_p', 'fun') == addArgToCalls(sb, ['fun','call'], '_p', 'fun')
    assert addArgToCalls(sb, ['call','fun'], '_p', 'call') == addArgToCalls(sb, ['fun','call'], '_p', 'call')
    assert addArgToCalls(sb, ['another','call','fun','other'], '_p', notFirst='fun') == 'def fun(x):\n\tcall(a, b, _p) + 4/6 - call(c, d, _p)/2'
    sc='1/(call(a, b)) + 4/6'
    print addArgToCalls(sc, ['call'], '_p', notFirst='call')
    print addArgToCalls('zeta(yrel(y,0),z)-1',['zeta','yrel'],'_p')
    assert addArgToCalls('zeta(yrel(y,0),z)-1',['zeta','yrel'],'_p')=='zeta(yrel(y,0, _p),z, _p)-1'
    print addArgToCalls("2*__A4__*__A1__/(BB_218(__A1__,__A2__,__A3__,__A4__)+sqrt(pow(BB_218(__A1__,__A2__,__A3__,__A4__),2)-4*(__A2__-__A1__)*__A4__*__A1__))",['BB_218'],'p_')
    assert addArgToCalls("2*__A4__*__A1__/(BB_218(__A1__,__A2__,__A3__,__A4__)+sqrt(pow(BB_218(__A1__,__A2__,__A3__,__A4__),2)-4*(__A2__-__A1__)*__A4__*__A1__))",['BB_218'],'p_')=="2*__A4__*__A1__/(BB_218(__A1__,__A2__,__A3__,__A4__, p_)+sqrt(pow(BB_218(__A1__,__A2__,__A3__,__A4__, p_),2)-4*(__A2__-__A1__)*__A4__*__A1__))"
    fs3 = """double __rhs_ig(int cond_, double e1_, double e2_, double *p_, double *wk_, double *xv_) {
  if (cond_) {return e1_;} else {return e2_;};
}"""
    fs3_p = addArgToCalls(fs3, ['getindex', 'Na_ma', 'Na_mb', \
                                'globalindepvar', 'g', 'getbound', \
                                'Ileak', 'heav', 'Na_hb', 'Na_ha', \
                                '__rhs_ig', 'Ina', 'initcond'], "OH_NO", notFirst='__rhs_ig')
    if fs3 != fs3_p:
        print fs3_p
        raise AssertionError, "Mismatch!"
    print "------------------"
    sc='initcond(x, 1) + 3 - 2/(initcond(y3))'
    print wrapArgInCall(sc, 'initcond', '"', notFirst=False)
    print wrapArgInCall(sc, 'initcond', '"', notFirst=True)
