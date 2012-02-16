# Tests on Pointset labelling class: PointInfo
# Robert Clewley, November 2005

from PyDSTool import *
from PyDSTool.parseUtils import symbolMapClass

sm = symbolMapClass({'a': 'A'})

pstr = """p = PointInfo()
p[3] = ('a', {'bif': 'sn'})
p[5] = ('a', {'bif': 'h'})
p['a'] = (5, {'bif': 'h'})
p[1] = ('b', {'bif': 'c'})"""

print pstr
exec pstr

print "\nprint p -->\n", p, "\n"
print "p['a'] -->", p['a'], "\n"
print "p.sortByIndex() -->", p.sortByIndex(), "\n"

pupstr1 = "p.update(3, 'a', {'blah': 'tangle'})"
print pupstr1, "\n"
exec pupstr1

print "p[3] -->", p[3]
print "(Notice how the result is not a list because only a single"
print "label exists for index 3)\n"

pupstr2 = "p.update(3, 'c', {'foo': 'bar'})"
print pupstr2, "\n"
exec pupstr2

pupstr3 = "p.update(3, 'h')"
print pupstr3, "\n"
exec pupstr3

print "p.getIndices() -->", p.getIndices(), "\n"

try:
    p[-1]
except IndexError:
    print "Successfully could not access an index out of range"

try:
    p[-1] = "wrong"
except IndexError:
    print "Successfully could not add an index out of range"

try:
    p[50]
except IndexError:
    print "Successfully could not access an index that does not exist"


try:
    # index 10 not present!
    p.remove('a', 3, 10)
except KeyError:
    print "Successfully could not remove index label for non-existent entry"

del p['b']
try:
    del p['not_there']
except KeyError:
    print "Successfully could not delete index label for non-existent entry"

p.remove('a')
print "p.remove('a') --> [removed all indices associated with label 'a']\n"
assert len(p) == 1, "p ended up the wrong size!"

p[0] = ['a', ('k', {'bif': 'H', 'otherinfo': 'hello'})]
print "Mapping 'a' labels to 'A' using mapNames(<symbolMapClass instance>) method of PointInfo"
p.mapNames(sm)
print "\nprint p -->\n", p
