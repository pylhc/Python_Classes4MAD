import itertools
a=range(1,13)
codigos = itertools.combinations(a,6)
x = list(itertools.islice(codigos,0,1000))
print "Combinations of 12 elements 6 by 6 is", len(x)
tot=0
for c in x:
    for i in range(1,8):
        count=0
        for el in c:
            #print el,i, el>=i, el<=i+6, count
            if el>=i and el<= (i+5): 
                count=count+1
        if count==3:
            tot=tot+1
            print c, "range", i,i+5
            break
print "Total", tot 
