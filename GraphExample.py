from numpy import *
from numpy.linalg import svd,eig,det,inv

a=array([[3,2,2],[2,2,1],[2,1,2]])

u,s,v = svd(a)
#print u
print "singular values,",s
#print v
#print dot(dot(u,diag(s)),v)
#print dot(u,v)

print "det", det(a)

e,v= eig(a)
print "eigen values,",e
print "eigen vectos,",v 
iv=inv(v)
print "inv vectors,", iv
print "a:", dot(dot(v,diag(e)),iv)



print dot(dot(iv,a),v)

print dot(dot(dot(dot(v,v),v),v),v)


###########################################################################
###########################################################################
### Codes for graphs and vertex and ants from http://www.boriel.com/?lang=en
############################################################################
############################################################################


class Vertex(object):
    def __init__(self, node, *nodeList):
        self.i = node
        self.nodeList = list(nodeList)
 
    def __hash__(self):
        return self.i
 
    def reaches(self, vertex):
        ''' Can receive an integer or a Vertex
        '''
        if isinstance(vertex, int):
            return vertex in self.nodeList
 
        return self.reaches(vertex.i)
 
    def __str__(self):
        return '< ' + str(self.i) + '>'
 
    def __repr__(self):
        return self.__str__()
 
  
class Graph(object):
    def __init__(self):
        self.vList = {}
 
    def add(self, node, *nodeList):
        vertex = Vertex(node, *nodeList)
        self.vList[node] = vertex
 
    def hamiltonian(self, current = None, pending = None, destiny = None):
        ''' Returns a list of nodes which represent
        a hamiltonian path, or None if not found
        '''
        if pending is None:
            pending = self.vList.values()
 
        result = None
 
        if current is None:
            for current in pending:
                result = self.hamiltonian(current, [x for x in pending if x is not current], current)
                if result is not None:
                    break
        else:
            if pending == []:
                if current.reaches(destiny):
                    return [current]
                else:
                    return None
 
            for x in [self.vList[v] for v in current.nodeList]:
                if x in pending:
                    result = self.hamiltonian(x, [y for y in pending if y is not x], destiny)
                    if result is not None:
                        result = [current] + result
                        break   
 
        return result
 


################## First problem from El Pais
if 1== 1:
    G = Graph()
    # Defining the following graph
    #    /  8
    #   / /6 7\     
    #   1-2 3 4-5
    #      9 10
    #       11
    # Well not so easy to draw, see
    #http://www.elpais.com/videos/sociedad/Primer/desafio/matematico/PAIS/solucion/hay/solucion/elpepusoc/20110322elpepusoc_1/Ves/
    G.add(1, 2, 8, 11)
    G.add(2, 1, 6, 9)
    G.add(3, 6, 7, 9, 10)
    G.add(4, 5, 7, 10)
    G.add(5, 4, 8, 11)
    G.add(6, 2, 3, 8)
    G.add(7, 3, 4, 8)
    G.add(8, 1, 6, 7, 5)
    G.add(9, 2, 3, 11)
    G.add(10, 3, 4, 11)
    G.add(11, 1, 9, 10, 5)
    # Finding Hamiltonian path, but this graph has no Hamiltonian
    # bor being bipartito and having odd number of elements.
    print G.hamiltonian()
   


###########################################################################
##########################################################################
####   2nd problem from El Pais
##############################################



from random import choice

class Point(Vertex):
    '''A vertex which can have a poison trap
    '''
    def __init__(self, *args):
        Vertex.__init__(self, *args)
        self.poisoned = False
 
    def __str__(self):
        tmp = '*' if self.poisoned else ''
        return '<' + tmp + str(self.i) + '>'
 
  
class Cube(Graph):
    ''' A Graph derived class in which each vertex (Point)
    may contain a poison trap or not
    '''
    def add(self, node, *nodelist):
        vertex = Point(node, *nodelist)
        self.vList[node] = vertex
 
    def run(self, N):
        ''' Puts an ant on vertex #1 and let it run for N
        random steps. Stops when the ant enters a poisoned
        vertex or it has run N steps. Returns ant's last position.
        '''
        i = 0
        pos = 1 # Initial ant position
        while i < N:
            i += 1
            pos = choice(self.vList[pos].nodeList)
            if self.vList[pos].poisoned:
                break
 
        return pos
 
  
if  1 == 1:
    G = Cube()
    # Graph which is a cube
    #
    G.add(1, 2, 4, 5)
    G.add(2, 1, 3, 6)
    G.add(3, 2, 4, 7)
    G.add(4, 1, 3, 8)
    G.add(5, 1, 6, 8)
    G.add(6, 2, 5, 7)
    G.add(7, 3, 6, 8)
    G.add(8, 4, 5, 7)
 
    G.vList[7].poisoned = True
    G.vList[8].poisoned = True
    # Carry out simple statistics
    stats = [0, 0, 0] # Counters
 
    for i in xrange(1000000): # Total number of experiments (runs)
        pos = G.run(100) # make the ant to run 100 steps
        if pos > 6: # Map 0..6 => 0, 7 => 1, 8 => 2
            pos -= 6 # 7 = 1, 8 = 2
        else:
            pos = 0
 
        stats[pos] += 1 # Increase corresponding counter
 
    tot = sum(stats) # Total #number of experiments (must be the same number as above)
    print "Number of survivals, dead in 7, dead in 8",stats # Absolute Frequencies
    print "Relative probabilities"
    print [100 * float(x) / tot for x in stats] # Relative frequencies (percentage)
    print "From theory numbers are: 0, 3/7, 4/7=",0,3./7,4./7
