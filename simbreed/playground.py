# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 10:29:11 2016

@author: domi89
"""
import time

ID = [1, 2, 3, 4, 5, 6, 7]
father = ['NA', 'NA', 1, 1, 4, 5, 4]
mother = ['NA', 'NA', 2, 'NA', 3, 2, 5]

# unsorted
ID = [4, 1, 7, 5, 2, 6, 3]
father = [1, 'NA', 4, 4, 'NA', 5, 1]
mother = ['NA', 'NA', 5, 3, 'NA', 2, 2]

pop = pedigreePopulation.fromPedigree(ID=ID, father=father, mother=mother)
print(np.array_str(pop.makeA(method='LDL').matrix, precision=3))
print(np.array_str(pop.makeA(method='tabular').matrix, precision=3))


pop = pedigreePopulation()
pop.createFounders(n=4)
for _ in range(20):
    pop.randomMating(size=5)
print(pop)

chromLengths = [100.0] * 10
cS = CrossoverSimulator(m = 0.0, p = 0.0, obligateChiasma = False)

skPop = skeletonPopulation.fromPedigreePopulation(population=pop,
                                                  chromLengths=chromLengths,
                                                  crossoverSimulator=cS)
print(skPop)
Ind = skPop[skPop.allID()[-1]]
%timeit -n100 -r5 simulateMeiosis(parent=Ind, crossoverSimulator=cS)




pop.randomMating(size = 1000)

%timeit -n10 -r3 pop.makeA()



pop._scramble()
pop.isSorted()
%timeit pop.sort()

pop.recode()
print(pop)

pop.isFounder(4)


print(pop.subset(ID=(7,6)))
print(pop._subset(Ind=(pop[7], pop[6])))
Ind = (pop[7],)

pop.pairExpGeneContr(IDa=7, IDb=1)

print(pop.expGeneContrMatrix().matrix)
print(pop.isAncestorMatrix().matrix)


print(pop)
pop.subset(ID=('EGkARifW',))
pop['EGkARifW'].parents

# Population
numChrom = 10
chromLengths = [100] * numChrom
#numLociPerChrom = [1000] * numChrom
#positions = [sorted(np.random.uniform(0.0, h, n)) for h,n in zip(chromLengths, numLociPerChrom)]
#genome = Genome(chromLengths, positions)
cS = CrossoverSimulator(m = 3, p = 0.5, obligateChiasma = True)

skPop = skeletonPopulation.fromPedigreePopulation(population=pop,
                                                  chromLengths=chromLengths,
                                                  crossoverSimulator=cS)

print(np.array_str(skPop.makeA().matrix, precision=2))
print(np.array_str(2.0*skPop.realIBDMatrix().matrix, precision=2))
print(np.array_str(skPop.expGeneContrMatrix().matrix, precision=2))
#print(np.array_str(skPop.realFounderGeneContrMatrix(includeUnknown=True).matrix, precision=2))
print(np.array_str(skPop.realGeneContrMatrix().matrix, precision=2))


print(skPop)
skPop.realGeneContr(ID=7, maxGenBack=2)


skPop.pairRealGeneContr(IDa=4, IDb=6)
skPop.pairExpGeneContr(IDa=4, IDb=6)
print(np.array_str(skPop.expGeneContrMatrix().matrix, precision=2))
print(np.array_str(skPop.realGeneContrMatrix().matrix, precision=2))


  ID father mother generation     DH     F
0  1   None   None          0  False   0.0
1  2   None   None          0  False   0.0
2  3      1      2          1  False  None
3  4      1   None          1  False  None
4  5      4      3          2  False  None
5  6      5      2          3  False  None
6  7      4      5          3  False  None



time.time() - now                                                  
print(skPop)
print(np.array_str(2.0 * skPop.realIBDMatrix().matrix, precision=3))

skPop.realFounderGeneContr('HXVp6ksY')
print(np.array_str(skPop.realFounderGeneContrMatrix(includeUnknown=True).matrix, precision=2))




print(skPop)

# Pedigree
pedigree = Pedigree()
pedigree.createFounder()

rawped = rawPedigree()
F1 = Individual(rawped, 'F1', None, None)
F2 = Individual(rawped, 'F2', None, None)
I1 = Individual(rawped, 'I1', F1, F2)
I2 = Individual(rawped, 'I2', F1, F2)
I3 = Individual(rawped, 'I3', F1, I1)
I4 = Individual(rawped, 'I4', I2, None)
I5 = Individual(rawped, 'I5', I3, I4)
self.extend((I1, F1, I5, I2, F2, I4, I3))
print(self)
self.sort()
rawped.isSorted()
rawped.sort()


rawped.toDataFrame

# Falconer
#DH = [False, False, False, True, True, False, False, True, True, False]
ID = ['D', 'E', 'F', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'X']
father = ['A', 'J', 'C', 'E', 'G', 'H', 'I', 'K', 'L', 'N', 'P']
mother = ['B', 'B', 'B', 'F', 'F', 'F', 'J', 'J', 'M', 'O', 'Q']
#F = dict(XX = 0.5, G = -100, K = None, P = 12)
F = None
#for i,f,m,d in zip(ID,father,mother,DH): print(i,f,m,d)
pop = pedigreePopulation.fromPedigree(ID=ID, father=father, mother=mother)
print(pop)
pop.makeA().matrix
pop.inbreedingCoef().F

temp = [[1.000000, 0.500000, 0.25000, 0.3750000, 0.125000, 0.5000000, 0.1875000, 0.2500000, 0.3437500, 0.1250000, 0.2343750],
[0.500000, 1.000000, 0.25000, 0.6250000, 0.125000, 0.2500000, 0.3125000, 0.3750000, 0.2812500, 0.1875000, 0.2343750],
[0.250000, 0.250000, 1.00000, 0.6250000, 0.500000, 0.1250000, 0.3125000, 0.5625000, 0.2187500, 0.2812500, 0.2500000],
[0.375000, 0.625000, 0.62500, 1.1250000, 0.312500, 0.1875000, 0.5625000, 0.7187500, 0.3750000, 0.3593750, 0.3671875],
[0.125000, 0.125000, 0.50000, 0.3125000, 1.000000, 0.0625000, 0.1562500, 0.6562500, 0.1093750, 0.3281250, 0.2187500],
[0.500000, 0.250000, 0.12500, 0.1875000, 0.062500, 1.0000000, 0.0937500, 0.1250000, 0.5468750, 0.0625000, 0.3046875],
[0.187500, 0.312500, 0.31250, 0.5625000, 0.156250, 0.0937500, 1.0000000, 0.3593750, 0.5468750, 0.1796875, 0.3632812],
[0.250000, 0.375000, 0.56250, 0.7187500, 0.656250, 0.1250000, 0.3593750, 1.1562500, 0.2421875, 0.5781250, 0.4101562],
[0.343750, 0.281250, 0.21875, 0.3750000, 0.109375, 0.5468750, 0.5468750, 0.2421875, 1.0468750, 0.1210938, 0.5839844],
[0.125000, 0.187500, 0.28125, 0.3593750, 0.328125, 0.0625000, 0.1796875, 0.5781250, 0.1210938, 1.0000000, 0.5605469],
[0.234375, 0.234375, 0.25000, 0.3671875, 0.218750, 0.3046875, 0.3632812, 0.4101562, 0.5839844, 0.5605469, 1.0605469]]

A_true = np.array(temp)
A_tab = ped.makeA(ID = list(ped.keys()))
A_LDL = ped.makeA(ID = list(ped.keys()), method = 'LDL')
# subseting might suck for LDL, check it!!!
np.allclose(A_tab, A_LDL)

A_inv = ped.makeAinv(ID = list(ped.keys()), method = 'LDL')
np.allclose(A_inv, np.linalg.inv(A_tab))

ped.makeD(ID = ['D', 'E', 'F', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'X'])


now = time.time()
for i in range(1):
    getNewID(3000)
time.time() - now


pop = pedigreePopulation()
pop.createFounders(2)
pop.randomMating(size=3, returnID = True)
pop.inbreeding()
pop.recode()
print(pop)
for c in pop['E'].children: print(c.ID)

population = pop
print(pop.makeA().matrix)
pop.inbreedingCoef().F
pop.expGeneContrMatrix()
pop.isAncestorMatrix()
print(pop)

now = time.time()
for _ in range(10):
    pop.randomMating(size=100, returnID = True)
time.time() - now
#pop.inbreeding()

print(pop.pedigree.inbreedingCoef().F)

print(pop.pedigree.makeA(ID=['QPEsUk4V', 'zZ3UMqjh']).matrix)
pop.pedigree.allID()



self=pop.pedigree
print(self)
self.recode()
self.toDataFrame()
print(self)
self[3].father.ID

Ind = (next(reversed(pop.pedigree.values())),)

print(pop)
print(pop.pedigree.subset(Ind))



pop.pedigree.makeA(ID=pop.pedigree.allID())

#pop.inbreeding()
pop.doubledHaploidization

print(pop)
pop.pedigree['W7VrYJzr'].gametes[0].__dict__



f = populationFounderIndividual(None, None, False, 0.0, ID=None)
f.__dict__
f.population = pop
f.gametes[0].__dict__
f._update()


import time

now = time.time()
for _ in range(1000):
#   ped.makeA(ID = list(ped.keys()))
#   ped.makeA(ID = list(ped.keys()), method = 'LDL')
   ped.makeD(ID = ['D', 'E', 'F', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'X'])

time.time() - now




F1 = Entry('F1', None, None, 0, 0.0, False)
F2 = Entry('F2', None, None, 0, 0.0, False)
C1 = Entry('C1', F1.ID, F2.ID, 1, 0.0, False)
C2 = Entry('C2', F1.ID, F2.ID, 1, 0.0, False)
C3 = Entry('C3', C1.ID, C2.ID, 2, 0.0, False)
C4 = Entry('C4', C2.ID, C3.ID, 2, 0.0, False)

ped = Pedigree()
ped.extend((C4, C3, C2, C1, F1, F2))
ped.isSorted()
ped.sort()
ped.pairExpGeneContr('C4', 'C4')
ped.expGeneContrMatrix(generationA = 1, generationP = 3)['matrix']


ped.isAncestor('C4', 'C3')


print(ped)

ped.tabulateInd()
ped.isPadded()




# Population
numChrom = 10
chromLengths = [100] * numChrom
numLociPerChrom = [1000] * numChrom
positions = [sorted(np.random.uniform(0.0, h, n)) for h,n in zip(chromLengths, numLociPerChrom)]
genome = Genome(chromLengths, positions)
cS = CrossoverSimulator(m = 3, p = 0.5, obligateChiasma = True)
# from pedigree
ID=[1, 2, 3, 4, 5]
fatherID=[None, None, 1, 1, 3]
motherID=[None, None, 2, 3, 4]
pop = Population.fromPedigree(cS, genome, ID, fatherID, motherID)

print(pop.pedigree)
#pop.pedigree.makeA(pop.allID())['matrix']
#pop.pedigree.makeA(pop.allID(), method='LDL')['matrix']
#pop.pedigree.makeAinv(pop.allID())['matrix']
#pop.pedigree.makeAinv(pop.allID(), method='LDL')['matrix']
print(pop.pedigree.expGeneContrMatrix(pop.allID())['matrix'])
print(pop.realGeneContrMatrix(pop.allID())['matrix'])

pop.pairRealGeneContr(4,5)


pop.pedigree.cross(ID=6, fatherID=3, motherID=5)
pop.cross(ID=6, fatherID=3, motherID=5)
print(pop.pedigree)
pop.createFounders(n=3)




pop.pedigree.inbreeding(pop.allID())['F']

print(pop.pedigree.subset(ID = [1,2,5]))





pedigree = Pedigree.constructor(,
                           DH=None, F=None)
print(pedigree)
pop = Population.fromPedigree(pedigree, crossoverSimulator, genome, 'mypop')



print(pedigree)


pop = Population(crossoverSimulator, genome, 'pop')

nF = 2
pop.addFounders(nFounder = nF)
#print(pop.pedigree)

for _ in range(2):
    pop.randomMating(size = 2)


ID = pop.allID()
#IDb = pop.pedigree.getID(2)
pop.pedigree.makeA

print(pop.realIBDMatrix(ID)['matrix'])
print(pop.pedigree.expGeneContrMatrix(ID)['matrix'])
print(pop.realGeneContrMatrix(ID)['matrix'])

crossoverSimulator.simulateCrossover(60)

L = np.identity(n=10, dtype=np.float_)
L[np.ix_([i for i in range(3)],[i for i in range(3)])]



a = np.arange(10).reshape(2, 5)
a


ixgrid = np.ix_([0,1], [2,4])
ixgrid


ixgrid[0].shape, ixgrid[1].shape

a[ixgrid]

















import ordered_set
import timeit

o = ordered_set.OrderedSet()

setup = """
import ordered_set
import collections
n = 1000
every = 100
objlist = [object() for _ in range(n)]
findobjlist = [objlist[i] for i in range(0,n,every)]

oset = ordered_set.OrderedSet(objlist)
uoset = set(objlist)
olist = collections.OrderedDict.fromkeys(objlist, None)

"""


min(timeit.Timer(""" 
for i in findobjlist:
    i in olist
""", setup=setup).repeat(7, 1000))

d = {'a':5, 'v':6}



while(d):
    i = d.popitem()
    print(i)

for i in d.popitem():
    print(i)
    
ct = collections.Counter([1,2,2,3,4,5,5,5])
ct.update([1]*3) 
ct.update({1:3})   
    
   FounderGamete(self.genome.chromLengths,
              lineage=getNewID(),
              ID=getNewID())
    
FounderGamete(self.genome.chromLengths,
              lineage=getNewID(),
              ID=getNewID())

    
    
    
    
    
    
    
    
    
    
    
    