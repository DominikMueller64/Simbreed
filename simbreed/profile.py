# math
import numpy as np
import math
## tools
import itertools  # isSorted
import operator  # isSorted
import warnings
import collections  # namedtuple
## random numbers
import uuid
import shortuuid
import random
#
from numba import jit

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 


# Test fromPedigree
ID = [1, 2, 3, 4, 5, 6, 7]
father = ['NA', 'NA', 1, 1, 4, 5, 4]
mother = ['NA', 'NA', 2, 'NA', 3, 2, 5]

pop = pedPop.fromPedigree(ID=ID, father=father, mother=mother)

pop.randomMating(ID=[3,4], size=4, reproGroupGen=0)
pop
pop.subset(ID=1)


pop = pedPop()
pop.createFounders(n=10)
pop.roundRobin(size=1)
pop

pop.getReproGroups(ID=pop.getID(22), reproGroupGen=1)
pop.randomMating(size=3, reproGroupGen=1)

plotMatrix(pop.makeA(generation=pop.lastGen()).matrix)
#plt.imshow(pop.(generation=pop.lastGen()).matrix, cmap='hot', interpolation='nearest')



pop = skPop(chromLengths=chromLengths, crossoverSimulator=cS)

nFounder = 5
tmp = pop.createFounders(n=nFounder, returnID=True)
for _ in range(3): pop.randomMating(size = 10)

#tmp = pop.cross(fatherID=tmp[0], motherID=tmp[1], returnID=True)
#pop.inbreeding()
#pop.doubledHaploidization()
#pop.recombinationNumber()
pop.getLineage(pos=50.0, chrom=0)

pop.inbreedingCoef().F
print(np.round(pop.makeA().matrix, decimals=3))
print(np.round(2.0*pop.realIBDMatrix().matrix, decimals=4))

pop.mendelianSamplingVariance().D
pop._scramble()
#g=copy(pop)
pop.sort()
print(pop)



pop.randomMating(size=1000)
pop = genoPop.fromSkPop(population=pop, genome=genome,
                        founderFlesh=founderFlesh)
%timeit -n1 -r1 pop.addFlesh(ID=15560240492579)
geno = pop.genoMatrix()

genoty.asd

pop[85994593567319].gametes[0].skeleton

pop.recode()
print(pop)
pop.sort()
print(pop)
pop.realFounderGeneContrMatrix().matrix

#for _ in range(1):
#    pop.randomMating(size=1)
#pop.inbreeding()
#pop.doubledHaploidization()
#
#print(pop)
#print(.matrix)
#pop.makeA(method='LDL').matrix
#
#
#pop.inbreedingCoef().F
#
##pop.doubledHaploidization()
#
#pop.inbreedingCoef(ID=pop.allID()).F


numChrom = 2
chromLengths = [100.0] * numChrom
numLociPerChrom = [1000] * numChrom
cS = CrossoverSimulator(m = 10, p = 0.5, obligateChiasma = False)
positions = [sorted(np.random.uniform(0.0, h, n)) for h, n in zip(chromLengths, numLociPerChrom)]
locusNames = [getNewID(n) for n in numLociPerChrom]

genome = Genome(chromLengths, positions, locusNames)

pop = pedPop()
pop.createFounders(n=3)
pop.roundRobin(size=1)
pop.inbreeding(size=5)
plotMatrix(pop.makeA().matrix)
pop.getReproGroups(generation=2, reproGroupGen=1)
#pop.synthetic(size=5)

pop = skPop.fromPedPop(population=pop, chromLengths=chromLengths,
                       crossoverSimulator=cS)
plotMatrix(pop.realIBD().matrix)
pop.getSegmentBorders(ID=pop.getID(pop.lastGen()))



founderFlesh = {fid: [[np.random.choice(a=(0, 1), size=n) for n in numLociPerChrom]
                                                         for i in range(2)]
                                                         for fid in pop.getID(0)}
pop = genoPop.fromSkPop(population=pop, genome=genome,
                        founderFlesh=founderFlesh)


#pop.LD(ID=pop.getID(2), distMin=0, distMax=5, metric='Dcor2')
plotMatrix(pop.GRM().matrix)
plotMatrix(pop.GRM(method='SMC').matrix)


pop.pairwiseLD(ID=pop.getID(2), metric='Dcor2')
pop.pairwiseLD(ID=pop.getID(2), metric='Dprime')

pop.LD(ID=pop.getID(2), metric='Dcor2', distMin=0, distMax=1)


pop.LPS(generation=1, generationB=2, distMin=0, distMax=1,
        method='direction')

plotMatrix(pop.genoMatrix().matrix)
plotMatrix(pop.makeA().matrix)








pop.doubledHaploidization()






l = list()
for d in range(0, 100):
    l.append(pop.LD(ID=pop.getID(2), distMin=d, distMax=d+1.0))

pop.LD(ID=[40564130215411], distMin=d, distMax=d+1.0)



import matplotlib.pyplot as plt
plt.plot(l)




geno = pop.genoMatrix(ID=pop.allID())
geno = genoMatrix(matrix=geno.matrix, rownames=geno.rownames, colnames=geno.colnames)

geno.colnames
geno.freq()




for ind in pop.values():
    for j in range(2):
        print(ind.gametes[j].individual.ID)

self = ind.gametes[j]
ind.population
    
len(founderFlesh[62982782467952][0])
    
    
    for i in range(2):
        


print(pop)



for f in (pop[fid] for fid in pop.allID()[:nFounder]):
    for i in range(2):
        f.gametes[i].flesh = 





import timeit
setup = '''
from simbreed import pedigreePopulation, CrossoverSimulator
chromLengths = [100.0] * 10
cS = CrossoverSimulator(m = 10, p = 0.5, obligateChiasma = False)
pop = pedigreePopulation()
pop.createFounders(n=4)
'''

min(timeit.Timer('for _ in range(3): pop.randomMating(size=1000)',
                 setup=setup).repeat(5, 5))





%timeit -n100 -r10    simulateMeiosis(skPop[skPop.allID()[-1]], cS)

#if __name__ == "__main__":

#%timeit -n100 -r10 simulateMeiosis(Ind, cS)
