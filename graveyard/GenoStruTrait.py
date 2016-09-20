# -*- coding: utf-8 -*-
"""
Created on Sun Aug  7 08:28:51 2016

@author: dominik
"""

#from scipy.stats import poisson, uniform
#from blist import sortedlist
import scipy.optimize
import numpy as np
import bisect
import math


#    def __allUnique(self, x):
#        seen = set()
#        return not any(i in seen or seen.add(i) for i in x)
def dpoisson(k, l):
        return l**k * math.exp(-l) / math.factorial(k)


class Locus:
    def __init__(self, locusName, position = None, alleles = None):
        self.locusName = locusName
        self.alleles = alleles
        self.position = position
        self.chrom = None
        
        
    @property
    def locusName(self):
        return self.__locusName
    
    @locusName.setter
    def locusName(self, value):
        self.__locusName = value    
        
        
    @property
    def alleles(self):
        return self.__alleles
    
    @alleles.setter
    def alleles(self, value):
        self.__alleles = value   
        
    @property
    def position(self):
        return self.__position
    
    @position.setter
    def position(self, value):
        self.__position = value
        
        
    @property
    def chrom(self):
        return self.__chrom
    
    @chrom.setter
    def chrom(self, value):
        self.__chrom = value

    def __lt__(self, other):
        return self.position < other.position
    
    def __eq__(self, other):
        return self.__dict__ == other.__dict__


    def __repr__(self):
        return self.locusName
        
    #    def __repr__(self):
#        return "Locus\nlocusName = %s\nposition = %f\nallele = %s\nallelesAllowed = %s" %(self.locusName, self.position, self.allele, self.allelesAllowed)

# test class locus
l1 = Locus('l1', position = 12.34, alleles = None)
l1_1 = Locus('l1', position = 12.34, alleles = None)
l2 = Locus('l2', position = 12.95, alleles = None)
l1 == l1_1
l1 < l2
    

class Chromosome:
    
    def __init__(self, chromName = None, headPos = None, tailPos = None):
        self.loci = list()
        self.chromName = chromName
        self.gen = None
        self.numLoci = 0
        self.headPos = headPos
        self.tailPos = tailPos
        self.positions = list()
        self.lociNames = list()
        self.hashtable = dict()
        self.__headProvided = bool(headPos)
        self.__tailProvided = bool(tailPos)
    
        
    @property
    def chromName(self):
        return self.__chromName
    @chromName.setter
    def chromName(self, value):
        self.__chromName = value    
    
        
    @property
    def loci(self):
#        print("get loci")
        return self.__loci
        
    @loci.setter
    def loci(self, value):
#        print("set loci")
        self.__loci = value
               
    
    @property
    def numLoci(self):
        return self.__numLoci
    # make this private!   
    @numLoci.setter
    def numLoci(self, value):
        self.__numLoci = value
        
        
    @property
    def gen(self):
        return self.__gen
    @gen.setter
    def gen(self, value):
        self.__gen = value
        
        
    def addLoci(self, *args):
        for locus in args:
            if locus.locusName in self.lociNames:
                raise ValueError("All loci must be unique.")
            locus.chrom = self 
            self.numLoci += 1
            index = bisect.bisect_right(self.loci, locus)
            self.loci.insert(index, locus)
            self.positions.insert(index, locus.position)
            self.lociNames.insert(index, locus.locusName)
            self.hashtable[locus.locusName] = locus
            if not self.__headProvided:
                self.headPos = self.positions[0]
            if not self.__tailProvided:
                self.tailPos = self.positions[-1]
        
           
    
    def removeLoci(self, *args):
        for locus in args:
            if isinstance(locus, Locus):
                locusName = locus.locusName
            elif isinstance(locus, str):
                locusName = locus
            else:
                raise ValueError("Arguments must be either names or of type 'Locus'")
            # removal of loci will leave the order unchanged!
            self.getLocus(locusName).chrom = None
            self.numLoci -= 1
            index = self.locusIndex(locusName)
            del self.loci[index]
            del self.positions[index]
            del self.lociNames[index]
            del self.hashtable[locusName]  
    
    
    def lociIndices(self, *args):
        tmp = [self.lociNames.index(locusName) for locusName in args]
        return tmp if len(tmp) > 1 else tmp[0]
    
        
    def lociNames(self, *args):
        tmp = [self.lociNames[locusIndex] for locusName in args]
        return tmp if len(tmp) > 1 else tmp[0]
    
     
    def getLocus(self, locus):
        if isinstance(locus, int):
            return self.loci[locus]
        elif isinstance(locus, str):
#            return self.loci[self.lociIndices(locus)]
            return self.hashtable[locus]
        else:
            raise ValueError("'locus' must be either an index or a name.")
     
    def lociPos(self, *args):
        tmp = [self.getLocus(locus).position for locus in args]
        return tmp if len(tmp) > 1 else tmp[0]
                           
    def locusDist(self, locus1, locus2):
        return self.getLocus(locus2).position - self.getLocus(locus1).position
    
        

#    def sort(self):          
#        order = np.argsort(self.loci).tolist()
#        newLoci = self.loci.copy()
#        newPositions = self.positions.copy()
#        newLociNames = self.lociNames.copy()
#        for i,j in enumerate(order):
#            newLoci[i] = self.loci[j]
#            newPositions[i] = self.positions[j]
#            newLociNames[i] = self.lociNames[j]
#        self.loci = newLoci
#        self.positions = newPositions
#        self.lociNames = newLociNames
#        self.__sorted = True
#        
#    def isSorted(self):
#        return self.__sorted  
#        
#    def chromName(self):
#        return self.__chromName           
#    def chromBegin(self):
#        if self.__sorted:
#            return 0
#        else:
#            tmp = [x.position for x in self.loci]
#            return tmp.index(min(tmp))        
#    def chromEnd(self):
#        if self.__sorted:
#            return self.numLoci
#        else:
#            tmp = [x.position for x in self.loci]
#            return tmp.index(max(tmp))        
#    def lociAllNames(self):
#        return [x.locusName for x in self.loci]
#    def lociAllPos(self):
#        return [x.position for x in self.loci] 
#    def numLoci(self):
#        return self.numLoc
rs16 = Locus(locusName = 'rs16', position = 16, alleles = ['A', 'D'])
rs23 = Locus(locusName = 'rs23', position = 23, alleles = ['A', 'B'])
rs54 = Locus(locusName = 'rs54', position = 54, alleles = ['A', 'B'])
rs78 = Locus(locusName = 'rs78', position = 78, alleles = ['A', 'D'])
ts98 = Locus(locusName = 'ts98', position = 98, alleles = ['A', 'D'])
ts122 = Locus(locusName = 'ts122', position = 122, alleles = ['A', 'D'])        
 

#chrom1 = Chromosome(chromName = 'chrom1')
#chrom1.addLoci(rs78, rs23, rs16)
#chrom1.loci
#chrom1.positions
#chrom1.lociNames
#chrom1.numLoci
#chrom1.hashtable.get('rs16')
#chrom1.removeLoci('rs23', 'rs16')
#chrom1.loci
#chrom1.positions
#chrom1.lociNames
#chrom1.addLoci(rs23, rs16)
#chrom1.locusIndex('rs23')
#chrom1.locusName(1)
#chrom1.getLocus('rs23')
#chrom1.getLocus(1)
#chrom1.locusPos('rs23')



chrom1 = Chromosome(chromName = 'chrom1')
chrom1.addLoci(rs16, rs23, rs54, rs78)
chrom2 = Chromosome(chromName = 'chrom2')   
chrom2.addLoci(ts98, ts122) 


class Genome:
    def __init__(self, genName = None):
        self.chromosomes = list()
        self.genName = genName
        self.numChrom = 0
        self.chromNames = list()
        self.hashtable = dict()

    @property
    def genName(self):
        return self.__genName
    @genName.setter
    def genName(self, value):
        self.__genName = value     
        
    @property
    def chromosomes(self):
        return self.__chromosomes
        
    @chromosomes.setter
    def chromosomes(self, value):
        self.__chromosomes = value
            
        
    @property
    def numChrom(self):
        return self.__numChrom
    # make this private!   
    @numChrom.setter
    def numChrom(self, value):
        self.__numChrom = value
        
#    @property
#    def infoHash(self):
#        return self.__infoHash
#    @infoHash.setter
#    def infoHash(self, value):
#        self._infoHash = value
        
 
    def addChrom(self, *args):
        for chrom in args:
            keys = self.hashtable.keys()
            for locusName in chrom.lociNames:
                if locusName in keys:
                    raise ValueError('Locus %s is already present in the genome!' % locusName)
            chrom.gen = self
            self.numChrom += 1
            self.chromosomes.append(chrom)
            self.chromNames.append(chrom.chromName)
            self.hashtable.update(chrom.hashtable)
            
    def removeChrom(self, *args):
        for chrom in args:
            if isinstance(chrom, Chromosome):
                chromName = chrom.chromName
            elif isinstance(chrom, str):
                chromName = chrom
            else:
                raise ValueError("Arguments must be either names or of type 'Chromosome'")
            
            chrom = self.getChrom(chromName)
            chrom.gen = None
            self.numChrom -= 1
            index = self.chromIndex(chromName)
            del self.chromosomes[index]
            del self.chromNames[index]
            for locusName in chrom.lociNames:
                del self.hashtable[locusName]
            
        
    def chromIndex(self, chromName):
        return self.chromNames.index(chromName)
        
    def chromName(self, chromIndex):
        return self.chromNames[chromIndex]
              
    def getChrom(self, chrom):
        if isinstance(chrom, int):
            return self.chromosomes[chrom]
        elif isinstance(chrom, str):
            return self.chromosomes[self.chromIndex(chrom)]
                
    def numLoci(self, *args):
        if not args:
            args = genome.chromNames
        tmp = [self.getChrom(chrom).numLoci for chrom in args]
        return tmp if len(tmp) > 1 else tmp[0]
            
    def lociRelIndices(self, *args):
        tmp = [self.hashtable[locusName].chrom.lociIndices(locusName) for locusName in args]
        return tmp if len(tmp) > 1 else tmp[0]    
          
        
    #    def relLocusIndex(self, locusName, chrom = None):
#        if chrom:
#            return self.getChrom(chrom).locusIndex(locusName)
#        else:
#            for chrom in self.chromosomes:
#                try:
#                    index = chrom.locusIndex(locusNames)
#    
#    def absLocusIndex(self, chrom, locusName):
#        if isinstance(chrom, int):
#            index = chrom
#        elif isinstance(chrom, str):
#            index = self.chromIndex(chrom)
#        return sum([self.chromosomes[i].numLoci for i in range(index)]) + self.relLocusIndex(chrom, locusName)     
#    def numLoci(self, chromName):
#        return self.getChrom(chromName).numLoci
    
    
        
    
# test Genome
 
#genome = Genome(genomeName = 'genome')
#genome.infoHash
#genome.addChromosomes(chrom1, chrom2)
#genome.infoHash
#genome.numChrom
#genome.getChrom('chrom1')
#genome.getChrom(1)
#genome.chromosomes

genome = Genome(genomeName = 'genome')
genome.addChromosomes(chrom1, chrom2)

class Gamete:
    
    def __init__(self, gameteName, sequence, motherGamete = None, fatherGamete = None):
        # add boilerplate code for properties (motherGamete, fatherGamete)
#        self.genome = None
        self.ind = None
        self.motherGamete = motherGamete
        self.fatherGamete = fatherGamete
        self.gameteName = gameteName
        self.sequence = sequence
        
        
    
    @property
    def gameteName(self):
        return self.__gameteName
    @gameteName.setter
    def gameteName(self, value):
        self.__gameteName = value
            
    @property
    def sequence(self):
        return self.__sequence
    @sequence.setter
    def sequence(self, value):
        if not all([isinstance(seq, np.ndarray) for seq in value]):
            raise ValueError("All sequence data must be provided as numpy.ndarray.")
        self.__sequence = value
    
    @property
    def ind(self):
        return self.__ind
    @ind.setter
    def ind(self, value):
        self.__ind = value
        
sequence = [['A', 'D', 'A', 'C'], ['T', 'C']]
gam1 = Gamete('gam1', sequence)
sequence = [['D', 'D', 'T', 'T'], ['A', 'C']]
gam2 = Gamete('gam1', sequence)


class Individual:
    
    def __init__(self, indName, motherInd = None, fatherInd = None):
        # add boilerplate properties
        self.motherInd = motherInd
        self.fatherInd = fatherInd
        self.indName = indName,
        self.gametes = list()
        self.pop = None
        
    @property
    def indName(self):
        return self.__indName
    @indName.setter
    def indName(self, value):
        self.__indName = value
        
    
    @property
    def gametes(self):
        return self.__gametes
    @gametes.setter
    def gametes(self, value):
        self.__gametes = value
    
    @property
    def pop(self):
        return self.__pop
    @pop.setter
    def pop(self, value):
        self.__pop = value
        
    def addGamete(self, gam):
#        gam.genome = self.genome
        gam.ind = self
        self.gametes = self.gametes + [gam]
        


ind1 = Individual('ind1')
ind1.addGamete(gam1)
ind1.addGamete(gam2)
ind1.gametes


class Population:
    # checks must be implemented here
    # implement lot more further message for finden, removing etc etc individuals
    def __init__(self, popName, gen):
        self.popName = popName,
        self.gen = gen
        self.individuals = list()
        
    @property
    def popName(self):
        return self.__popName
    @popName.setter
    def popName(self, value):
        self.__popName = value
    
    @property
    def gen(self):
        return self.__gen
    @gen.setter
    def gen(self, value):
        self.__gen = value
     
    
    @property
    def individuals(self):
        return self.__individuals
    @individuals.setter
    def individuals(self, value):
        self.__individuals = value
    
        
    def addIndividual(self, ind):
        ind.pop = self
        self.individuals = self.individuals + [ind]  

pop = Population('pop', genome)
pop.addIndividual(ind1)
ind1.pop
ind = ind1
























class CrossingSimulator:
    
    def __init__(self, meiosisSimulator):
        self.meiosisSimulator = meiosisSimulator
    
    def cross(self, ind1, ind2, newIndName):       
        gamete1 = self.meiosisSimulator.meiosis(ind1, getNewId())
        gamete2 = self.meiosisSimulator.meiosis(ind2, getNewId())
        newInd = Individual(indName = newIndName,
                            motherInd = ind1,
                            fatherInd = ind2)
        newInd.addGamete(gamete1)
        newInd.addGamete(gamete2)
        
        return newInd
        

class RandomMatingSimulator:
    
    def __init__(self, crossingSimulator, allowSelfing = True):
        self.crossingSimulator = crossingSimulator
        self.allowSelfing = allowSelfing
        
    def randomMating(self, pop, size):
#        indAllNames = [ind.indName for ind in pop.individuals]
        
        for i in range(size):
            newIndName = getNewId()
            ind1 = pop.individuals[random.choice(range(len(pop.individuals)))]
            ind2 = pop.individuals[random.choice(range(len(pop.individuals)))]
            pop.addIndividual(self.crossingSimulator.cross(ind1, ind2, newIndName))
        



import uuid
import random
import time

def getNewId():
    return uuid.uuid4().hex.upper()


def createRandomGenome(numChrom, numLociPerChrom, sizesChrom):
    genome = Genome()
    for i in range(numChrom):
        numLoci = numLociPerChrom[i]
        sizeChrom = sizesChrom[i]        
        positions = np.random.uniform(low = 0.0, high = sizeChrom, size = numLoci).tolist()
        lociNames = [uuid.uuid4().hex.upper() for _ in range(numLoci)]
        loci = [Locus(locusName, position) for locusName, position in zip(lociNames, positions)]
        chrom = Chromosome(chromName = "chrom" + str(i), headPos = 0.0, tailPos = sizeChrom)
        chrom.addLoci(*loci)
        genome.addChrom(chrom)

    return genome


def createRandomPopulation(numInd, alleles, genome):
    
    pop = Population(popName = 'pop', gen = genome)
    numChrom = genome.numChrom
    for iInd in range(numInd):
        indName = uuid.uuid4().hex.upper()
        gameteName1 = uuid.uuid4().hex.upper()
        gameteName2 = uuid.uuid4().hex.upper()
        sequence1 = [[] for _ in range(numChrom)]
        sequence2 = [[] for _ in range(numChrom)]
        for iChrom in range(numChrom):
            numLoci = genome.numLoci()[iChrom]
            sequence1[iChrom] = np.array([random.choice(seq = alleles) for _ in range(numLoci)])
            sequence2[iChrom] = np.array([random.choice(seq = alleles) for _ in range(numLoci)])
        
        gamete1 = Gamete(gameteName1, sequence1)
        gamete2 = Gamete(gameteName2, sequence2)
        ind = Individual(indName)
        ind.addGamete(gamete1)
        ind.addGamete(gamete2)
        pop.addIndividual(ind)
    
    return pop
        

numChrom = 10
numLociPerChrom = [5000]*numChrom
sizesChrom = [100] * numChrom
numInd = 5
alleles = ['A', 'C']

genome = createRandomGenome(numChrom, numLociPerChrom, sizesChrom)
pop = createRandomPopulation(numInd, alleles, genome)

ind1 = pop.individuals[0]
ind2 = pop.individuals[1]

meiSim = MeiosisSimulator()
crossSim = CrossingSimulator(meiSim)
rMSim = RandomMatingSimulator(crossSim)


len(pop.individuals)

now = time.time()
rMSim.randomMating(pop, 1000)
time.time() - now

meiosisSimulator = MeiosisSimulator()
now = time.time()
for i in range(1000): meiosisSimulator.meiosis(ind1, 'gg')
time.time() - now

now = time.time()
sorted(uniform(0, 100).rvs(3).tolist())
time.time() - now

import numpy as np
now = time.time()
for _ in range(1000) :
#    poisson.rvs(5, size = 1)
    np.random.poisson(5, 1)
time.time() - now

newInd = crossingSimulator.cross(ind1, ind2, getNewId())
newInd.pop
    
    
gamete = meiosisSimulator().meiosis(pop.individuals[0], 'gg')
gamete.motherGamete
gamete.sequence
newGam
meiosisSimulator().meiosis(ind, 'gg').sequence





###########################################################

class Gamete:
    
    def __init__(self,
                 gameteName,
                 lineages = [],
                 locations = [],
                 motherGamete = None,
                 fatherGamete = None):

        self.ind = None
        self.lineages = lineages
        self.locations = locations
        self.motherGamete = motherGamete
        self.fatherGamete = fatherGamete
        self.gameteName = gameteName
        
        def addLineage(self, lineage):
            self.lineages.append(lineage)
            
        def addLocation(self, location):
            self.locatins.append(location)


class Individual:
    
    def __init__(self, indName, motherInd = None, fatherInd = None):
        self.motherInd = motherInd
        self.fatherInd = fatherInd
        self.indName = indName,
        self.gametes = list()
        self.pop = None
        
        
    def addGamete(self, gamete):
#        gam.genome = self.genome
        gamete.ind = self
        self.gametes.append(gamete)
        
class Population:
 
    def __init__(self, popName, genome):
        self.popName = popName,
        self.genome = genome
        self.individuals = list()
        
    def addIndividual(self, ind):
        ind.pop = self
        self.individuals.append(ind)
        
    def createParent(self, lineage1 = 1, lineage2 = 2):
        lineages1 = [[lineage1] for _ in range(genome.numChrom)]
        lineages2 = [[lineage2] for _ in range(genome.numChrom)]
        locations = [[chrom.tailPos] for chrom in genome.chromosomes]
        gamete1 = Gamete(getNewId(), lineages1, locations)
        gamete2 = Gamete(getNewId(), lineages2, locations)
        ind = Individual(indName = getNewId())
        ind.addGamete(gamete1)
        ind.addGamete(gamete2)
        
        return ind


   
            
            
class MeiosisSimulator:
    
    def __init__(self):
        pass
    
    
    def meiosis(self, ind, newGameteName):
        
        genome = ind.pop.gen
        gam1 = ind.gametes[0]
        gam2 = ind.gametes[1]
        newSequence = [[] for _ in range(genome.numChrom)]
        for i in range(genome.numChrom):
            chrom = genome.getChrom(i)
            pos = chrom.positions
            ln = chrom.tailPos - chrom.headPos
            nExpRec = ln * 1/100
            nRec = np.random.poisson(nExpRec)
#            chiaPos = sorted(np.random.uniform(chrom.headPos, chrom.tailPos, nRec).tolist())
            chiaPos = sorted([random.random() for _ in range(nRec)])
            
            s1 = gam1.sequence[i]
            s2 = gam2.sequence[i]

            if not chiaPos:
                if random.random() > 0.5:
                   newSequence[i] = s1
                else:
                   newSequence[i] = s2
            else:
                for chia in chiaPos:
                    idx = bisect.bisect_left(pos, chia)
                    if random.random() > 0.5:
                        newSequence[i] = np.concatenate((s1[:idx], s2[idx:]))
#                        newSequence[i] = s1[:idx] + s2[idx:]
                    else:
                        newSequence[i] = np.concatenate((s2[:idx], s1[idx:]))
#                        newSequence[i] = s2[:idx] + s1[idx:] 
                
        return Gamete(gameteName = newGameteName, 
                      sequence = newSequence,
                      motherGamete = gam1,
                      fatherGamete = gam2)     



#
#        curpos = lastpos+1;
#        loc = loc_clean;
#        alle = alle_clean;



parent = Population('pop1', genome).createParent()
child = simulateMeiosis(parent, cS)
child.









numChrom = 2
numLociPerChrom = [20]*numChrom
sizesChrom = [100] * numChrom

genome = createRandomGenome(numChrom, numLociPerChrom, sizesChrom)
cS = crossoverSimulator(m = 0, p = 0.5, obligateChiasma = True)






