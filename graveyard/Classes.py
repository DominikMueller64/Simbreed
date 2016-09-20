# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 15:23:04 2016

@author: domi89
"""

import scipy.optimize
import numpy as np
import bisect
import math
import uuid
import random
import time
from collections import OrderedDict
import heapq
import operator
import itertools
import sys

def getNewId():
    return uuid.uuid4().int


def dpoisson(k, l):
        return l**k * math.exp(-l) / math.factorial(k)
        
    
def is_sorted(iterable, reverse = False):
    """Check if the iterable is sorted, possibly reverse sorted."""

    def pairwise(iterable):
        """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)

    my_operator = operator.ge if reverse else operator.le

    return all(my_operator(current_element, next_element)
                   for current_element, next_element in pairwise(iterable))


def order(iterable, reverse = False):
    return sorted(range(len(iterable)), key=lambda k: iterable[k], reverse = reverse)
    
def reorder(iterable, order):
    """Reorder a list according to the indices specified in 'order'"""
    new = iterable.copy()
    for i, o in enumerate(order):
        new[i] = iterable[o]
    return new  



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

### Genome must be simplified, because there is a lot of boilerplate stuff that
### will not be used in 99% of the cases, like adding individual chromosomes or
### even loci. Also, creating a genome by doing it incrementally is nice for 
### practicing programming, but practically useless, because nobody will ever do it
### I need to keep the data structures as simple and plain as possible and avoid
### to much boilerplate code shit.
### Simply: For the abstract simulation, a genome does not need to be more
### than a list/vector of chromosomal lengths, and this is everything. So why
### make it more than than?


        
        
class crossoverSimulator:
    
    def __init__(self, m, p, obligateChiasma = False):
        self.m = m
        self.p = p
        self.obligateChiasma = obligateChiasma
        
    
    
    def calculateLStar(self, L):
        m = self.m
        p = self.p
        
        if L <= 50:
            raise ValueError("Chromosome must have length > 50 cM")
        if m < 0 or not isinstance(m, int):
            raise ValueError("Parameter 'm' must be a non-negative integer.")
        
        if p < 0 or p > 1:
            raise ValueError("Parameter 'p' must be in [0, 1]")
        if p == 1:
            # if p == 1, might as well take m = 0, p = 0
            m = 0
            p = 0
            
    
            
        def funcToZero(Lstar, L, m, p):           
            if m == 0:
                denom = 1 - math.exp(-Lstar / 50)
            else:
                lambda1 = Lstar/50 * (m + 1) * (1 - p)
                lambda2 = Lstar/50 * p
                
                sm = 0.0
                for i in range(m + 1):
                    sm += dpoisson(i, lambda1) * (m + 1 - i)/ (m + 1)
                
                denom = 1 - sm * math.exp(-lambda2)
            
            return 2 * L - 2 * Lstar / denom
     
     
     
    
        return scipy.optimize.brentq(f = funcToZero, a = math.exp(-8), b = L, args = (L, m, p))  
                           
    
    def simulateCrossover(self, L):
        m = self.m
        p = self.p                       
        obligateChiasma = self.obligateChiasma
#        print(m,p,obligateChiasma)
        
        if obligateChiasma:
            Lstar = self.calculateLStar(L)
        else:
            Lstar = L
            
        # no-interference model is a lot easier
        if m == 0:
            
            if obligateChiasma:
                # rejection sampling to get at least one chiasma
                while True:
                    nXO = np.random.poisson(Lstar/50.0)
                    if nXO > 0:
                        break
                nXO = np.random.binomial(nXO, 0.5)
            else:
                nXO = np.random.poisson(L/100.0)
            
            return sorted(np.random.uniform(0.0, L, nXO).tolist())
            
        
        lambda1 = Lstar/50.0 * (m + 1) * (1.0 - p)
        lambda2 = Lstar/50.0 * p
        
        while True:
            # chiasma and intermediate points
            nPoints = np.random.poisson(lambda1)
            
            # which point is the first chiasma?
            first = np.random.choice(range(m + 1))
            if first > nPoints:
                nIchi = 0
            else:
                nIchi = nPoints / (m + 1) + int(first < (nPoints % (m + 1)))

            # no. chiasma from no interference process
            if p > 0:
                nNIchi = np.random.poisson(lambda2)
            else:
                nNIchi = 0
                
            if not obligateChiasma or nIchi + nNIchi > 0:
                break
        
        # locations of chiasmata and intermediate points for process w/ interference
        pointLocations = sorted(np.random.uniform(0.0, L, nPoints).tolist())
        
        # move every (m+1)st point back to front
        nChi = 0
        for j in range(first, nPoints - 1, m + 1):
            # Here, j > nChi is required, otherwise a mess will happen
            pointLocations[nChi] = pointLocations[j]
            nChi += 1
        
        # chiasma locations from non-interference process
        NIchiLocations = np.random.uniform(0.0, L, nNIchi).tolist()
        
        # combine interference and no interference chiasma locations
        chiLocations = sorted(pointLocations[:nChi] + NIchiLocations)
        
        # thin by 1/2
        nXO = 0
        XOLocations = list()
        for i in range(len(chiLocations)):
            if random.random() < 0.5: # flip coin -> chiasma
                nXO += 1
                XOLocations.append(chiLocations[i])
                
        return XOLocations



def simulateMeiosis(parent, cS):
    
    pat = parent.gametes[0].skeleton
    mat = parent.gametes[1].skeleton 
    numChrom = len(mat)
    skeleton = [{} for _ in range(numChrom)]

    for iChrom in range(numChrom):
        
        matalle = mat[iChrom]['lineages']
        matloc = mat[iChrom]['locations']
    
        patalle = pat[iChrom]['lineages']
        patloc = pat[iChrom]['locations']
        
        L = matloc[-1]
        
        # simulate crossover locations; add -1 to the beginning
        tmp = cS.simulateCrossover(L)
        product = [-1] + tmp
    

        cur_allele = random.choice((0,1)) # first allele (0 or 1)
    
        loc = list()
        alle = list()
        curpos = 0
        
        if len(product) == 1:
            if cur_allele == 0:
                skeleton[iChrom]['lineages'] = matalle
                skeleton[iChrom]['locations'] = matloc
            else:
                skeleton[iChrom]['lineages'] = patalle
                skeleton[iChrom]['locations'] = patloc
                
        else:
            for i in range(1, len(product)):
                leftLoc = product[i - 1]
                rightLoc = product[i]
                
#                if cur_allele == 0:
#                    for j in range(len(matloc)):
#                        ml = matloc[j]
#                        ma = matalle[j]
#                        if leftLoc <= ml < rightLoc:
#                            loc.append(ml)
#                            alle.append(ma)
#                        elif rightLoc < ml:
#                            break
#                    loc.append(rightLoc)
#                    alle.append(ma)
#                        
#                else:
#                    for j in range(len(patloc)):
#                        pl = patloc[j]
#                        pa = patalle[j]
#                        if leftLoc <= pl < rightLoc:
#                            loc.append(pl)
#                            alle.append(pa)
#                        elif rightLoc < pl:
#                            break  
#                    loc.append(rightLoc)
#                    alle.append(pa)
                    
                if cur_allele == 0:
                    parloc = matloc
                    paralle = matalle
                else:
                    parloc = patloc
                    paralle = patalle
                    
                for j in range(len(parloc)):
                    pl = parloc[j]
                    pa = paralle[j]
                    if leftLoc <= pl < rightLoc:
                        loc.append(pl)
                        alle.append(pa)
                    elif rightLoc < pl:
                        break  
                loc.append(rightLoc)
                alle.append(pa)
    
                cur_allele = 1 - cur_allele
            
            lastxo = product[-1];
                 
            
#            if cur_allele == 0:
#                for j in range(len(matloc)):
#                    if matloc[j] > lastxo:
#                        loc.append(matloc[j])
#                        alle.append(matalle[j])
#            else:
#                for j in range(len(patloc)):
#                    if patloc[j] > lastxo:
#                        loc.append(patloc[j])
#                        alle.append(patalle[j])
                        
            if cur_allele == 0:
                parloc = matloc
                paralle = matalle
            else:
                parloc = patloc
                paralle = patalle            
            
            for j in range(len(parloc)):
                if parloc[j] > lastxo:
                    loc.append(parloc[j])
                    alle.append(paralle[j])
            
            
                        
            curpos = len(loc)            
            if curpos > 1: #clean up repeated alleles
                loc_clean = [None] * curpos
                alle_clean = [None] * curpos
                loc_clean[0] = loc[0]
                alle_clean[0] = alle[0]
                lastpos = 0
                
                for i in range(1, len(loc)):
                    if alle_clean[lastpos] == alle[i]:
                        loc_clean[lastpos] = loc[i]
                    else:
                        lastpos += 1
                        loc_clean[lastpos] = loc[i]
                        alle_clean[lastpos] = alle[i]
                
                curpos = lastpos + 1
                loc = loc_clean[0:curpos]
                alle = alle_clean[0:curpos]
                
            skeleton[iChrom]['lineages'] = alle
            skeleton[iChrom]['locations'] = loc
                        
    return Gamete(name = getNewId(), father = parent.gametes[0],
                  mother = parent.gametes[0], skeleton = skeleton)
                  
            


def cross(father, mother, cS, name = None):
    gamete1 = simulateMeiosis(father, cS)
    gamete2 = simulateMeiosis(mother, cS)
    if not name:
        name = getNewId() 
    return Individual(name = name, father = father,
                      mother = mother, gametes = [gamete1, gamete2])
    

class Skeleton:
    
    def __init__(self, ):
        self.value = value


#def createFounders(nFounder):
                                
class Gamete:
    
    def __init__(self, name = None, father = None, mother = None,
                 skeleton = None, flesh = None):
        self.name = name
        self.father = father
        self.mother = mother
        self.skeleton = skeleton
        self.flesh = flesh
              
class FounderGamete(Gamete):
    
    def __init__(self, chromLengths, lineage, name = None, flesh = None):
        skeleton = [ {'lineages' : [lineage], 'locations' : [chromLen]} for chromLen in chromLengths]
        Gamete.__init__(self, name = name, skeleton = skeleton, flesh = flesh)
        self.lineage = lineage 

class Individual:
    
    def __init__(self, name = None, father = None, mother = None,
                 gametes = None):
        self.name = name
        self.father = father
        self.mother = mother
        self.gametes = gametes
#        self.pop = None


class Subpopulation:
    
    def __init__(self, name = None, genome = None, pedigree = None):
        self.name = name
        self.founderGametes = dict()
        self.individuals = dict()
        self.genome = genome
        self.pedigree = pedigree
        
    def __addFounderGamete(self, gamete):
        if isinstance(gamete, FounderGamete):
            self.founderGametes[gamete.lineage] = gamete
        else:
            raise ValueError('Not a founder gamete.')
        
    def add(self, obj):
        if isinstance(obj, FounderGamete):
            self.__addFounderGamete(obj)
        elif isinstance(obj, Individual):
            gamete1 = obj.gametes[0]
            gamete2 = obj.gametes[1]
            if isinstance(gamete1, FounderGamete):
                self.__addFounderGamete(gamete1)
            if isinstance(gamete2, FounderGamete):
                self.__addFounderGamete(gamete2)
            self.individuals[obj.name] = obj
    
    def getInd(self, name):
        return self.individuals[name]
        
    def realIBD(self, generation):
        """Computes the realized IBD for all pairs"""
        Id = self.pedigree.getId(generation = generation)
        nId = len(Id)
        IBD = np.empty(shape = (nId, nId), dtype = np.float_)
        
        for i in range(nId):
            for j in range(i):
                IBD[i,j] = self.pairwiseRealIBD(Id[i], Id[j])
                IBD[j,i] = IBD[i,j]
            IBD[i,i] = self.pairwiseRealIBD(Id[i])
            
        return IBD


        
    def pairwiseRealIBD(self, id1, id2 = None):
        ind1 = subpop.getInd(id1)
        if id2:
            ind2 = subpop.getInd(id2)
        else:
            ind2 = subpop.getInd(id1)
            
        IBD = 0.25 * (self.realIBDGam(ind1.gametes[0], ind2.gametes[0])
                      + self.realIBDGam(ind1.gametes[0], ind2.gametes[1])
                      + self.realIBDGam(ind1.gametes[1], ind2.gametes[0])
                      + self.realIBDGam(ind1.gametes[1], ind2.gametes[1]))
        return IBD
        
                
    def realIBDGam(self, gam1, gam2):
        gam1 = gam1.skeleton
        gam2 = gam2.skeleton
        numChrom = len(gam1)
        chromLengths = [chrom['locations'][-1] for chrom in gam1]
        cumlen = 0.0
        for i, chrom1, chrom2 in zip(range(numChrom), gam1, gam2):

            loc1, loc2 = chrom1['locations'], chrom2['locations']
            lin1, lin2 = chrom1['lineages'], chrom2['lineages']
#            loc1 = [10, 15, 30, 55, 100]
#            loc2 = [20, 40, 50, 70, 100]
#            lin1 = [1,2,1,1,2]
#            lin2 = [1,1,2,1,2]
            
            # alternative code, simpler algorithm, should be correct
            ixi1, ixi2 = 0, 0
            left = 0.0
            l1, l2 = lin1[ixi1], lin2[ixi2]
            while ixi1 < len(loc1) and ixi2 < len(loc2):
                # find the interval
                if loc1[ixi1] <= loc2[ixi2]:
                    right = loc1[ixi1]
                    l1 = lin1[ixi1]
                    l2 = lin2[ixi2]
                    ixi1 += 1
                else:
                    right = loc2[ixi2]
                    l1 = lin1[ixi1]
                    l2 = lin2[ixi2]
                    ixi2 += 1
                    
                if l1 == l2:
                    cumlen += right - left
                 
                left = right
                  
#            # This should work as well, but terribly obfuscated
#            # initialize
#            cumlen = 0.0
#            ix1, ix2 = 0, 0
#            left = 0.0
#            if loc1[ix1] <= loc2[ix2]:
#                which = 0
#                right = loc1[ix1]
#            else:
#                which = 1
#                right = loc2[ix2]
# 
#                
#            flag = True    
#            while flag:
#                if which == 0:
#                    while loc1[ix1] <= loc2[ix2]:
#                        if lin1[ix1] == lin2[ix2]:
#                            cumlen += right - left
#                        ix1 += 1
#                        left = right
#                        if ix1 == len(loc1):
#                            flag = False
#                            break
#                        
#                        if loc1[ix1] <= loc2[ix2]:
#                            right = loc1[ix1]
#                        else:
#                            right = loc2[ix2]
#                                           
#                else:
#                    while loc2[ix2] <= loc1[ix1]:
#                        if lin1[ix1] == lin2[ix2]:
#                            cumlen += right - left
#                        ix2 += 1
#                        left = right
#                        if ix2 == len(loc2):
#                            flag = False
#                            break
#                        
#                        if loc2[ix2] <= loc1[ix1]:
#                            right = loc2[ix2]
#                        else:
#                            right = loc1[ix1]
#                        
#                which = 1 - which
        
#            sum(x  for x, y in zip(inCommon, chromLengths)) / sum(chromLengths)
        #           list(OrderedDict.fromkeys(heapq.merge(chr1['locations'], chr2['locations'])))

        return cumlen / sum(chromLengths)
    

    def addFleshGamete(self, gamete):
        """Add Genotypes to a single gamete"""
        numChrom = len(gamete.skeleton)
        flesh = [None] * numChrom
        for ichrom in range(numChrom):
            positions = self.genome.chromosomes[ichrom].positions
            lin = gamete.skeleton[ichrom]['lineages']
            loc = gamete.skeleton[ichrom]['locations']
            tmp = [None] * len(positions)
            ixl = 0
            for lc, li in zip(loc, lin):
                ixr = bisect.bisect(positions, lc)
                # find founder
                tmp[ixl:ixr] = self.founderGametes[li].flesh[ichrom][ixl:ixr]
                ixl = ixr
            flesh[ichrom] = tmp
        gamete.flesh = flesh
    
    
    def addFlesh(self, generation = None):
        ## might also add flesh to founder gametes, but should
        ## change nothing. Cater for later on.
        if not generation:
            allGenerations = list(set(self.pedigree.generation))
        else:
            allGenerations = [generation]
            
        for generation in allGenerations:
            names = self.pedigree.getId(generation = generation)
            for nm in names:
                ind = subpop.getInd(nm)
                self.addFleshGamete(ind.gametes[0])
                self.addFleshGamete(ind.gametes[1])
                
    def genoMatrixDict(self, generation = None):
        names = self.pedigree.getId(generation)   
        tmp = {}
        for nm in names:
            ind = self.getInd(nm)                
            tmp[nm] = np.vstack((np.hstack(ind.gametes[0].flesh),
                                 np.hstack(ind.gametes[1].flesh)))
        return tmp   
                
        
        
        


class Pedigree:
    
    def __init__(self, nFounders):
        """Instantiate a 'Pedigree' object with a specified number of founder individuals."""
        self.nFounders = nFounders
        self.id = [getNewId() for _ in range(nFounders)]
        self.mother = [None] * nFounders
        self.father = [None] * nFounders
        self.generation = [0] * nFounders

#        self.curid = self.id.copy()
#        self.curgen = 0
        
    # Accessor methods
    def numAllInd(self):
        """Get the number of all individuals in the pedigree"""
        # retrun len(set(self.id)) if no trust in uniqueness
        return len(self.id)
    
#    def numCurInd(self):
#        """Get the number of the current individuals in the pedigree"""
#        return len(self.curid)
        
    def lastGen(self):
        """Get the number of the last generation in the pedigree"""
        # more efficient if assumed sorted: self.generation[-1]
        return max(set(self.generation))
        
    def numGen(self):
        """Get the number of generations in the pedigree"""
        return len(set(self.generation))
        
        
    def getId(self, generation = None):
        """ Get the IDs of all individuals in 'generation'.
        If generation == None, all IDs are returned.
        Default is the last generation"""
        if not generation:
            return self.id
        else:
            ret = []
            for g, i in zip(self.generation, self.id):
                if g == generation:
                    ret.append(i)
        return ret
        
    def isSorted(self):
        """Check if the pedigree is sorted according to generations."""
        return is_sorted(self.generation)
    
    def sort(self):
        """Sort the pedigree according to generations."""
        index = order(self.generations)
        self.generations = reorder(self.generations, index)
        self.id = reorder(self.id, index)
        self.father = reorder(self.father, index)
        self.mother = reorder(self.mother, index)
         
    # These mating methods modify the Pedigree object
        
    def randomMating(self, size, allowSelfing = True):
        """Perform random mating of the current individuals."""
        curgen = self.lastGen()
        curid = self.getId(curgen).copy()
        if not allowSelfing and len(curid) < 2:
            raise ValueError("It is not possible to have 'allowSelfing = False' if there are less than 2 individuals in the population.")
        
        newId = [getNewId() for _ in range(size)]
        self.id.extend(newId)
        
        for i in range(size):
            id1 = random.choice(curid)
            id2 = random.choice(curid)
            if not allowSelfing:
                while id1 == id2:
                    id2 = random.choice(curid)
            self.mother.append(id1)
            self.father.append(id2)
        
        self.generation.extend([curgen + 1] * len(newId))
        
        
    def selfing(self):
        """Perform self-fertilization of all current individuals."""
        curgen = self.lastGen()
        curid = self.getId(curgen).copy()
        newId = [getNewId() for _ in range(len(curid))]
        self.id.extend(newId)
        self.mother.extend(curid)
        self.father.extend(curid)
        self.generation.extend([curgen + 1] * len(newId))
        
    def roundRobin(self):
        """Perform mating according to a round robin design"""
        curgen = self.lastGen()
        curid = self.getId(curgen).copy()
        newId = [getNewId() for _ in range(len(curid))]
        self.id.extend(newId)            
        for id1, id2 in zip(curid, curid[1:] + curid[:1]):
            self.mother.append(id1)
            self.father.append(id2)
        self.generation.extend([curgen + 1] * len(newId))
        

     
        
def simulatePedigree(pedigree, chromLengths, cS, name = None):
    if not isinstance(pedigree, Pedigree):
        raise ValueError("Argument 'pedigree' must be of class 'Pedigree'.")
    if not pedigree.isSorted():
        raise ValueError("The 'pedigree' must be sorted.")
        
    pop = Subpopulation(name = name)
    pop.pedigree = pedigree

    ctFGam = 0
    for i in range(len(pedigree.id)):
        ID = pedigree.id[i]
        father = pedigree.father[i]
        mother = pedigree.mother[i]
        gen = pedigree.generation[i]
             
        if father:
            fatherGamete = simulateMeiosis(pop.getInd(father), cS)
        else:
            fatherGamete = FounderGamete(chromLengths, lineage = ctFGam, name = getNewId())
            ctFGam += 1
            
        if mother:
            motherGamete = simulateMeiosis(pop.getInd(mother), cS)
        else:
            motherGamete = FounderGamete(chromLengths, lineage = ctFGam, name = getNewId())
            ctFGam += 1
        
        newInd = Individual(name = ID, father = father, mother = mother,
                            gametes = [fatherGamete, motherGamete])
        
        pop.add(newInd)
    
    return pop                


def createRandomGenome(numChrom, numLociPerChrom, sizesChrom):
    genome = Genome()
    for i in range(numChrom):
        numLoci = numLociPerChrom[i]
        sizeChrom = sizesChrom[i]        
        positions = np.random.uniform(low = 0.0, high = sizeChrom, size = numLoci).tolist()
        lociNames = [getNewId() for _ in range(numLoci)]
        loci = [Locus(locusName, position) for locusName, position in zip(lociNames, positions)]
        chrom = Chromosome(chromName = "chrom" + str(i), headPos = 0.0, tailPos = sizeChrom)
        chrom.addLoci(*loci)
        genome.addChrom(chrom)

    return genome
    
def createRandomFlesh(numLociPerChrom):
    return [np.random.choice(a = [0,1], size = numLoci, replace = True).tolist() for numLoci in numLociPerChrom]
     




numChrom = 10
chromLengths = [100] * numChrom
cS = crossoverSimulator(m = 4, p = 0.5, obligateChiasma = False) 
pedigree = Pedigree(2)

for i in range(10):
    pedigree.randomMating(size = 1000)
#    pedigree.selfing()
    
subpop = simulatePedigree(pedigree, chromLengths, cS)
subpop.getInd(pedigree.getId(pedigree.lastGen())[0]).gametes[0].skeleton
subpop.getInd(pedigree.getId(pedigree.lastGen())[0]).gametes[1].skeleton

subpop.pairwiseRealIBD(pedigree.getId(pedigree.lastGen())[0])

now = time.time()
print(subpop.realIBD(generation = pedigree.lastGen()))
time.time() - now

numLociPerChrom = [1000] * numChrom
genome = createRandomGenome(numChrom = numChrom, 
                            numLociPerChrom = numLociPerChrom,
                            sizesChrom = chromLengths)

subpop.founderGametes[0].flesh = createRandomFlesh(numLociPerChrom)
subpop.founderGametes[1].flesh = createRandomFlesh(numLociPerChrom)
subpop.founderGametes[2].flesh = createRandomFlesh(numLociPerChrom)
subpop.founderGametes[3].flesh = createRandomFlesh(numLociPerChrom)
subpop.genome = genome

now = time.time()
subpop.addFlesh()
#subpop.addFlesh(generation = pedigree.lastGen())
time.time() - now

#subpop.getInd(random.choice(pedigree.getId(pedigree.lastGen()))).gametes[0].flesh

G = subpop.genoMatrixDict()
len(G)
type(G.values())
G = np.vstack(G.values())



import bitmath
bitmath.Byte(G.nbytes).best_prefix()
Gp = np.packbits(G, axis = 1)
Gp.shape
bitmath.Byte(Gp.nbytes).best_prefix()
Gup = np.unpackbits(Gp, axis = 1)
Gup.shape

np.all(Gup == G) 

np.unpackbits(Gp, a).shape

bitmath.Byte(G.astype(np.uint8).nbytes).best_prefix()



np.unpackbits(np.packbits(np.array([1,0,1,0,1,0,0,1])))

bitmath.Byte(G.astype(np.bool_).nbytes).best_prefix()

len(gmd)
sys.getsizeof(gmd)


subpop.founderGametes


[subpop.realIBD(id) for id in pedigree.getId(pedigree.lastGen())]


#
subpop.individuals[random.choice(list(subpop.individuals.keys()))].gametes[0].skeleton







cS = crossoverSimulator(m = 4, p = 0.5, obligateChiasma = False)
cS.simulateCrossover(100)

        

        

# create founders

gametes = [FounderGamete([100, 100], 1), FounderGamete([100, 100], 2)]
founder1 = Individual(name = "founder1", gametes = gametes)
gametes = [FounderGamete([100, 100], 3), FounderGamete([100, 100], 4)]
founder2 = Individual(name = "founder2", gametes = gametes)
gametes = [FounderGamete([100, 100], 5), FounderGamete([100, 100], 6)]
founder3 = Individual(name = "founder3", gametes = gametes)


subpop = Subpopulation('subpop')
subpop.add(founder1)
subpop.add(founder2)
subpop.add(founder3)

progam = simulateMeiosis(parent, cS) 

now = time.time()
for i in range(10000):
    for i in range(2):
        cS.simulateCrossover(100)
time.time() - now

now = time.time()
for i in range(10000):
    np.random.poisson(4, 9)  
time.time() - now

progam.skeleton

import cProfile
cProfile.run('simulateMeiosis(parent, cS)')






















