# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 15:20:12 2016

@author: domi89
"""

import numpy as np
import sys
from collections import Counter
import pandas as pd
import bisect
import string
import random
import math

def isstring(s):
    # if we use Python 3
    if (sys.version_info[0] == 3):
        return isinstance(s, str)
    # we use Python 2
    return isinstance(s, basestring)


def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


# this should be large enought to create truely unique labels
ID_HAPLOTYPE_LENGTH = 30
     
# the object argument is necessary for Python 2.x, not for Python 3.x   
class chromosome(object):

    def __init__(self, alleles):
        self.alleles = alleles
    
    def get_chromosome(self):
#        print("GETTER")
        return self._alleles
        
    def set_chromosome(self, value):
#        print("SETTER")
        self._alleles = value
        self.nLoci = len(value)
    
    alleles = property(get_chromosome, set_chromosome, None, 'Chromosome')
    
    def __repr__(self):
        return repr(self.alleles)
        
    def __len__(self):
        return self.nLoci
    
    def __getitem__(self, key):
        return self.alleles[key]
    
#chr = chromosome(np.array([1,2,3]))
#chr.alleles
#chr.nLoci
#chr.alleles = np.array([5,6,4,3,4,5])
#chr.alleles
#len(chr)
#chr
        
class haplotype(object):

    def __init__(self, haplo, id = None, motherId = None, fatherId = None):
        self.haplo = haplo
        self.id = id
        self.motherId = motherId
        self.fatherId = fatherId

    def get_haplo(self):
        #print("GETTER")
        return self._haplo

    def set_haplo(self, value):
        #print("SETTER")
        self._haplo = value
        self.nChr = len(value)
        self.nLociPerChr = [len(ch) for ch in value]
        self.nLoci = sum(self.nLociPerChr)
        

    haplo = property(get_haplo, set_haplo, None, 'Haplotype')
        
    def __repr__(self):
        return repr(self.haplo)
        
        
    def __getitem__(self, key):
        return self.haplo[key]


class genotype(object):
    
    def __init__(self, geno, id = None):
        self.geno = geno
        self.id = id
    
    def get_geno(self):
        return self._geno
    
    def set_geno(self, value):
        self._geno = value
        self.ploidy = len(value)
        self.nLoci = value[0].nLoci
        self.nChr = value[0].nChr
    
    geno = property(get_geno, set_geno, None, 'Genotype')
    
    def __repr__(self):
        return repr(self.geno)
    
    def __getitem__(self, key):
        return self.geno[key]
        
        
    def stack(self):
        tmp = [None for _ in range(self.nChr)]
        for i_chr in range(self.nChr):
            tmp[i_chr] = np.vstack([self[i_haplo][i_chr].alleles for i_haplo in range(self.ploidy)])
        return tmp
        
    def count_alleles(self):
        return [np.apply_along_axis(lambda x: len(set(x)), 0, ch) for ch in self.stack()]
                
    
    def which_polymorphic(self):
        return [ch > 1 for ch in self.count_alleles()]
        
    def recombine(self, gM):
        
        gamete = [None for _ in range(self.nChr)]
        for i in range(self.nChr):         
            ch = [haplo[i] for haplo in self.geno] #!!!            
            sgM = mymap[i] #!!!
            pos = sgM.pos.tolist()
            minPos = min(pos)
            maxPos = max(pos)
            ln = maxPos - minPos
            nExpRec = ln * 1/100
            nRec = np.random.poisson(nExpRec, 1)[0]
            chiaPos = sorted(np.random.uniform(low = minPos, high = maxPos, size = nRec).tolist())
            
            a,b = ch[0].alleles.copy(), ch[1].alleles.copy()
            for chia in chiaPos:
                idx = bisect.bisect_left(pos, chia)
                # rd = np.random.permutation(range(g1.ploidy)).tolist()#!!!
                an = np.concatenate((a[:idx], b[idx:]))
                b = np.concatenate((b[:idx], a[idx:]))
                a = an
                
            coin = np.random.choice(self.ploidy)
            if coin == 0:
                gam = a
            else:
                gam = b
                
            gamete[i] = gam
        
        return haplotype(gamete, 
                         id = id_generator(size = ID_HAPLOTYPE_LENGTH),
                         motherId = self[0].id, 
                         fatherId = self[1].id)
        
        


class subpopulation(object):
    
    def __init__(self, pop, id = None, gM = None):
        self.pop = pop
        self.id = id
        self.gM = gM
        
    def get_pop(self):
        return self._pop
    
    def set_pop(self, value):
        self._pop = value
        self.ploidy = value[0].ploidy
        self.nLoci = value[0][0].nLoci
        self.nChr = value[0].nChr
        self.nGeno = len(value)
        
    pop = property(get_pop, set_pop, None, 'Subpopulation')
    
    def __repr__(self):
        pass
    
    def size(self):
        return len(self.pop)
    
    
    def get_ids(self):
        return [_.id for _ in self.pop]
    
    def __getitem__(self, key):
        if isinstance(key, int):
            return self.pop[key]
        else:
            return self.pop[self.get_ids().index(key)]
            
    
    def stack(self):
        return [np.vstack([self[i_geno].stack()[i_chr] for i_geno in range(self.nGeno)]) for i_chr in range(self.nChr)]
    

    def count_alleles(self):
        # return [np.apply_along_axis(lambda x: len(set(x)), 0, ch) for ch in self.stack()]
        return [[len (__) for __ in _] for _ in self.unique_alleles()]
  
    def unique_alleles(self):
        return [[list(set(ch[..., j])) for j in range(ch.shape[1])] for ch in self.stack()]
        
    def allele_frequency(self):
        return [[dict(Counter(ch[..., j]))  for j in range(ch.shape[1])] for ch in self.stack()]
        
    def reproduce(self, mateOps):
        
        sp1.


class matingOperator(object):
    def __init__(self, selectedFraction = 1.0)
        self.selectedFraction = selectedFraction

class randomMating(matingOperator):
    
    def __init__(self, allowSelfing = True):
        self.allowSelfing = allowSelfing
        
    def get_matings(self, subpop, newPopSize):
        nSelected = math.ceil(self.selectedFraction * subpop.size())
        selected = np.random.choice(a = subpop.get_ids(), size = nSelected, replace = False).tolist()

        if self.allowSelfing:
                matings = np.random.choice(a = selected, size = 2 * newPopSize, replace = True).reshape((newPopSize, 2)).tolist()
            
        else:
                # Only possible if population size > 1
                matings = [None for _ in range(newPopSize)]
                for i in range(newPopSize):
                    first = np.random.choice(a = selected)
                    tmp = c=selected.copy()
                    tmp.remove(first)
                    second = np.random.choice(a = tmp)
                    matings[i] = [first, second]
        return matings
        
        
        
        
        
        
        
class matingOperator(object):
        
    def __init__(self, mode, selectedFraction = 1.0, allowSelfing = True):
                     
        modes = ['randomMating', 'inbreeding', 'roundRobin', 'synthetic', 'forcedOutbreeding']  
        if mode not in modes: 
            raise ValueError('Invalid mode. Expecting one of: %s' % modes)
                     
        
        self.selectedFraction = selectedFraction
        self.randomMating = randomMating
        self.allowSelfing = allowSelfing
    
    def get_matings(self, subpop, newPopSize):
        nSelected = math.ceil(self.selectedFraction * subpop.size())
#        print(nSelected)
        # Selected individuals that reproduce
        selected = np.random.choice(a = subpop.get_ids(), size = nSelected, replace = False).tolist()        
#        print(str(selected))
        if self.mode == 'randomMating':
            
            if self.allowSelfing:
                matings = np.random.choice(a = selected, size = 2 * newPopSize, replace = True).reshape((newPopSize, 2)).tolist()
            
            else:
                # Only possible if population size > 1
                matings = [None for _ in range(newPopSize)]
                for i in range(newPopSize):
                    first = np.random.choice(a = selected)
                    tmp = c=selected.copy()
                    tmp.remove(first)
                    second = np.random.choice(a = tmp)
                    matings[i] = [first, second]
                    
        if self.mode == 'inbreeding':
            matings = [[_,_] for _ in selected]
            
        if self.mode = 'roundRobin':
            matings = [(selected + selected[0:1])[i:(i+2)] for i in range(nSelected)]
            
        if self.mode == 'forcedOutbreeding':
            print("forcedOutbreeding might be implemented in future")
            return None
                
                
                    
        return matings
            
    
        
    
#    def subset(self, **kwargs) :
#        if 'genotype' in kwargs.keys():
                    
    
    
    
class genMap(object):
    
    def __init__(self, gM):
        self.gM = gM
        
        
    def get_gM(self):
        return self._gM
    
    def set_gM(self, value):
        self._gM = value
        self.nChr = len(value)
        self.nLociPerChr = [ch.pos.count() for ch in value]
        self.nLoci = sum(self.nLociPerChr)
        
        
    gM = property(get_gM, set_gM, None, 'Subpopulation')
    
    def concat(self):
        return pd.concat(self.gM)
        
    def __getitem__(self, key):
        return self.gM[key]
        
    def subset(self, **kwargs):
        tmp = mymap.concat()
        if 'id' in kwargs.keys():
            tmp = tmp.loc[kwargs['id'], ]
        if 'chr' in kwargs.keys():
            tmp = tmp[tmp.chr.isin(kwargs['chr'])]
            
        self.gM = [tmp[tmp.chr == ch] for ch in set(tmp.chr)]
        
#    def list_indices



import random
import timeit
import time

nLoci = 5
nChr = 3

pos = np.linspace(start = 0, stop = 100, num = nLoci)

gM = [pd.DataFrame({"chr" : pd.Series([i + 1] * nLoci, 
                                      index = list(map('L{}'.format, range(i * nLoci + 1, (i + 1) * nLoci + 1)))),
                    "pos" : pd.Series(pos, index = list(map('L{}'.format, range(i * nLoci + 1, (i + 1) * nLoci + 1))))})
                    for i in range(nChr)]
                        
mymap = genMap(gM)

tmp = ['A', 'B', 'C', 'D']
h1 = haplotype([chromosome(np.array([random.choice(tmp) for _ in range(nLoci)])) for _ in range(nChr)], id = "asldkjf")
h2 = haplotype([chromosome(np.array([random.choice(tmp) for _ in range(nLoci)])) for _ in range(nChr)], id = "lasjdf")

g1 = genotype([h1, h2], "G1")

gamete1 = g1.recombine(mymap)
gamete2 = g1.recombine(mymap)

g2 = genotype([h1, h2], "G2")
g3 = genotype([h1, h2], "G3")
sp = subpopulation([g1, g2, g3])

matOps = matingOperator(allowSelfing = False, selectedFraction = 0.65)
matOps.get_matings(sp, 10)



sp1.get_ids()
sp1.nGeno
sp1["G2"]
sp1.stack()
sp1.unique_alleles()
sp1.count_alleles()
sp1.allele_frequency()


[ch for ch in sp1.unique_alleles()]
sp1.count_alleles()

    for i_chr in range(sp1.nChr):       
        tmp = [sp1[i_geno][i_haplo][i_chr].alleles for i_haplo in range(sp1.ploidy) for i_geno in range(sp1.nGeno)]
    
    
    





c1 = chromosome(np.array(['A', 'B','A', 'B']))
c2 = chromosome(np.array(['B', 'A','A', 'A']))
h1 = haplotype([c1, c2])
c3 = chromosome(np.array(['A', 'D','B', 'B']))
c4 = chromosome(np.array(['C', 'A','D', 'D']))
h2 = haplotype([c3, c4])

g = genotype([h1, h2])
g.nLoci
g.ploidy
g.nChr

allele_count = g.count_alleles()

[g.geno[0][i] for i in range(g.nChr)]

g.geno[0].haplo[i].alleles
g.geno[1].haplo[i].alleles




    






        
      