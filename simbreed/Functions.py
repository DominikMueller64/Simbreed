# math
import numpy as np
import math
## tools
import itertools  # isSorted
import copy
import operator  # isSorted
import warnings
import collections  # namedtuple
from collections import Iterable
## random numbers
import uuid
import shortuuid
import random
#
from numba import jit

import matplotlib.pyplot as plt

# These are functions that should not need other package moduls.


# %%
def allUnique(x):
    """Check whether all elements are unique.

    Parameters
    ----------
    x : iterable
      Any iterable.

    Returns
    -------
    out : boolean
      Returns ``True`` is all elements in are unique and ``False`` if not.

    """
    seen = set()
    return not any(i in seen or seen.add(i) for i in x)
# %%

def getNewID_old(n=1, nChar=8):
    """Get a new UUID (Universally Unique Identifier).

    Parameters
    ----------
    n : natural number (default = 1)
      The number of IDs that should be generated.
    nChar : int (default = 8)
      The number of characters of the ID.
    Returns
    -------
    out : integer or list of integers
      | If n = 1, corresponding to the default, then a single ID is returned.
      | If n > 1, a list of IDs is returned.

    If nChar < 22, the ID is not guaranteed to be universally unique (UUID).
    However, for the default (nChar = 8), there are 57^8 =~ 100 Billon
    possible sequences, which still guarantees an essentially zero collision
    probability within a simulation session while greatly increasing
    readability of IDs for introspection.
    """
    def fnc():
        uid = shortuuid.uuid()
        if nChar < 1 or nChar > len(uid):
            return uid
        else:
            temp = np.random.choice([*uid], size = nChar, replace = False)
            return ''.join(temp.tolist())

    temp = [fnc() for _ in range(n)]
    return temp if n > 1 else temp[0]
# %%


@jit(nopython=True)
def getNewIDNumbaArray(n, digits):
    tmp = 10**(digits-1)
    ret = np.zeros(n, np.int64)
    for i in range(n):
        ret[i] = random.randint(tmp, 10 * tmp - 1)
    return ret

@jit(nopython=True)
def getNewIDNumbaSingle(digits):
    tmp = 10**(digits-1)
    return random.randint(tmp, 10 * tmp - 1)

def getNewID(n=1, digits=14, sequence=False):
    if n == 1:
        if sequence:
            return [getNewIDNumbaSingle(digits=digits)]
        return getNewIDNumbaSingle(digits=digits)
    else:
        return getNewIDNumbaArray(n=n, digits=digits).tolist()

#%timeit -n10000 -r10 for _ in range(5): getNewIDNumbaSingle(14)
#%timeit -n10000 -r10 getNewID(n=5)
#%timeit -n1000 -r10 getNewIDNumba(100, 14)

def isSorted(iterable, reverse=False, abs_tol=1e-6):
    """Check if an iterable is sorted, possibly reverse sorted.

    See the answer by hughdbrown
    `here <http://stackoverflow.com/questions/3755136/
    pythonic-way-to-check-if-a-list-is-sorted-or-not/>`_.

    Parameters
    ----------
    iterable : iterable
      An iterable.
    reverse : boolean (default = False)
      Should the iterable be checked for reverse sorting?
    abs_tol : float (default = 1e-6)
      The absolute tolerance used to check if values are distinct.

    Returns
    -------
    out : boolean
      If the iterable is (reverse) sorted ``True``, otherwise ``False``.

    """
    def pairwise(iterable):
        """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)

    def comp_func(x, y, reverse=reverse, abs_tol=abs_tol):
        my_operator = operator.ge if reverse else operator.le
        if math.isclose(x, y, abs_tol=abs_tol):
            raise ValueError("Values %f and %f are not distinct at a "
                             "tolerance level of %f." % (x, y, abs_tol))
        return my_operator(x, y)

    return all(comp_func(current_element, next_element)
               for current_element, next_element in pairwise(iterable))
# %%


def anyClose(iterable, abs_tol=1e-6):
    """Check if adjacent floats in an iterable are numerically identical

    Parameters
    ----------
    iterable : iterable
      An iterable.
    abs_tol : float (default = 1e-6)
      The absolute tolerance used to check if values are distinct.

    Returns
    -------
    out : boolean
      True if any pair is close, False otherwise.

    """
    def pairwise(iterable):
        """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)

    return any(math.isclose(x, y, abs_tol=abs_tol)
               for x, y in pairwise(iterable))
# %%%


def order(iterable, reverse=False):
    """Determine the order of an iterable, possibly reverse sorted.

    This is analogous to R's ``order()`` function.

    Parameters
    ----------
    iterable : iterable
      An iterable for which the order of the elements is sought.
      Must support indexing.
    reverse : boolean (default = False)
      Should the order be that of the reverse sorted iterable?

    Returns
    -------
    out : list
      List of integers with the (reverse) order of elements in ``iterable``.

    """
    return sorted(range(len(iterable)), key=lambda k: iterable[k],
                  reverse=reverse)
# %%


def reorder(iterable, order):
    """Reorder a list.

    Parameters
    ----------
    iterable : iterable
      An iterable that should be reordered.
      Must support indexing.
    order : list
      A list that contains the new order as permuted indices.
      Must support ``enumerate()``

    Returns
    -------
    out : list
      A reordered iterable of the same type as ``iterable``.

    """
    new = iterable.copy()
    for i, o in enumerate(order):
        new[i] = iterable[o]
    return new
# %%


def checkEqual(iterator):
    """ Check if all elements are identical. """
    return iterator is None or len(set(iterator)) == 1
# %%

def isSingle(obj):
    return isinstance(obj, str) or not isinstance(obj, collections.Iterable)
# %%

@jit(nopython=True)
def myBisect_left(a, x, lo=0):
    """ My pruned bisection algorithm. """
    hi = a.size
    while lo < hi:
        mid = (lo + hi) // 2
        if a[mid] < x:
            lo = mid + 1
        else:
            hi = mid
    return lo
# %%

#def cross(father, mother, crossoverSimulator, ID=None):
#    """Cross two individuals
#
#    Parameters
#    ----------
#    father : Individual
#      The father of the cross
#    mother : Individual
#      The mother of the cross
#    crossoverSimulator : CrossoverSimulator
#      An object of type 'CrossoverSimulator'
#    ID : string/integer (default = None)
#      A new identifier for the child. If 'None', a new UUID is obtained.
#
#    Return
#    ------
#    out : Individual
#      The child of the cross.
#
#    """
#    if ID is None:
#        ID = getNewID()
#    elif not isinstance(ID, (str, int)):
#        raise ValueError("If ID is provided, must be integer or string.")
#
#    fatherGamete = simulateMeiosis(father, crossoverSimulator)
#    motherGamete = simulateMeiosis(mother, crossoverSimulator)
#    return Individual(ID, father, mother, [fatherGamete, motherGamete])
## %%
#
#
#def selfing(parent, crossoverSimulator, ID=None):
#    """Self-fertilize an individual."""
#    return cross(parent, parent, crossoverSimulator, ID)
## %%
#
#
#def DH(parent, crossoverSimulator, ID=None):
#    """Produce a doubled haploid."""
#    if ID is None:
#        ID = getNewID()
#    elif not isinstance(ID, (str, int)):
#        raise ValueError("If ID is provided, must be integer or string.")
#
#    gamete = simulateMeiosis(parent, crossoverSimulator)
#    return Individual(ID, parent, parent, [gamete, gamete])
# %%

#def createRandomGenome(chromLengths, numLociPerChrom):
#    """Create a "random" genome
#
#    Parameters
#    ----------
#    chromLengths : list
#     A list with the lengths of chromosomes.
#    numLociPerChrom : list
#     A list with the number of loci per chromosome
#
#    Return
#    ------
#    out : Genome
#      A new genome
#
#    Loci will be uniformly distributed on the chromosomes.
#
#    """
#    numChrom = len(chromLengths)
#    if not numLociPerChrom == numChrom:
#        raise ValueError("'chromLengths' and 'numLociPerChrom' must have "
#                         "equal length.")
#    positions = [sorted(np.random.uniform(0.0, h, n)) for h, n in
#                 zip(chromLengths, numLociPerChrom)]
#    return Genome.Genome(chromLengths, positions)
# %%


def IDGenerationHandler(idGetter, ID, generation):
        if ID is None:
            if generation is None:
                ID=idGetter()  # Get all ID.
#                raise ValueError("Either 'ID' or 'generation' must be given.")
            else:
                ID = idGetter(generation)
#                ID.sort()
        else:
            if generation is not None:
                warnings.warn("'ID' was given, so 'Generation' is ignored")
            # TODO: replace by isSingle after writing unit tests.
            if isinstance(ID, (str, int)):
                ID = (ID,)
        return ID
# %%


def dualIDGenerationHandler(idGetter, ID, generation, IDb, generationB):
    ID = IDGenerationHandler(idGetter, ID, generation)
    if IDb is None and generationB is None:
        IDb = ID
    else:
        IDb = IDGenerationHandler(idGetter, IDb, generationB)
    rtuple = collections.namedtuple('rtuple', ['ID', 'IDb'])
    return rtuple(ID=ID, IDb=IDb)
# %%


def matrixComputer(idGetter, function,
                   ID=None, generation=None,
                   IDb=None, generationB=None,
                   symmetric=True, dtype=np.float_):
    rtuple = collections.namedtuple('rtuple', ['matrix', 'rownames', 'colnames'])
    ID = IDGenerationHandler(idGetter, ID, generation)
    if IDb is None and generationB is None:
        mat = np.zeros(shape=(len(ID), len(ID)), dtype=dtype)
        for ixa, ida in enumerate(ID):
            for ixb, idb in enumerate(ID[ixa + 1:], start=ixa + 1):
                temp = function(ida, idb)
                mat[ixa, ixb] = temp
                if symmetric:
                    mat[ixb, ixa] = temp
            mat[ixa, ixa] = function(ida)
        return rtuple(matrix=mat, rownames=ID, colnames=ID)
#        pkg = {'matrix': mat, 'rownames': ID, 'colnames':ID}
    else:
        IDb = IDGenerationHandler(idGetter, IDb, generationB)
        mat = np.zeros(shape=(len(ID), len(IDb)), dtype=np.float_)
        for ixa, ida in enumerate(ID):
            for ixb, idb in enumerate(IDb):
                mat[ixa, ixb] = function(ida, idb)
        return rtuple(matrix=mat, rownames=ID, colnames=IDb)
#        pkg = {'matrix': mat, 'rownames': ID, 'colnames': IDb}

#    return pkg


# %%

def simulateMeiosis(parent, crossoverSimulator):
    """Simulate Meiosis

    Parameters
    ----------
    parent : Individual
      An Individual for which meiosis should be simulated
    crossoverSimulator : CrossoverSimulator
      A CrossoverSimulator previously instantiated

    Return
    ------
    out : Gamete
      A Gamete as meiotic product

    """

    father = parent.gametes[0]; mother = parent.gametes[1]
    pat = father.skeleton; mat = mother.skeleton
    numChrom = len(mat)
    skeleton = [{} for _ in range(numChrom)]
    for iChrom in range(numChrom):
        matalle = mat[iChrom]['lineages']; matloc = mat[iChrom]['locations']
        patalle = pat[iChrom]['lineages']; patloc = pat[iChrom]['locations']
        L = matloc[-1]
        # simulate crossover locations; add -1 to the beginning
        Xlocations = crossoverSimulator.simulateCrossover(L)
#        product = [-1] + Xlocations
        product = np.hstack(([-1], Xlocations))
        cur_allele = int(random.random() > 0.5)  # much faster
        # first allele (0 : father, 1 : mother)
        parentStart = father if cur_allele == 0 else mother
        loc = list(); alle = list(); curpos = 0
        if len(product) == 1:
            if cur_allele == 0:
                skeleton[iChrom]['lineages'] = patalle
                skeleton[iChrom]['locations'] = patloc
            else:
                skeleton[iChrom]['lineages'] = matalle
                skeleton[iChrom]['locations'] = matloc
        else:
            for i in range(1, len(product)):  # Could use generator here.
                leftLoc = product[i - 1]; rightLoc = product[i]
                if cur_allele == 0:
                    parloc = patloc; paralle = patalle
                else:
                    parloc = matloc; paralle = matalle
                for j in range(len(parloc)):  # Could use generator here.
                    pl = parloc[j]
                    pa = paralle[j]
                    if leftLoc <= pl < rightLoc:
                        loc.append(pl); alle.append(pa)
                    elif rightLoc < pl:
                        break
                loc.append(rightLoc); alle.append(pa)
                cur_allele = 1 - cur_allele

            lastxo = product[-1]

            if cur_allele == 0:
                parloc = patloc; paralle = patalle
            else:
                parloc = matloc; paralle = matalle

            for j in range(len(parloc)):  # Could use generator here.
                if parloc[j] > lastxo:
                    loc.append(parloc[j]); alle.append(paralle[j])

            curpos = len(loc)
            if curpos > 1:  # clean up repeated alleles
                loc_clean = [None] * curpos; alle_clean = [None] * curpos
                loc_clean[0] = loc[0]; alle_clean[0] = alle[0]
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

        skeleton[iChrom]['Xlocations'] = Xlocations
        skeleton[iChrom]['parentStart'] = parentStart

    return skeleton


# %%

@jit(nopython=True)
def subroutineMeiosis(Xlocations, cur_allele, m, patloc, patalle,
                      matloc, matalle):
    loc = np.zeros(m, dtype=np.float64)
    alle = np.zeros(m, dtype=np.int64)
    curpos = 0
    ct = 0
    leftLoc = -1.0
    for i in range(len(Xlocations)):
        rightLoc = Xlocations[i]
        if cur_allele == 0:
            parloc = patloc; paralle = patalle
        else:
            parloc = matloc; paralle = matalle
        # TODO: Here is still large optimization potential. Store the last
        # position and only start searching from there.
        for j in range(len(parloc)):
            pl = parloc[j]
            pa = paralle[j]
            if leftLoc <= pl < rightLoc:
                loc[ct] = pl
                alle[ct] = pa
                ct += 1
            elif rightLoc < pl:
                break
        loc[ct] = rightLoc
        alle[ct] = pa
        cur_allele = 1 - cur_allele
        ct += 1
        leftLoc = rightLoc

    lastxo = Xlocations[-1]

    if cur_allele == 0:
        parloc = patloc; paralle = patalle
    else:
        parloc = matloc; paralle = matalle

    for j in range(len(parloc)):  # Could use generator here.
        if parloc[j] > lastxo:
            loc[ct] = parloc[j]
            alle[ct] = paralle[j]
            ct += 1

    curpos = ct
    if curpos > 1:  # clean up repeated alleles
        loc_clean = np.zeros(curpos, dtype=loc.dtype)
        alle_clean = np.zeros(curpos, dtype=alle.dtype)
        loc_clean[0] = loc[0]
        alle_clean[0] = alle[0]
        lastpos = 0
        for i in range(1, ct):
            if alle_clean[lastpos] == alle[i]:
                loc_clean[lastpos] = loc[i]
            else:
                lastpos += 1
                loc_clean[lastpos] = loc[i]
                alle_clean[lastpos] = alle[i]

        curpos = lastpos + 1
        return loc_clean[:curpos], alle_clean[:curpos]



def simulateMeiosis(parent, crossoverSimulator):
    """Simulate Meiosis

    Parameters
    ----------
    parent : Individual
      An Individual for which meiosis should be simulated
    crossoverSimulator : CrossoverSimulator
      A CrossoverSimulator previously instantiated

    Return
    ------
    out : Gamete
      A Gamete as meiotic product

    """

    father = parent.gametes[0]; mother = parent.gametes[1]
    pat = father.skeleton; mat = mother.skeleton
    numChrom = len(mat)
    skeleton = [{} for _ in range(numChrom)]
    for iChrom in range(numChrom):
            matalle = mat[iChrom]['lineages']; matloc = mat[iChrom]['locations']
            patalle = pat[iChrom]['lineages']; patloc = pat[iChrom]['locations']
            L = matloc[-1]
            # simulate crossover locations; add -1 to the beginning
            Xlocations = crossoverSimulator.simulateCrossover(L)
            cur_allele = int(random.random() > 0.5)  # much faster
            # first allele (0 : father, 1 : mother)
            parentStart = father if cur_allele == 0 else mother
            if len(Xlocations) == 0:
                if cur_allele == 0:
                    skeleton[iChrom]['lineages'] = patalle
                    skeleton[iChrom]['locations'] = patloc
                else:
                    skeleton[iChrom]['lineages'] = matalle
                    skeleton[iChrom]['locations'] = matloc
            else:
                m = len(Xlocations) + len(patloc) + len(matloc)
                loc, alle = subroutineMeiosis(Xlocations=Xlocations, cur_allele=cur_allele,
                                              m=m, patloc=patloc, patalle=patalle,
                                              matloc=matloc, matalle=matalle)
                skeleton[iChrom]['lineages'] = alle
                skeleton[iChrom]['locations'] = loc

            skeleton[iChrom]['Xlocations'] = Xlocations
            skeleton[iChrom]['parentStart'] = parentStart

    return skeleton

#    return Gamete(skeleton=skeleton, ID=getNewID(),
#                  father=father, mother=mother)
# %%

def nrow(x):
    """ Return number of rows of a 2D numpy array. """
    return x.shape[0]
# %%

def ncol(x):
    """ Return number of columns of a 2D numpy array. """
    return x.shape[1]
# %%

def isCyclic(graph):
            """Return True if the directed graph g has a cycle.
            g must be represented as a dictionary mapping vertices to
            iterables of neighbouring vertices. For example:

            >>> cyclic({1: (2,), 2: (3,), 3: (1,)})
            True
            >>> cyclic({1: (2,), 2: (3,), 3: (4,)})
            False

            """
            # http://codereview.stackexchange.com/questions/86021/check-if-a-directed
            #-graph-contains-a-cycle
    
            # The algorithm here is recursive, which might cause problems for very
            # large pedigrees. In this case, increase the recursion limit or
            # rewrite iteratively using a stack.
            path = set()
            visited = set()

            def visit(vertex):
                if vertex in visited:
                    return False
                visited.add(vertex)
                path.add(vertex)
                for neighbour in getattr(graph.get(vertex), 'parents', ()):
                    if neighbour in path or visit(neighbour):
                        return True
                path.remove(vertex)
                return False

            return any(visit(v) for v in graph)
# %%


def sort(graph, placeFounderFirst=True):
            """Sort the pedigree.

            If no data on generations are available, the algorithm provided in
            `Zhang et al. (2009) <http://www.medwelljournals.com/fulltext/
            ?doi=javaa.2009.177.182>`_ is applied to sort from scratch.

            """
            n = len(graph)  # We are sure that all IDs are unique.
            for v in graph.values():  # Set all generation data to zero.
                v.generation = 0
            ng = type(graph)()
            while graph:
                shc = copy.copy(graph)
                for v1 in shc.values():
                    for v2 in shc.values():
                        if not v1 is v2:
                            if v1.ID in v2.parents:
                                v1.generation += 1
                                break
                    else:
                        ng.update({v1.ID: graph.pop(v1.ID)})

#            for uu in ng.values():print(uu.generation)
#            maxgen = next(reversed(ng.values())).generation
            maxgen = max(ind.generation for ind in ng.values())
            for v in ng.values():
                if placeFounderFirst and v.founder:
                    v.generation = 0
                else:
                    v.generation = abs(v.generation - maxgen)
            graph.update(ng)
#            graph.update(sorted(ng.items(), key = lambda x: x[1].generation))
# %%


def plotMatrix(X):
    plt.imshow(X=X, cmap='YlOrRd', interpolation='nearest')
# %%


@jit(nopython=True)
def pearson(x, y):
    mx = np.mean(x)
    my = np.mean(y)
    xy = yy = xx = 0.0
    for i in range(len(x)):
        xic = x[i] - mx
        yic = y[i] - my
        xy += xic * yic
        yy += yic**2
        xx += xic**2
    if xx == 0.0 or yy == 0.0:
        return np.nan
    return xy / math.sqrt(xx * yy)


@jit(nopython=True)
def helper(x, y, ssq1, ssq2):
    if ssq1 == 0.0 or ssq2 == 0.0:
        return np.nan
    else:
        return np.dot(x, y) / math.sqrt(ssq1 * ssq2)




@jit(nopython=True)
def LPS_routine(chromA, chromB, pos, distMin, distMax):
    corA = np.zeros(len(pos) * (len(pos) - 1) // 2)
    corB = np.zeros_like(corA)
    ct = 0

    for ix1 in range(len(pos) - 1):
        p1 = pos[ix1]
        ix2 = ix1 + 1
        p2 = pos[ix2]
        lim2 = p1 + distMax
        while p2 <= lim2:
            if distMin <= p2 - p1:
                corA[ct] = pearson(x=chromA[:, ix1], y=chromA[:, ix2])
                corB[ct] = pearson(x=chromB[:, ix1], y=chromB[:, ix2])
                ct += 1
            ix2 += 1
            p2 = pos[ix2]

#    chromA = chromA - chromA.mean(axis=0)
#    chromB = chromA - chromB.mean(axis=0)
#
#    ssqA = np.zeros(ncol)
#    ssqB = np.zeros_like(ssqA)
##    for j in range(ncol):
##        tmp = chromA[:, j]
##        ssqA[j] = np.dot(tmp, tmp)
##        tmp = chromB[:, j]
##        ssqB[j] = np.dot(tmp, tmp)
##
##    for ix1 in range(len(pos) - 1):
##        p1 = pos[ix1]
##        ix2 = ix1 + 1
##        p2 = pos[ix2]
##        lim2 = p1 + distMax
##        while p2 <= lim2:
##            if distMin <= p2 - p1:
##                corA[ct] = helper(x=chromA[:, ix1], y=chromA[:, ix2], ssq1=ssqA[ix1], ssq2=ssqA[ix2])
##                corB[ct] = helper(x=chromB[:, ix1], y=chromB[:, ix2], ssq1=ssqB[ix1], ssq2=ssqB[ix2])
##                ct += 1
##            ix2 += 1
##            try:
##                p2 = pos[ix2]
##            except IndexError:
##                break
#
#
#
    return corA[:ct], corB[:ct]
# %%%

@jit(nopython=True)
def LD_r2_routine(chrom, pos, distMin, distMax):
    cor = np.zeros(len(pos) * (len(pos) - 1) // 2)
    ct = 0
    for ix1 in range(len(pos) - 1):
        p1 = pos[ix1]
        ix2 = ix1 + 1
        p2 = pos[ix2]
        lim2 = p1 + distMax
        while p2 <= lim2:
            if distMin <= p2 - p1:
                cor[ct] = pearson(x=chrom[:, ix1], y=chrom[:, ix2])
                ct += 1
            ix2 += 1
            p2 = pos[ix2]

    return cor[:ct]
# %%

@jit(nopython=True)
def LD_routine(chrom, p, pos, distMin, distMax, metric=1):
    LD = np.zeros(len(pos) * (len(pos) - 1) // 2)
    ct = 0
    for ix1 in range(len(pos) - 1):
        p1 = pos[ix1]
        x1 = chrom[:, ix1]
        f1 = p[ix1]

        ix2 = ix1 + 1
        p2 = pos[ix2]
        lim2 = p1 + distMax
        while p2 <= lim2:
            if distMin <= p2 - p1:
                x2 = chrom[:, ix2]
                f2 = p[ix2]
                # Avoids checking for a specific allele coding.
                LD[ct] = np.mean(x1 == x2) / 2.0 - f1 * f2 - 0.5 * (1.0 - f1 - f2)

#                s = 0
#                for i in range(len(x1)):
#                    if x1[i] == x2[i] == 1:
#                        s += 1
#                D[ct] = s / len(x1) - f1 * f2

                if metric in (2, 3):
                    # Dcor
                    if metric == 2:
                        div = math.sqrt(f1 * (1.0 - f1) * f2 * (1.0 - f2))
                    # Dprime 
                    elif metric == 3:
                        if LD[ct] < 0.0:
                            # -min produces absolute values in numerator.
                            div = -min(f1 * f2, (1.0 - f1) * (1.0 - f2))
                        elif LD[ct] > 0.0:
                            div = min(f1 * (1.0 - f2), (1.0 - f1) * f2)

                    if div == 0.0:
                        LD[ct] = np.nan
                    else:
                        LD[ct] /= div


                ct += 1
            ix2 += 1
            p2 = pos[ix2]

    return LD[:ct]
# %%

@jit(nopython=True)
def no_unique(x):
    if estimate is None:
        z = np.empty_like(x)

    ct = 1
    z[0] = x[0]
    for i in range(1, len(x)):
        c = x[i]
        for j in range(ct):
            if c == z[j]:
                break
        else:
            z[ct] = c
            ct += 1

    return ct
# %%

@jit(nopython=True)
def SMC_routine(X):
    m = X.shape[0] // 2           
    GRM = np.eye(m)
    for i in range(m - 1):
        ixi1 = i * 2
        ixi2 = i * 2 + 1
        for j in range(i + 1, m):
            ixj1 = j * 2
            ixj2 = j * 2 + 1
            
            
            
            GRM[i, j] = GRM[j, i] = np.mean(X[ixi1:ixi2] == X[ixj1:ixj2])

    return GRM






