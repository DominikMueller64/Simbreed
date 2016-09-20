# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 20:19:32 2016

@author: dominik
"""

import scipy.optimize
import numpy as np
import numpy.matlib
import bisect
import math
import uuid
import random
# import time
import collections
from collections.abc import Mapping, Sequence  # Pedigree, Population
# import heapq
import operator  # isSorted
import itertools  # isSorted
# import sys
import copy
import pandas as pd
import warnings

"""Demonstrate high quality docstrings.

Module-level docstrings appear as the first "statement" in a module. Remember,
that while strings are regular Python statements, comments are not, so an
inline comment may precede the module-level docstring.

After importing a module, you can access this special string object through the
``__doc__`` attribute; yes, it's actually available as a runtime attribute,
despite not being given an explicit name! The ``__doc__`` attribute is also
what is rendered when you call ``help()`` on a module, or really any other
object in Python.

You can also document a package using the module-level docstring in the
package's ``__init__.py`` file.

"""


# def average(a, b):
#    """
#    Given two numbers a and b, return their average value.
#
#    Parameters
#    ----------
#    a : number
#      A number
#    b : number
#      Another number
#
#    Returns
#    -------
#    res : number
#      The average of a and b, computed using 0.5*(a + b)
#
#    Example
#    -------
#    >>> average(5, 10)
#    7.5
#
#    """
#
#    return (a + b) * 0.5

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


def getNewID(n=1):
    """Get a new UUID (Universally Unique Identifier).

    Parameters
    ----------
    n : natural number (default = 1)
      The number of UUID that should be generated.

    Returns
    -------
    out : integer or list of integers
      | If n = 1, corresponding to the default, then a single UUID in integer
        format is returned.
      | If n > 1, a list of UUID in integer format is returned.

    """
    if n == 1:
        return uuid.uuid4().int
    else:
        return [uuid.uuid4().int for _ in range(n)]
# %%


def dpoisson(k, lambda_):
    """Density function of the Poisson distribution.

    The Poisson distribution is the limit of the binomial distribution
    for large N.

    Parameters
    ----------
    k : non-negative integer
      The non-negative integer for which the discrete density function should
      be evaluated
    lambda_ : positive real number
      The expectation/variance.

    Returns
    -------
    out : float
      Value of the Possion density function.

    """
    return lambda_**k * math.exp(-lambda_) / math.factorial(k)
# %%


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
# %%


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


class Genome(object):

    """Genome object

    Instance variables
    ------------------
    chromLengths : A list with chromosome lengths in **Centimorgan**. This is
    everything that is required for starting to simulate.

    positions : A list containing lists of genetic positions of loci. This is
    only required for creating genotypes after simulation.

    locusNames : A list containing lists of names of loci.

    numChrom : The number of chromosomes.

    Note that ``chromLengths`` can only be set once during instantiation, but
    not later on, as opposed to ``positions`` and ``locusNames``. If the length
    of chromosomes changes, a new ``Genome`` object must be instantiated.

    Note than ``LocusNames`` can only be specified if ``positions`` are given.

    Note that all chromosomes start at 0.0 Centrimorgan.

    """

    def __init__(self, chromLengths, positions=None, locusNames=None):
        """Create a new genome

        Parameters
        ----------
        chromLengths : list of floats
          A list with chromosome lengths in **Centimorgan**. Practically, these
          can be the positions of the respective last locus on the chromosome,
          pluse some delta (e.g., 1e-6).
        positions : list of lists of floats (default = None)
          A list containing lists of **sorted** genetic positions of loci on
          the chromosomes.
        locusNames : list of lists of integers or strings (default = None)
          A list containing lists of names of loci on the chromosomes in the
          same order as ``positions``.
        """

        # check chromLengths
        for ix, cl in enumerate(chromLengths):
            if cl <= 0.0:
                raise ValueError("Chromosome %d has nonpositve length (%f)."
                                 % (ix, cl))
            elif cl < 10.0:
                warnings.warn("Chromosome %d has only a length of %f "
                              "Centimorgan" % (ix, cl))

        self.__numChrom = len(chromLengths)
        self.__chromLengths = chromLengths
        self.positions = positions
        self.locusNames = locusNames

    @property
    def numChrom(self):
        return self.__numChrom

    @property
    def chromLengths(self):
        return self.__chromLengths

    @property
    def positions(self):
        return self.__positions

    @positions.setter
    def positions(self, value):
        abs_tol = 1e-6

        if not len(value) == self.__numChrom:
            raise ValueError(
                "The length of 'positions' (%d) does not match the number of "
                "chromosomes (%d)." % (len(value), self.__numChrom))

        for ix, pos in enumerate(value):
            if not isSorted(pos, abs_tol=0.0):
                raise ValueError("Positions on chromosome %d are not sorted "
                                 % (ix))
            # We can now assume pos to be sorted
            if pos[0] < 0.0:
                raise ValueError("The first position on chromosome %d (%f) "
                                 "must be >= 0.0." % (ix, pos[0]))

            if pos[-1] > self.__chromLengths[ix]:
                raise ValueError("The last position on chromosome %d (%f), "
                                 "must be <= %f."
                                 % (ix, pos[0], self.__chromLengths[ix]))

            if anyClose(pos, abs_tol=abs_tol):
                raise ValueError("Adjacent positions on chromosome %d are "
                                 "too close (tolerance = %f)" % (ix, abs_tol))

        self.__positions = value

    @property
    def locusNames(self):
        return self.__locusNames

    @locusNames.setter
    def locusNames(self, value):
        if value is None:
            self.__locusNames = value
            return None  # exit

        if not self.__positions:
            raise AttributeError("Before setting 'locusNames', 'positions' "
                                 "must be specified.")

        if not len(value) == self.__numChrom:
            raise ValueError(
                "The length of 'locusNames' (%d) does not match the number of "
                "chromosomes (%d)." % (len(value), self.__numChrom))

        for ix, names in enumerate(value):
            if not len(names) == len(self.positions[ix]):
                raise ValueError(
                    "The number of loci names (%d) on chromosome %d does not "
                    "match the number of positions (%d)."
                    % (len(names), ix, len(self.positions[ix])))

        if not allUnique((name for chrom in value for name in chrom)):
            raise ValueError("There are duplicated locus names.")

        self.__locusNames = value
# %%


class Gamete:

    """Gamete

    Instance variables
    ------------------
    skeleton : structure of the genome with reference to founder gametes

    ID : unique identifier

    father : father gamete

    mother : mother gamete

    flesh : genotype with reference to founder alleles

    Note that the words 'father' and 'mother' are only used to distinguish
    the two parental gametes.

    """
    # TODO: Check how the concept of an parent attribute for gametes can be
    # useful for later computaions
    def __init__(self, skeleton, ID=None, father=None,
                 mother=None, flesh=None, parent=None):
        """Instantiate a new Gamete

        Parameters
        ----------
        skeleton : list of dictionaries
          A list containing four-element dictionaries for each chromosome.
          One key is 'lineages' and refers to a list of integers indicating
          the origin of chromosomal segments. The second key is 'locations'
          and refers to a list of genetic positions, indicating the start of
          a new chromosomal segment. The fourth key is 'Xlocations' and harbors
          the positions of the crossover events that generated the chromosome.
          The last key is 'parentalStart' and tells from which parental gamete
          the head of the chromosome originates.
        ID : integer of string
          A unique identifier.
        father : Gamete
          First gamete from which this gamete was meiotically derived.
        mother : Gamete
          Second gamete from which this gamete was meiotically derived.

        Usually the user does not need to specify the 'skeleton', because this
        automatically happens when new founder individuals are added to the
        population. Also, 'father' and 'mother' do not need to be set.

        The user only needs to set 'flesh' for gametes of founder individuals
        (founder gametes) if genotypic data should be generated for the
        progeny. This can be done after the simulate is finished.

        """
        if not ID:
            ID = getNewID()  # each gamete get a UUID

        self.ID = ID
        self.father = father
        self.mother = mother
        self.skeleton = skeleton
        self.flesh = flesh
        self.parent = parent

    @property
    def ID(self):
        return self.__ID

    @ID.setter
    def ID(self, value):
        self.__ID = value

    @property
    def father(self):
        return self.__father

    @father.setter
    def father(self, value):
        if value is None:
            self.__father = value
            return None
        elif not isinstance(value, Gamete):
            raise ValueError("Attribute 'father' must be a 'Gamete'.")
        self.__father = value

    @property
    def mother(self):
        return self.__mother

    @mother.setter
    def mother(self, value):
        if value is None:
            self.__mother = value
            return None
        elif not isinstance(value, Gamete):
            raise ValueError("Attribute 'mother' must be a 'Gamete'.")
        self.__mother = value

    @property
    def skeleton(self):
        return self.__skeleton

    @skeleton.setter
    def skeleton(self, value):
        for ix, chrom in enumerate(value):
            if not set(chrom.keys()).issuperset(['lineages', 'locations',
                                                 'Xlocations', 'parentStart']):
                raise ValueError("The provided 'skeleton' on chromosome %d "
                                 "does not contain the necessary keys." % (ix))
        self.__skeleton = value

    @property
    def flesh(self):
        return self.__flesh

    @flesh.setter
    def flesh(self, value):
        # Check if flesh fits to the skeletton.
        # This might be a later performance bottleneck and direct access to
        # __flesh should be pursued in such circumstances.
        if value is None:
            self.__flesh = value
            return None
        elif not len(value) == len(self.skeleton):
            raise ValueError("The length of 'flesh' (%d) does not match the "
                             "length of 'skeleton' (%d) in gamete %s."
                             % (len(value), len(self.skeleton), str(self.ID)))

        for ix in range(len(value)):
            if not len(value[ix]) == len(self.skeleton[ix]):
                raise ValueError("The number of loci of the flesh does not "
                                 "match the number of loci of the skeleton on "
                                 "chromosome %d." % (ix))

        self.__flesh = value

    @property
    def parent(self):
        return self.__parent

    @parent.setter
    def parent(self, value):
        if value is None:
            self.__parent = value
            return None
        elif not isinstance(value, Individual):
            raise ValueError("Attribute 'parent' must be an 'Individual'.")
        self.__parent = value
# %%


class FounderGamete(Gamete):

    """FounderGamete

    Subclasses ``Gamete``

    Instance variables
    ------------------
    lineage : lineage of the founder gamete

    The user does not need to manually instantiate new FounderGametes.

    """

    def __init__(self, chromLengths, lineage, ID=None, flesh=None):
        """Instantiate a new FounderGamete

        Parameters
        ----------
        lineage : integer
          The lineage of a founder gamete is an integer identifier labelling
          the gamete. This is used to track the origin of chromosomal
          segments in later generations to founder gametes.

        """
        skeleton = [{'lineages': [lineage],
                     'locations': [chromLen],
                     'Xlocations': [None],
                     'parentStart': None} for chromLen in chromLengths]

        Gamete.__init__(self, ID=ID, skeleton=skeleton, flesh=flesh)

        self.__lineage = lineage

    @property
    def lineage(self):
        return self.__lineage
# %%


class Individual:

    """Individual

    Instance variables
    ------------------
    ID : Identifier

    father : father individual

    motehr : mother individual

    gametes : gametes carried by the individual


    """

    def __init__(self, ID=None, father=None, mother=None, gametes=None):
        """Instantiate an Individual

        Parameters
        ---------
        ID : integer or string
          A unique identifier
        father : Individual
          The father of the individual
        mother : Individual
          The mother of the individual
        gametes : list
          A list with the two Gametes carried by the individual


        """
        self.ID = ID
        self.father = father
        self.mother = mother
        self.gametes = gametes

        @property
        def ID(self):
            return self.__ID

        @ID.setter
        def ID(self, value):
            self.__ID = value

        @property
        def father(self):
            return self.__father

        @father.setter
        def father(self, value):
            if value is None:
                self.__father = value
                return None
            elif not isinstance(value, Individual):
                raise ValueError("Attribute 'father' must be "
                                 "an 'Individual'.")
            self.__father = value

        @property
        def mother(self):
            return self.__mother

        @mother.setter
        def mother(self, value):
            if value is None:
                self.__mother = value
                return None
            elif not isinstance(value, Individual):
                raise ValueError("Attribute 'mother' must be "
                                 "an 'Individual'.")
            self.__mother = value

# %%


class Population(Mapping):

    """Population

    Instance variables
    ------------------
    name : population name

    founderGametes : gametes that founded the population

    individuals : individuals of the population

    founderIndividuals: individuals that founded the population

    genome : Genome object

    pedigree : Pedigree object

    crossoverSimulator : CrossoverSimulator object

    Instance methods
    ----------------

    """

    def __init__(self, crossoverSimulator, genome=None, ID=None):
        """Instantiate a new Population

        Parameters
        ----------
        crossoverSimulator : CrossoverSimulator
          Used to simulate crossover events.

        genome: Genome
          A Genome that has at least chromosome lengths specified.

        name : integer or string
          Name or identifier of the population


        """
        self.__ID = ID
        self.crossoverSimulator = crossoverSimulator
        self.genome = genome
        self.__pedigree = Pedigree([], [], [])  # Pedigree is empty
        self.__individuals = dict()
        self.__founderIndividuals = dict()
        self.__gametes = dict()
        self.__founderGametes = dict()

    # ID
    @property
    def ID(self):
        return self.__ID

    # crossoverSimulator
    @property
    def crossoverSimulator(self):
        return self.__crossoverSimulator

    @crossoverSimulator.setter
    def crossoverSimulator(self, value):
        if not isinstance(value, CrossoverSimulator):
            raise ValueError("A valid CrossoverSimulator must be provided.")
        self.__crossoverSimulator = value

    # genome
    @property
    def genome(self):
        return self.__genome

    @genome.setter
    def genome(self, value):
        if not isinstance(value, Genome):
            raise ValueError("A valid Genome must be provided.")

        if not value.chromLengths:
            raise ValueError("The genome must at minimum have chromLengths.")

        self.__genome = value

    # pedigree
    @property
    def pedigree(self):
        return self.__pedigree

    # individuals
    @property
    def individuals(self):
        return self.__individuals

    # founderIndividuals
    @property
    def founderIndividuals(self):
        return self.__founderIndividuals

    # gametes
    @property
    def gametes(self):
        return self.__gametes

    # founderGametes
    @property
    def founderGametes(self):
        return self.__founderGametes

    def __getitem__(self, item):
        return self.__individuals[item]

    def __len__(self, item):
        return self.__individuals.__len__

    def __iter__(self):
        return self.__individuals.items().__iter__()

    def _addIndividual(self, ind):

        if not isinstance(ind, Individual):
            raise ValueError("Must be of class 'Individual'.")

        if ind.ID in self.individuals:
            raise ValueError("An individual with the same ID is "
                             "already present.")

        fatherGamete = ind.gametes[0]
        motherGamete = ind.gametes[1]
        for gamete in (fatherGamete, motherGamete):
            if gamete.ID in self.gametes:
                raise ValueError("A gamete with the same ID is "
                                 "already present.")

        isFatherFounder = isinstance(fatherGamete, FounderGamete)
        isMotherFounder = isinstance(motherGamete, FounderGamete)
        if isFatherFounder or isMotherFounder:
            self.founderIndividuals[ind.ID] = ind

        if isFatherFounder:
            self.founderGametes[fatherGamete.ID] = fatherGamete

        if isMotherFounder:
            self.founderGametes[motherGamete.ID] = motherGamete

        self.individuals[ind.ID] = ind
        self.gametes[fatherGamete.ID] = fatherGamete
        self.gametes[motherGamete.ID] = motherGamete

    def addFounders(self, nFounder, ID=None, founderInbreeding=None):
        """ Add founders to the population and the pedigree

        Parameters
        ----------
        nFounder : integer
          The number of founder individuals that should be added

        ID : list (default = None)
          Unique identifiers for the founders (integer or strings)

        """
        if ID is None:
            ID = [getNewID() for _ in range(nFounder)]
        elif nFounder != len(ID):
            raise ValueError("If the optional argument 'ID' is provided, "
                             "it must have the correct length.")
        elif not allUnique(ID):
            raise ValueError("If the optimal argument 'ID' is provided, "
                             "all IDs must be unique.")
        if founderInbreeding is None:
            founderInbreeding = dict.fromkeys(ID, 0.0)
        
        self.pedigree.ID.extend(ID)
        self.pedigree.father.extend([None] * nFounder)
        self.pedigree.father.extend([None] * nFounder)
#        self.pedigree.generation.extend([0] * nFounder)  # per convention!
        self.pedigree.doubledHaploid.extend([False] * nFounder)
        self.pedigree.founderInbreeding.update(founderInbreeding)
        self.pedigree.sort()

#        ctFGam = 0
        for id_ in ID:
            gametes = [None, None]
            for i in range(2):
                gametes[i] = FounderGamete(
                                    chromLengths=self.genome.chromLengths,
                                    lineage=getNewID(),
                                    ID=getNewID())
#                ctFGam += 1
            founder = Individual(ID=id_, gametes=gametes)
            self._addIndividual(founder)

    def __update(self):
        """ Update the population according to the pedigree """
        # Determine which individuals exist and which are missing
        IDpop = self.getID()
        IDped = self.pedigree.ID
        ID = list(set(IDped) - set(IDpop))

        for id_ in ID:
            details = self.pedigree.getDetails(id_)
            father = details['father']
            mother = details['mother']
#            generation = details['generation']

            if not father or not mother:
                raise ValueError("Individual has no father or mother.")

            fatherGamete = simulateMeiosis(self.getInd(father),
                                           self.crossoverSimulator)
            motherGamete = simulateMeiosis(self.getInd(mother),
                                           self.crossoverSimulator)

            newInd = Individual(ID=id_,
                                father=self.getInd(father),
                                mother=self.getInd(mother),
                                gametes=[fatherGamete, motherGamete])

            self._addIndividual(newInd)

    def randomMating(self, size, allowSelfing=True):
        """Perform random mating

        Parameters
        ----------
        size : integer
          Number of progeny to be produced.
        allowSelfing : boolean
          Should self-fertilization be allowed?

        """
        self.pedigree.randomMating(size=size, allowSelfing=allowSelfing)
        self.__update()

    def selfing(self, size):
        """Perform self-fertilization

        Parameters
        ----------
        size : integer or list
          A single integer or a list of integers with the same length as the
          current number of individuals

        If size is a single integer, a number of 'size' selfing progeny are
        produced from each individual. If it is a list, a variable number of
        progeny is produced as specified in the list.

        """
        self.pedigree.selfing(size)
        self.__update()

    def roundRobin(self):
        self.pedigree.roundRobin()
        self.__update()

    def getID(self):
        """Get the IDs of all individuals in the population."""
        return list(self.individuals.keys())

    def getInd(self, ID):
        """Get an Individual.

        Parameters
        ----------
        ID : string/integer
          Identifier

        Return
        ------
        out : Individual

        """
        return self.individuals[ID]

    # TODO: check the argument. Maybe change default and allow many generations
    def realIBD(self, generation=None):
        """Compute the realized IBD coefficient between all individuals in a
        generation

        Parameters
        ----------
        generation : integer (defaults to the last generation)
          The generation for which realized IBD coefficients should be
          calculated.

        """
        if generation is None:
            generation = self.pedigree.lastGen()
        ID = self.pedigree.getID(generation)
        nID = len(ID)
        IBD = np.empty(shape=(nID, nID), dtype=np.float_)

        for i in range(nID):
            for j in range(i):
                IBD[i, j] = self.pairRealIBD(ID[i], ID[j])
                IBD[j, i] = IBD[i, j]
            IBD[i, i] = self.pairRealIBD(ID[i])

        return IBD

    def pairRealIBD(self, ID1, ID2=None):
        """Compute the realized IBD coefficient between two individuals

        Parameters
        ----------
        ID1 : string/integer
          Identifier of the first individual
        ID2 : string/integer (default = None)
          Identifier of the second individual

        Return
        ------
        out : float
          The realized IBD coefficient between the given individiuals.

        """
        ind1 = self.getInd(ID1)
        if ID2:
            ind2 = self.getInd(ID2)
        else:
            ind2 = self.getInd(ID1)

        IBD = 0.25 * (self.__realIBDGam(ind1.gametes[0], ind2.gametes[0]) +
                      self.__realIBDGam(ind1.gametes[0], ind2.gametes[1]) +
                      self.__realIBDGam(ind1.gametes[1], ind2.gametes[0]) +
                      self.__realIBDGam(ind1.gametes[1], ind2.gametes[1]))
        return IBD

    def __realIBDGam(self, gam1, gam2):
        """Compute the realized IBD coefficient between two gametes."""
        gam1 = gam1.skeleton
        gam2 = gam2.skeleton
        chromLengths = self.genome.chromLengths
        numChrom = self.genome.numChrom
#        chromLengths = [chrom['locations'][-1] for chrom in gam1]  # dreadful

        cumlen = 0.0
        for i, chrom1, chrom2 in zip(range(numChrom), gam1, gam2):
            loc1, loc2 = chrom1['locations'], chrom2['locations']
            lin1, lin2 = chrom1['lineages'], chrom2['lineages']
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

        return cumlen / sum(chromLengths)

    def __addFleshGamete(self, gamete):
        """Add Genotypes to a single gamete"""
        numChrom = self.genome.numChrom
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

    def addFlesh(self, generation=None):
        """Add genotypes to individuals

        Parameters
        ----------
        generation : integer (default = all generations)
          The generation for which genotypes should be computed.

        """
        if generation is None:
            generation = list(set(self.pedigree.generation))
        else:
            generation = [generation]

        for gen in generation:
            ID = self.pedigree.getID(gen)
            for id_ in ID:
                ind = self.getInd(id_)
                self.__addFleshGamete(ind.gametes[0])
                self.__addFleshGamete(ind.gametes[1])

    def genoMatrixDict(self, generation=None):
        """Compute genotype matrices for individuals

        Parameters
        ----------
        generation : integer (default = last generation)
          The generation for which genotype matrices should be computed.

        Return
        ------
        out : dictionary
          A dictionary containing the genotypes.

        """
        if generation is None:
            generation = self.pedigree.lastGen()

        ID = self.pedigree.getID(generation)
        dictionary = dict()
        for id_ in ID:
            ind = self.getInd(id_)
            dictionary[id_] = np.vstack((np.hstack(ind.gametes[0].flesh),
                                         np.hstack(ind.gametes[1].flesh)))
        return dictionary

    def genoMatrix(self, generation=None):
        """Compute a genotype matrix

        Parameters
        ----------
        generation : integer (default = last generation)
          The generation for which the genotype matrix should be computed.

        Return
        ------
        out : dictionary
          A dictionary holding the the IDs of the individuals, locus names if
          provided and the genotype matrix.

        """
        pass

    def constructPedigree(self):
        """ Construct a new Pedigree from the individuals in the population.

        Return
        ------
        out : Pedigree
          A new pedigree.

        """
        Ind = self.individuals
        ID = list(Ind.keys())
        father = [None] * len(ID)
        mother = [None] * len(ID)
        for ix, id1 in enumerate(ID):
            ind1 = Ind[id1]
            for id2 in ID:
                ind2 = Ind[id2]
                if ind1.father is ind2:
                    father[ix] = ind2.name
                if ind1.mother is ind2:
                    mother[ix] = ind2.name

        pedigree = Pedigree()
        pedigree.id = ID
        pedigree.father = father
        pedigree.mother = mother

        pedigree.sort()
        return pedigree

    def realFounderGeneContr(self, ID):
        """Calculate realized gene contribution of founders

        Parameters
        ----------
        ID : string/integer or list (default = all individuals)

        Return
        ------
        out : dictionary
          A dictionary containing a matrix with founder gene contributions,
          a list with column names (founders) and row names (individuals).

        """
        if ID is None:
            ID = self.pedigree.ID()
        elif isinstance(ID, (str, int)):
            ID = [ID]

        lookup = {}
        for i, ind in self.founderIndividuals.items():
            for gamete in ind.gametes:
                lookup[gamete.lineage] = i
#                try:
#                    lookup[gamete.lineage].append(i)
#                except KeyError:
#                    lookup[gamete.lineage] = [i]

        founderID = list(self.founderIndividuals.keys())
        mat = np.zeros(shape=(len(ID), len(self.founderIndividuals)),
                       dtype=np.float_)

        for i, id_ in enumerate(ID):
            contr = self.__singleRealFounderGeneContribution(self, id_)
            for lineage, value in contr.items():
                mat[i, founderID.index(lookup[lineage])] += value

        return {'colnames': founderID, 'rownames': ID, 'matrix': mat}

    def __singleRealFounderGeneContr(self, ID):
        """Calculate realized gene contribution of founders to a single
        individual."""
        ind = self.getInd(ID)
        d = {}
        totalLen = 2.0 * sum(self.genome.chromLengths)
        for gamete in ind.gametes:
            for ix, chrom in enumerate(gamete.skeleton):
                lineages = chrom['lineages']
                locations = [loc/totalLen for loc in chrom['locations']]
                left = 0.0
                for lin, loc in zip(lineages, locations):
                    try:
                        d[lin] += loc - left
                    except KeyError:  # lineage is not yet in the dictionary
                        d[lin] = loc - left
                    left = loc  # shift forwards

        if not math.isclose(sum(d.values()), 1.0):
            raise RuntimeError("Chunk lengths do not add up! BUG")

        return d

    def realGeneContrMatrix(self, IDa, IDp):
        # This implementation is horribly inefficient
        # TODO: Group IDa into generations, then compute
        """Calculate a matrix of realized gene contributions.

        Parameters
        ----------
        IDa : list
          List of parent individuals
        IDp : list
          List of progeny individuals

        Return
        ------
        out : dictionary
          A dictionary containing a matrix of realized gene contributions,
          a list with row names (IDp) and a list with column names (IDa).

        """
        mat = np.matlib.zeros(shape=(len(IDp), len(IDa)), dtype=np.float_)
        for ixa, ida in enumerate(IDa):
            for ixp, idp in enumerate(IDp):
                mat[ixp, ixa] = self.pairRealGeneContr(IDa=ida, IDp=idp)

        return {'matrix': mat, 'rownames': IDp, 'colnames': IDa}

    def pairRealGeneContr(self, IDa, IDp):
        """Calculate realized gene contribution of one individual to another.

        Parameters
        ----------
        IDa : string/integer
          Identifier of parent individual
        IDp : string/integer
          Identifier of progeny individual

        Return
        ------
        out : float
          Realized gene contribution of IDa to IDp.

        """
        ga = self.pedigree.getGeneration(IDa)
        gp = self.pedigree.getGeneration(IDp)
        if not ga < gp:
            raise ValueError("Individual %s is younger than %s" % (IDa, IDp))

        if not self.pedigree.isAncestor(IDa, IDp):
            return 0.0

        contr = self.realGeneContr(ID=IDp, maxGenBack=gp - ga)
        try:
            return(contr[IDa])
        except KeyError:
            return 0.0

    def realGeneContr(self, ID, maxGenBack):
        """Calculate the realized gene contribution of ancestors of a single
        individual.

        Parameter
        ---------
        ID : string/integer
          Identifier of the individual
        maxGenBack : integer
          The maximum number of generations to look back. This is internally
          restricted to the generation of the individual itself, because
          looking farther back than to the founders is not possible.

        Return
        ------
        out : dictionary
          A dictionary containing the gene contributions of the respective
          ancestors.

        """
        class Chunk:
            def __init__(self, individual, gamete, chrom, parent, start, stop):
                self.individual = individual  # individual holding the chunk
                self.gamete = gamete  # gamete the chunk refers to
                self.chrom = chrom  # chromosome index
                self.parent = parent  # individual that produced the chunk
                self.start = start
                self.stop = stop

        ind = self.getInd(ID)
        maxGenBack = min(maxGenBack, self.pedigree.getGeneration(ID))
        chromLengths = self.genome.chromLengths
        allChunks = set()

        for chrom, chromLen in enumerate(chromLengths):
            chunkSet = set([Chunk(ind, ind.gametes[0], chrom, ind.father,
                                  0.0, chromLen),
                            Chunk(ind, ind.gametes[1], chrom, ind.mother,
                                  0.0, chromLen)])

            for counter in range(maxGenBack):
                newChunkSet = set()
#                IS = iter(chunkSet)
#                chunk = next(IS)
                for chunk in chunkSet:
                    gamete = chunk.gamete
#                    individual = chunk.individual
                    parent = chunk.parent
                    start = chunk.start
                    stop = chunk.stop
                    if parent is None:  # reached end of pedigree, place back
                        newChunkSet.add(chunk)
                        continue

                    tmp = gamete.skeleton[chunk.chrom]  # info stored here
                    Xlocations = tmp['Xlocations'].copy()
                    grandparents = (parent.father, parent.mother)
                    parentAtStart = tmp['parentStart']

                    parentalGametes = (gamete.father, gamete.mother)
                    cur_par = parentalGametes.index(parentAtStart)
#                    IX = iter(zip([0.0] + Xlocations,
#                                  Xlocations + [chromLen]))
#                    xol, xor = next(IX)
                    for xol, xor in zip([0.0] + Xlocations,
                                        Xlocations + [chromLen]):
                        # distinguish all 6 possible mutually exlusive cases.
                        if xor <= start:
                            continue

                        elif xol < start and start < xor <= stop:
                            nstart = start
                            nstop = xor
#                            print("xol < start and start < xor <= stop")
                        elif xol < start and stop < xor:
                            nstart = start
                            nstop = stop
#                            print("xol < start and stop < xor")
                        elif start <= xol < stop and start < xor <= stop:
                            nstart = xol
                            nstop = xor
#                            print("start <= xol < xor <= stop")
                        elif start <= xol < stop and stop < xor:
                            nstart = xol
                            nstop = stop
#                            print("start <= xol < stop and stop < xor")
                        elif stop < xol:
                            break

#                        print("start: %.2f \t stop: %.2f \t xol: %.2f \t "
#                              "xor: %.2f" % (start, stop, xol, xor))
#                        print("%.5s" % gamete.ID + "\tLevel " +
#                              str(counter + 1) + "\tnew Chunk " +
#                              "\tnstart:%.2f " % (nstart),
#                              "\tnstop:%.2f" % (nstop))
#                        print("\n")
                        newChunkSet.add(Chunk(gamete=parentalGametes[cur_par],
                                              individual=parent,
                                              chrom=chrom,
                                              parent=grandparents[cur_par],
                                              start=nstart,
                                              stop=nstop))

                        cur_par = 1 - cur_par

                chunkSet = newChunkSet

            allChunks.update(chunkSet)
        # analyze allChunks for contributions
        contr = dict()
#        IT = iter(allChunks)
#        chunk = next(IT)
        totalLength = 2.0 * sum(chromLengths)
        for chunk in allChunks:
            start, stop = chunk.start, chunk.stop
            try:
                contr[chunk.individual.ID] += (stop - start) / totalLength
            except KeyError:
                contr[chunk.individual.ID] = (stop - start) / totalLength

        return contr
# %%


class Pedigree(Sequence):

    """Pedigree

    Instance variables
    -----------------
    ID : Indentfiers of individuals

    father : Identifiers of the the male parents

    mother : Identifiers of the female parents

    generation : Generation of the individuals

    Usually, pedigrees do not need to be created by the user, but it is
    possible.


    """
    _vAttr = ['ID', 'father', 'mother', 'generation', 'doubledHaploid']

    def __init__(self, ID, father, mother,
                 generation=None,
                 doubledHaploid=None,
                 founderInbreeding=None):
        """Instantiate a new Pedigree

        Parameters
        ----------
        ID : list
          Identfiers of individuals
        father : list
          Identifiers of the the male parents
        mother : list
          Identifiers of the female parents
        generation : list
          Generation of the individuals

        """

        self.ID = ID
        self.checkIDNone()
        self.checkIDUnique()

        self.mother = mother
        self.father = father
        self.checkIDLength()

        if generation is not None:
            self.generation = generation
            self.checkGenLength()

        if doubledHaploid is None:
            self.doubledHaploid = [False] * len(ID)
        else:
            self.doubledHaploid = doubledHaploid
            self.checkDHSelfing()

        if founderInbreeding is None:
            founderID = [i for i in self.ID if self.isFounder(i)]
            self.founderInbreeding = dict.fromkeys(founderID, 0.0)
        else:
            self.founderInbreeding = founderInbreeding
            self.checkFounderInbreeding()
            # TODO: Special check that inbreeding coefficient for doubled
            # haploid founders are 2.0.

    def checkIDNone(self):
        """Check if any ID is None."""
        if any(i is None for i in self.ID):
#            return False
            raise ValueError("IDs cannot be 'None'.")
#        return True

    def checkIDUnique(self):
        """Check if all IDs are unique."""
        if not allUnique(self.ID):
#            return False
            raise ValueError("Not all IDs are unique")
#        return True

    def checkIDLength(self):
        """Check if ID, father and mother have the same length."""
        if not (len(self.ID) == len(self.father) == len(self.mother)):
#            return False
            raise ValueError("'ID', 'father' and 'mother' must have "
                             "same length.")
#        return True

    def checkGenLength(self):
        """Check if generation has the correct length."""
        if not len(self.generation) == len(self.ID):
            raise ValueError("'generation' must have appropriate length.")

    def checkOwnAncestor(self):
        """Check if all Individuals are not their own ancestor."""
        if any((i == f or i == m for i, f, m in
                zip(self.ID, self.father, self.mother))):
#            return False
            raise ValueError("An individual has itself as ancestor.")
#        return True

    def checkDHSelfing(self):
        """Check if all DH are selfing progeny."""
        if any(not f == m for d, f, m in
               zip(self.doubledHaploid, self.father, self.mother) if d):
#            return False
            raise ValueError("Not all doubled haploids are "
                             "selfing progeny.")
#        return True

    def checkFounderInbreeding(self):
        """Check if inbreeding coefficients of all founders are provided."""
        founderID = [i for i in self.ID if self.isFounder(i)]
        diff = set(founderID) - self.founderInbreeding.keys()
        if diff:
#            return False
            raise ValueError(
                    "If 'founderInbreeding' is provided, inbreeding "
                    "coefficients for all founder individuals must "
                    "be specified. They are missing for the following "
                    "individuals:\n\n " + '{}'.format(diff))
#        return True
#    def check(self):
#        """ Check a pedigree.
#
#        Performs a series of checks on a pedigree.
#
#        """
#        checkList = ['_isIDNotNone', '_isIDUnique',
#                     '_isIDFatherMotherSameLength', '_isNotOwnAncestor',
#                     '_isDHSelfingProgeny', '_isFounderInbreedingProvided']
#
#        {getattr(self, cn)() for cn in checkList}
#
#        # Check if pedigree is sorted
#        if not self.isSorted():
#            raise ValueError("The pedigree is not sorted.")
#
#        hasattr(self, 'doubledHaploid')
#
#
#        hasattr(self, 'founderInbreeding')

    def __getitem__(self, item):
        return self.getDetails(item)

    def __len__(self):
        return self.ID.__len__()

    def isExtended(self):
        """Check if the pedigree is extended.

        Return
        ------
        out : boolean
          True, if the pedigree is extended, otherwise, False.

        A pedigree is extended if all parents are themselfs in the pedigree
        with unknown origin.

        """
        # Note that bool(0) is False
        for temp in zip(self.father, self.mother):
            for p in temp:
                if not (p is None or p in self.ID):
                    return False
        return True

    def extend(self):
        """ Extend the pedigree.

        Inbreeding of founders is assumed to be 1.0.

        """
        for temp in zip(self.father, self.mother):
            for p in temp:
                if not (p is None or p in self.ID):
                    # append to start
                    self.ID.insert(0, p)
                    self.father.insert(0, None)
                    self.mother.insert(0, None)
                    if hasattr(self, 'generation'):
                        self.generation.insert(0, 0)
                    if hasattr(self, 'doubledHaploid'):
                        self.doubledHaploid.insert(0, False)
                    self.founderInbreeding[p] = 1.0

    def checkHasGeneration(self):
        if not hasattr(self, 'generation'):
            raise AttributeError("Data on generations is needed.")

    def checkHasDoubledHaploid(self):
        if not hasattr(self, 'doubledHaploid'):
            raise AttributeError("Data on doubled haploidy is needed.")

    def checkIsExtended(self):
        if not self.isExtended():
            raise RuntimeError("Pedigree is not extended.")

    def numInd(self, generation):
        """Get the number of all individuals in a given generation

        Parameters
        ----------
        generation : integer
          The generation

        Return
        ------
        out : integer
          Number of individuals in generation 'generation'

        """
        self.checkHasGeneration()
        return self.generation.count(generation)

    def lastGen(self):
        """Get the number of the last generation

        Return
        ------
        out : integer
          The number of the last generation.

        Counting starts at 0 for the generation of the founders.

        """
        self.checkHasGeneration()
        return max(self.generation)

    def numGen(self):
        """Get the number of generations in the pedigree."""
        self.checkHasGeneration()
        return len(set(self.generation))

    def getID(self, generation):
        """Get IDs of individuals in a given generation

        Parameters
        ----------
        generation : integer
          The generation for which IDs should be returned.

        Return
        ------
        out : list
          A list with IDs.

        """
        self.checkHasGeneration()
        return [i for i, g in zip(self.ID, self.generation) if g == generation]

    def getDoubledHaploids(self, generation):
        """Get IDs of doubled haploids in a given generation

        Parameters
        ----------
        generation : integer
          The generation for which IDs should be returned.

        Return
        ------
        out : list
          A list with IDs.

        """
        self.checkHasGeneration()
        return [i for i, g, d in zip(self.ID, self.generation,
                self.doubledHaploid) if g == generation and d]

    def getIndex(self, ID):
        """Get index of an individual in the pedgiree

        Parameters
        ----------
        ID : integer or string
          The identifier of the individual

        Return
        ------
        out : integer
          Index in the pedigree

        """
        return self.ID.index(ID)

    def getFather(self, ID):
        """Get the ID of the father of an individual

        Parameters
        ----------
        ID : integer or string
          The identifier of the individual

        Return
        ------
        out : integer or string
          The identifier of the father

        """
        return self.father[self.getIndex(ID)]

    def getMother(self, ID):
        """Get the ID of the mother of an individual

        Parameters
        ----------
        ID : integer or string
          The identifier of the individual

        Return
        ------
        out : integer or string
          The identifier of the mother

        """
        return self.mother[self.getIndex(ID)]

    def getGeneration(self, ID):
        """Get the genration of an individual

        Parameters
        ----------
        ID : integer or string
          The identifier of the individual

        Return
        ------
        out : integer
          The generation of the individual

        """
        self.checkHasGeneration()
        return self.generation[self.getIndex(ID)]

    def isDoubledHaploid(self, ID):
        """Check if individual is a doubled haploid

        Parameters
        ----------
        ID : integer or string
          The identifier of the individual

        Return
        ------
        out : boolean
          True if the individual is doubled haploid, otherwise False.

        """
        self.checkHasDoubledHaploid()
        return self.doubledHaploid[self.getIndex(ID)]

    def getDetails(self, ID):
        """Obtain details for a pedigreed individual.

        Parameters
        ----------
        ID : string/integer
          Identifier of the individual

        Return
        ------
        out : dictionary
          A dictionary containing the father, mother and generation number of
          the individual.

        """
        index = self.getIndex(ID)
#        keys = set(self.__dict__.keys())
#        keys.discard('founderInbreeding')
        out = collections.OrderedDict((name, getattr(self, name)[index]) for name in Pedigree._vAttr)
        if ID in self.founderInbreeding:
            out['founderInbreeding'] = self.founderInbreeding[ID]
        return out

    # TODO: Thoroughly test this method
    def isSorted(self):
        """Check if the pedigree is sorted."""
        seen = set([None])
        return not any(f not in seen or m not in seen or seen.add(i)
                       for i, f, m in zip(self.ID, self.father, self.mother))

#        seen = set([None])
#        for i, f, m in zip(self.ID, self.father, self.mother):
#            print(i, f, m)
#            if f not in seen or m not in seen:
#                print("not in seen")
#                break
#            seen.add(i)

    def sort(self, placeFounderFirst=True):
        """Sort the pedigree.

        If no data on generations are available, the algorithm provided in
        `Zhang et al. (2009) <http://www.medwelljournals.com/fulltext/
        ?doi=javaa.2009.177.182>`_ is applied to sort from scratch.

        """
        # If data on generations is availabe, the task is easy
        if hasattr(self, 'generation'):
            index = order(self.generation)
            for name in Pedigree._vAttr:
                if hasattr(self, name):
                    setattr(self, name, reorder(getattr(self, name), index))

        else:
            n = len(self.ID)  # We are sure that all IDs are unique.
            stru = dict(zip(self.ID,
                        map(list,
                            zip(self.father, self.mother, [0]*n))))
            ID, mother, father, generation = [], [], [], []
            while stru:
                notparents = []
                for id1 in list(stru.keys()):  # list(), because it changes.
                    for id2, (f, m, _) in stru.items():
                        if not id1 == id2:  # Skip this comparison.
                            # Check if i is a parent and break if necessary.
                            if id1 == f or id1 == m:
                                break
                    # If the inner for loop ends normally, i.e., if i is not
                    # a parent, else is conducted.
                    else:
                        notparents.append(id1)
                        ID.append(id1)
                        temp = stru[id1]
                        father.append(temp[0])
                        mother.append(temp[1])
                        generation.append(temp[2])

                # Remove nonparents and add +1 to generation
                for id_ in list(stru.keys()):
                    if id_ in notparents:
                        del stru[id_]
                    else:
                        stru[id_][2] += 1

            maxgen = max(generation)
            generation = [abs(gen - maxgen) for gen in generation]
            self.ID, self.father = ID, father
            self.mother, self.generation = mother, generation

            if placeFounderFirst:
                for id_ in self.ID:
                    ix = self.getIndex(id_)
                    if self.isFounder(id_):
                        self.generation[ix] = 0

            # The pedigree has no generation information. A recursive call
            # sort() will use this information for final sorting!
            # TODO: Test the functionality here.
            self.sort()

    def randomMating(self, size, allowSelfing=True):
        """Perform random mating."""
        curgen = self.lastGen()
        curID = self.getID(curgen).copy()
        if not allowSelfing and len(curID) < 2:
            raise ValueError(
                "It is not possible to have 'allowSelfing = False' if there "
                "are less than 2 individuals in the population.")

        self.ID.extend(getNewID() for _ in range(size))
        for i in range(size):
            id1 = random.choice(curID)
            id2 = random.choice(curID)
            if not allowSelfing:
                while id1 == id2:
                    id2 = random.choice(curID)
            self.mother.append(id1)
            self.father.append(id2)

        self.generation.extend([curgen + 1] * size)
        self.doubledHaploid.extend([False] * size)

    def selfing(self, size=1):
        """Perform self-fertilization."""
        curgen = self.lastGen()
        curID = self.getID(curgen).copy()

        if size == 1:
            size = [1] * len(curID)
        elif not len(size) == len(curID):
            raise ValueError("'size' must be either 1 or a list with "
                             "length equal to the number of individuals in "
                             "the current generation.")

        parentID = [i for i, n in zip(curID, size) for _ in range(n)]
        self.ID.extend(getNewID() for _ in range(sum(size)))
        self.mother.extend(parentID)
        self.father.extend(parentID)
        self.generation.extend([curgen + 1] * sum(size))
        self.doubledHaploid.extend([False] * sum(size))

    # TODO: Add 'size' like for selfing, think about integration with
    # family memberships.
    def roundRobin(self):
        """Perform mating according to a round robin design."""
        curgen = self.lastGen()
        curID = self.getID(curgen).copy()
        if len(curID) < 3:
            raise ValueError("RoundRobin with less than 3 individuals "
                             "is not possible.")
        self.ID.extend(getNewID() for _ in range(len(curID)))
        for id1, id2 in zip(curID, curID[1:] + curID[:1]):
            self.mother.append(id1)
            self.father.append(id2)
        self.generation.extend([curgen + 1] * len(curID))
        self.doubledHaploid.extend([False] * len(curID))
    
    
    def synthetic(self, size=1, reciprocals=False):
        """Perform mating according to a synthetic design."""
        curgen = self.lastGen()
        curID = self.getID(curgen).copy()
        n = len(curid)  # number of current individuals
        if n < 2:
            raise ValueError("A Synthetic with less than 2 individuals "
                             "is not possible.")
        nmat = n * (n - 1) / (2 - int(reciprocals))  # number of matings
        if isinstance(size, int):
            size = [1] * nmat
        elif not len(size) == nmat:
            raise ValueError("'size' must be either an integer or a list with "
                             "length equal to the number of mating in.")
        self.ID.extend(getNewID() for _ in range(sum(size)))
        # TODO: This solution looks kinda crappy and must be revised.
        ct = 0
        for ix1, id1 in enumerate(curID):
            for ix2, id2 in enumerate(curID[ix1 + 1:], start = ix1 + 1):
                self.mother.extend([id1] * size[ct])
                ct += 1
                self.father.extend([id2] * size[ct])
                ct += 1
                if reciprocals:
                    self.mother.extend([id2] * size[ct])
                    ct += 1
                    self.father.extend([id1] * size[ct])
                    ct += 1
        self.generation.extend([curgen + 1] * nmat)
        self.doubledHaploid.extend([False] * nmat)         
    
    
    

    def countOffspring(self, ID):
        """ Count the number of offs

        Parameters
        ----------
        ID : string/integer or list (default = None)
          Can be a single identifier or a list of identifiers.

        Return
        ------
        out : integer/list
          If ID is a single identifier, an integer. Otherwise, a list.
          For the default (ID = None), all individuals in the pedigree are
          considered (according to their order).

        """
        if isinstance(ID, (str, int)):
            ID = [ID]
        if not set(self.ID).issuperset(set(ID)):
            raise ValueError("Not all individuals in ID are contained in "
                             "the pedigree")

        noff = [0] * len(ID)
        for ix, i in enumerate(ID):
            for f, m in zip(self.father, self.mother):
                noff[ix] += (f == i or m == i)

        if len(ID) == 1:
            return noff[0]
        return noff

    def isFounder(self, ID):
        """Test if individuals are founders

        Parameters
        ----------
        ID : string/integer or list (default = None)
          Can be a single identifier or a list of identifiers.

        Return
        ------
        out : integer/list
          If ID is a single identifier, an integer. Otherwise, a list.
          For the default (ID = None), all individuals in the pedigree are
          considered (according to their order).


        Founder individuals have both parents unkown

        """
        allID = self.ID
        if isinstance(ID, (str, int)):
            ID = [ID]

        isfounder = [None] * len(ID)
        for ix, i in enumerate(ID):
            try:
                index = allID.index(i)
            except ValueError:
                raise ValueError("The individual %s is not "
                                 "present in the pedigree." % (str(i)))
            else:
                father = self.father[index]
                mother = self.mother[index]
                isfounder[ix] = not father and not mother

        if len(ID) == 1:
            return isfounder[0]
        return isfounder


    def IDGenerationHandler(self, ID, generation):
        if ID is None:
            if generation is None:
                raise ValueError("Either 'ID' or 'generation' must be given.")
            ID = self.getID(generation)
        else:
            if generation is not None:
                warnings.warn("'ID' was given, so 'Generation' is ignored")
            if isinstance(ID, (str, int)):
                ID = [ID]
        return ID

    # TODO: revise, add option for having A already. But beware: parents must
    # always be present!
    # TODO: Add check for founder individuals, because here, it is impossible
    # to compute D.
    # TODO: check for basic funcitonality
    def makeD(self, ID=None, generation=None, method='tabular'):
        ID = self.IDGenerationHandler(ID, generation)
        A = self.makeA(self.ID, generation = None, method = method)
        dmat = numpy.matlib.identity(n=len(ID), dtype=np.float_)
        for indexi, i in enumerate(ID):
            ixi = self.ID.index(i)
            ixfi = self.ID.index(self.father[ixi])
            ixmi = self.ID.index(self.mother[ixi])
            for indexj, j in enumerate(ID[indexi + 1:], start = indexi + 1):
                ixj = self.ID.index(j)
                ixfj = self.ID.index(self.father[ixj])
                ixmj = self.ID.index(self.mother[ixj])
                dom = 0.25 * (A[ixfi, ixfj] + A[ixmi, ixmj] +
                              A[ixfi, ixmj] + A[ixmi, ixfj])
                dmat[indexi, indexj] = dom
                dmat[indexj, indexi] = dom

            dmat[indexi, indexi] = 0.25 * (A[ixfi, ixfi] + A[ixmi, ixmi] +
                                           2.0 * A[ixfi, ixmi])
        return dmat

    def makeAinv(self, ID=None, generation=None, method='tabular'):
        if not self.isSorted():
            raise RuntimeError("The pedigree must be sorted.")
        ID = self.IDGenerationHandler(ID, generation)
        if method == 'tabular':
            A = self._tabularMethod(ID)
            return np.linalg.inv(A)
            # Cholesky approach not faster even for large (5000) matrices.
#            n = A.shape[0]
#            L = np.linalg.cholesky(A)
#            Y = scipy.linalg.solve_triangular(L, np.eye(n), lower=True, check_finite=False)
#            X = scipy.linalg.solve_triangular(L, Y, trans=1, lower=True, overwrite_b=True, check_finite=False)
        elif method == 'LDL':
            LDL = self._LDLMethod(ID)
            Linv = np.linalg.inv(LDL['L'])
            return Linv.T * np.diag(1.0 / LDL['D']) * Linv
        else:
            raise ValueError("'method' must be either 'tabular' "
                             "or 'LDL'.")

    def makeA(self, ID=None, generation=None, method='tabular'):
        if not self.isSorted():
            raise RuntimeError("The pedigree must be sorted.")
        ID = self.IDGenerationHandler(ID, generation)
        if method == 'tabular':
            return self._tabularMethod(ID)
        elif method == 'LDL':
            LDL = self._LDLMethod(ID)
            return LDL['L'] * np.diag(LDL['D']) * LDL['L'].T
        else:
            raise ValueError("'method' must be either 'tabular' "
                             "or 'LDL'.")

    def inbreeding(self, ID=None, generation=None, method='tabular'):
        if not self.isSorted():
            raise RuntimeError("The pedigree must be sorted.")

        ID = self.IDGenerationHandler(ID, generation)
        if method == 'tabular':
            A = self._tabularMethod(ID)
            return (np.diag(A) - 1.0).tolist()
        elif method == 'LDL':
            return self._LDLMethod(ID)['F'].tolist()
        else:
            raise ValueError("'method' must be either 'tabular' "
                             "or 'LDL'.")
        return A

    # TODO: This is considerable slower than the tabular method. Profile for
    # finding the reason!
    def _LDLMethod(self, ID):
        n = len(self.ID)
        L = numpy.matlib.identity(n=n, dtype=np.float_)
        indices = [None] * len(ID)
        for ix, i in enumerate(ID):
            try:
                indices[ix] = self.ID.index(i)
            except ValueError:
                raise ValueError("The individual %s is not present."
                                 % (str(i)))
        for j in range(n):
            for i in range(j + 1, n):
                fi = self.father[i]
                mi = self.mother[i]
                if fi is not None:
                    if mi is not None:
                        ixfi = self.ID.index(fi)
                        ixmi = self.ID.index(mi)
                        L[i, j] = 0.5 * (L[ixfi, j] + L[ixmi, j])
                    else:
                        ixfi = self.ID.index(fi)
                        L[i, j] = 0.5 * L[ixfi, j]
                else:
                    if mi is not None:
                        ixmi = self.ID.index(mi)
                        L[i, j] = 0.5 * L[ixmi, j]
                    else:
                        L[i, j] = 0

        # Relies on the first individual beeing a founder.
        D = np.zeros(n)
        F = np.zeros(n)
        for j in range(n):
            fj = self.father[j]
            mj = self.mother[j]
            if fj is not None:
                if mj is not None:
                    if self.doubledHaploid[j]:
                        D[j] = 1.0  # set directly to zero
                    else:
                        ixfj = self.ID.index(fj)
                        ixmj = self.ID.index(mj)
                        D[j] = 0.5 - 0.25 * (F[ixfj] + F[ixmj])
                else:
                    ixfj = self.ID.index(fj)
                    D[j] = 0.75 - 0.25 * F[ixfj]
            else:
                if mj is not None:
                    ixmj = self.ID.index(mj)
                    D[j] = 0.75 - 0.25 * F[ixmj]
                else:  # The individual must be a founder!
                    id_ = self.ID[j]
                    ## Check this thourougly!
                    D[j] = 1.0 + self.founderInbreeding[id_]
#                    D[j] = 1 # normal code

#            temp = 0.0
#            for i in range(j + 1):
#                temp += L[j, i]**2 * D[i]
#            I[j] = temp
#            id_ = self.ID[j]
#            if id_ in self.founderInbreeding:
#                F[j] = self.founderInbreeding[id_]
#            else:
            F[j] = np.square(L[j,]) * D.reshape(n,1) - 1.0

        L = L[np.ix_(indices, indices)]
        D = D[np.ix_(indices)]
        F = F[np.ix_(indices)]
        return {'L': L, 'D': D, 'F': F}

    # TODO: compute IBD matrix only till the youngest specified individual.
    def _tabularMethod(self, ID):
        """ Compute the numberator relationship (IBD-) matrix.

        Parameters
        ----------
        ID : string/integer or list (default = all individuals)
          Identifiers for the individuals.

        Return
        ------
        out : dictionary
          A dictionary with the IBD matrix and the column/row names.

        """
#        if not self.isSorted():
#            raise RuntimeError("The pedigree must be sorted.")
        mat = numpy.matlib.identity(n=len(self.ID), dtype=np.float_)
        indices = [None] * len(ID)
        for ix, i in enumerate(ID):
            try:
                indices[ix] = self.ID.index(i)
            except ValueError:
                raise ValueError("The individual %s is not present."
                                 % (str(i)))
        asc1 = self.father
        asc2 = self.mother
        testAsc1 = list(map(bool, asc1))
        testAsc2 = list(map(bool, asc2))
        # subjects where both parents are known
        testAsc = [x and y for x, y in zip(testAsc1, testAsc2)]
        # alternative: list(map(lambda x,y: x and y, testAsc1, testAsc2))
        # subjects with at least one ascendant known
        set_ = [i for i, t in enumerate(zip(testAsc1, testAsc2)) if any(t)]

        for i in set(range(len(self.ID))) - set(set_):
            id_ = self.ID[i]
            if id_ in self.founderInbreeding:  # This is likely redundant !!!
                    mat[i, i] = 1.0 + self.founderInbreeding[id_]
        # --- Core ---
        for i in set_:
            # Diagonal
            if testAsc[i]:
                # test for doubled haploidy
                # TODO: confirm and test this concept
                if self.doubledHaploid[i]:
                    mat[i, i] = 2.0
                else:
                    mat[i, i] = 1.0 + 0.5 * mat[self.ID.index(asc1[i]),
                                                self.ID.index(asc2[i])]
            # Off-diagonal
            j = range(i)  # working with lower triangle
            if testAsc1[i]:
                tmp1 = 0.5 * mat[self.ID.index(asc1[i]), j]
            else:
                tmp1 = 0.0

            if testAsc2[i]:
                tmp2 = 0.5 * mat[self.ID.index(asc2[i]), j]
            else:
                tmp2 = 0.0

            mat[i, j] = tmp1 + tmp2
            mat[j, i] = mat[i, j]

        return mat[np.ix_(indices, indices)]

    # TODO: test for basic functionality
    # TODO: make more efficient
    # TODO: include handler
    def subset(self, ID=None, generation=None):
        """ Subset a pedigree

        Parameters
        ----------
        ID : string/integer or list
          Identifiers for individuals

        Return
        ------
        out : Pedigree
          A new pedigree instance containing only individuals in ID and their
          ancestors.

        """
        if not self.isSorted():
            raise ValueError("The pedigree must be sorted.")

        if ID is None:
            if generation is None:
                raise ValueError("Either 'ID' or 'generation' must be given.")
            ID = self.getID(generation)
        else:
            if generation is not None:
                warnings.warn("'ID' was given, so 'Generation' is ignored")
            if isinstance(ID, (str, int)):
                ID = [ID]

        father, mother, doubledHaploid, generation = [], [], [], []
        for i in reversed(self.ID):
            if i in ID:
                details = self.getDetails(i)
                f = details['father']
                if f not in ID and f is not None:
                    ID.append(f)
                m = details['mother']
                if m not in ID and m is not None:
                    ID.append(m)

        ID.reverse()
        # This can be solved much more efficiently
        for i in ID:
            index = self.ID.index(i)
            father.append(self.father[index])
            mother.append(self.mother[index])
            doubledHaploid.append(self.doubledHaploid[index])
            generation.append(self.generation[index])

        return Pedigree(ID, father, mother, generation, doubledHaploid)

    # TODO: Check whether (if possible) the matrix L in the decomposion of A
    # already contains this information somehow.
    def pairExpGeneContr(self, IDa, IDp):
        """Calculate the expected gene contribution form one individual to
        another

        Parameters
        ----------
        IDa : string/integer
          Identifier of the ancestral individual
        IDp : string/integer
          Identifier of the progeny individual

        Return
        ------
        out : float
          The expected gene contribution from IDa to IDp

        """
        # TODO: write more efficient version of the test
        self.checkHasGeneration()

        ixa = self.ID.index(IDa)
        ixstart = self.ID.index(IDp)
        ga = self.generation[ixa]
        gp = self.generation[ixstart]
        if not ga < gp:
            raise ValueError("Individual %s is younger than %s" % (IDa, IDp))
        # That's a freaking brain teaser. Bottom to top approach will likely
        # be much more efficient if many comparisions are sought
        # A dictionary counter that tell which individuals we have to look for.
        # It keeps the count for duplicated parents (arising from selfings
        # or common lineages)
        qdict = collections.Counter([IDp])
        contr = 0.0  # contribution of IDa to IDp
        c = 0.5  # weighting factor, halved for each step to the ancestors
        while qdict:  # as long as there is some lineage to query
            nqdict = collections.Counter()  # new counter dict for next round
            ixmax = ixstart
            for i in qdict.keys():
                try:
                    ix = self.ID.index(i, ixa, ixstart + 1)
                except:
                    continue  # if i cannot possibly be a progeny of IDa, skip

                ixmax = min(ixmax, ixstart)  # keep track of the minimum index
                f = self.father[ix]
                m = self.mother[ix]
                for p in (f, m):
                    if p == IDa:  # if parent is IDa, contribute weighted c
                        contr += c*qdict[i]
                    elif p:  # p is not IDa and not None, so step down
                        nqdict.update([p]*qdict[i])  # transmit weight

            ixstart = ixmax  # no need to traverse whole list, start here
            c *= 0.5
            qdict = nqdict

        return contr

    # TODO: Handler for IDs
    def expGeneContrMatrix(self, IDa, IDp):
        """Calculate a matrix of expected gene contributions.

        Parameters
        ----------
        IDa : list
          List of parent individuals
        IDp : list
          List of progeny individuals

        Return
        ------
        out : dictionary
          A dictionary containing a matrix of expected gene contributions,
          a list with row names (IDp) and a list with column names (IDa).

        """
        if isinstance(IDa, (str, int)):
            IDa = [IDa]
        if isinstance(IDp, (str, int)):
            IDp = [IDp]

        mat = np.matlib.zeros(shape=(len(IDp), len(IDa)), dtype=np.float_)
        for ixa, ida in enumerate(IDa):
            for ixp, idp in enumerate(IDp):
                mat[ixp, ixa] = self.pairExpGeneContr(ida, idp)

        return {'matrix': mat, 'rownames': IDp, 'colnames': IDa}

    def isAncestor(self, IDa, IDp):
        """Check if an individual as an ancestor of another one

        Parameters
        ----------
        IDa : string/integer
          Identifier of the suspected ancestor
        IDp : string/integer
          Identifier of the individual

        Return
        ------
        out : boolean
          True, if IDa is an ancestor of IDp, False otherwise.

        """
        if isinstance(IDa, (str, int)):
            IDa = [IDa]
        if isinstance(IDp, (str, int)):
            IDp = [IDp]

        ga = self.getGeneration(IDa)
        compset = set([IDp])
        while compset:
            new_compset = set()
            for cand in compset:
                if cand == IDa:
                    return True
                gcand = self.getGeneration(cand)
                if ga < gcand:
                    dt = self.getDetails(cand)
                    new_compset.update((dt['father'], dt['mother']))
                    new_compset.discard(None)

            compset = new_compset
        return False

    # TODO: Handler for IDs
    def isAncestorMatrix(self, IDa, IDp):
        """Calculate a matrix of ancestry relations

        Parameters
        ----------
        IDa : list
          List of parent individuals
        IDp : list
          List of progeny individuals

        Return
        ------
        out : dictionary
          A dictionary containing a matrix of ancestry relations,
          a list with row names (IDp) and a list with column names (IDa).

        This function returns a matrix indicating which individual in IDa
        are ancestors of which individuals in IDp.

        """
        mat = np.matlib.zeros(shape=(len(IDp), len(IDa)), dtype=np.bool_)
        for ixa, ida in enumerate(IDa):
            for ixp, idp in enumerate(IDp):
                mat[ixp, ixa] = self.isAncestor(IDa=ida, IDp=idp)

        return {'matrix': mat, 'rownames': IDp, 'colnames': IDa}

    def toDataFrame(self):
        """Convert pedigree to a pandas DataFrame

        Return
        ------
        out : DataFrame
          A pandas DataFrame holding the pedigree.

        """
#        if self.isEmpty():
#            return pd.DataFrame()
        temp = collections.OrderedDict((name, getattr(self, name))
                                       for name in Pedigree._vAttr
                                       if hasattr(self, name))
        return pd.DataFrame(temp)

    def recode(self):
        self.checkIsExtended()
        newPed = copy.deepcopy(self)
        look = {k:v for v, k in enumerate(newPed.ID, start = 1)}
        for i in range(len(newPed)):
            newPed.ID[i] = look[newPed.ID[i]]
            try:
                newPed.father[i] = look[newPed.father[i]]
            except KeyError:
                pass
            try:
                newPed.mother[i] = look[newPed.mother[i]]
            except KeyError:
                pass
        newPed.founderInbreeding =\
            {look[k]: v for k, v in newPed.founderInbreeding.items()}
        return newPed

    def __repr__(self):
        return str(self.toDataFrame())

    def write(self, file, sep=' ', na_rep='NA'):
        """Write a pedgiree to a file

        Parameters
        ----------
        file : stringer
          Path to the target file
        sep : single character (default = ' ')
          field separator
        na_rep = string
          How unknown parents should be coded

        This function merely creates a pandas DataFrame (via toDataFrame())
        and writes it to the file (via to_csv()). For more flexibility, do
        these two steps yourself.

        """
        temp = self.pedigree.toDataFrame()
        temp.to_csv(file, sep=sep, na_rep=na_rep, index=False)

        #  TODO: This function is rubbish, revise thoroughly
#    def inferGenerations(self, start=0):
#        """Infer generations in a pedigree
#
#        Parameters
#        ----------
#        start = integer (default = 0)
#          The generation that should be assigned to the founders.
#
#
#        Return
#        ------
#        out : dictionary
#          A dictionary holding the assigned generation for each individual.
#
#        """
#        ID = self.ID
#        n = len(ID)
#        generation = [None] * n
#
#        for i in range(n):
#            founder = self.isFounder(ID[i])
#
#            if founder:
#                generation[i] = start
#            else:
#                genlist = []
#                for parent in [self.father[i], self.mother[i]]:
#                    if not parent:
#                        genlist.append(start)
#                    else:
#                        generation.append(generation[ID.index(parent)])
#
#        return dict(zip(ID, generation))

        # HINT: This function is redundant if empty pedigrees are impossible.
#    def isEmpty(self):
#        """Check if the Pedigree is empty."""
#        if not any(self.__dict__.values()):
#            return True
#        return False

    # HINT: This function is redundant if empty pedigrees are impossible.
#    def _checkEmpty(self):
#        """Throw an error if the pedigree is not empty."""
#        if not self.isEmpty():
#            raise RuntimeError("The pedigree is not empty.")

    # HINT: This method is redundant under altered API
#    def addFounders(self, nFounders, ID=None):
#        """Add founders to an empty pedigree.
#
#        Parameters
#        ----------
#        nFounders : integer
#          The number of founder individuals to be added.
#        ID : list
#          Identifiers of the founder individuals
#
#        If IDs are not provided, UUID are generated.
#
#        """
#        self._checkEmpty()
#        ID = copy.deepcopy(ID)  # Make sure the user won't get bit.
#        if ID is None:
#            ID = [getNewID() for _ in range(nFounders)]
#        elif nFounders != len(ID):
#            raise ValueError("'ID' must have the length 'nFounders'.")
#
#        self.ID = ID
#        self.mother = [None] * nFounders  # references to single None object!
#        self.father = [None] * nFounders
#        self.generation = [0] * nFounders
#        self.doubledHaploid = [False] * nFounders
#        self.nFounders = nFounders

# %%


class CrossoverSimulator:

    """Simulate crossover locations on a single meiotic product using
    the Stahl model.

    Instance Variables
    ------------------
    m : Interference parameter

    p : Proportion of chiasmata from non-interference mechanism

    obligateChiasma : Boolean indicating if obligate chiasma is required

    Chiasma locations are a superposition of two processes: a proportion p
    exhibiting no interference, and a proportion (1-p) following the
    chi-square model with interference parameter m. Crossover locations are
    derived by thinning the chiasma locations with probability 1/2.

    Simulations are under the Stahl model with the interference parameter
    being an integer. This is an extension of the chi-square model, but with
    chiasmata being the superposition of two processes, one following the
    chi-square model and the other exhibiting no interference.

    The source code is a translation from
    `R/simcross <https://github.com/kbroman/simcross/>`_ package

    References
    ----------

    Copenhaver, G. P., Housworth, E. A. and Stahl, F. W. (2002) Crossover
    interference in arabidopsis. Genetics 160, 1631–1639.

    Foss, E., Lande, R., Stahl, F. W. and Steinberg, C. M. (1993) Chiasma
    interference as a function of genetic distance. Genetics 133, 681–691.

    Zhao, H., Speed, T. P. and McPeek, M. S. (1995) Statistical analysis of
    crossover interference using the chi-square model. Genetics 139, 1045–1056.


    """

    def __init__(self, m=10, p=0, obligateChiasma=False):
        """ Instantiate a CrossoverSimulator

        Parameters
        ----------
        m : nonnegative integer (default = 10)
          Interference parameter (m = 0 corresponds to no interference)
        p : float between 0 and 1 (default = 0)
          Proportion of chiasmata from no-interference mechanism
          (p = 0 gives a pure chi-square model)
        obligateChiasma : boolean
          If TRUE, require an obligate chiasma on the 4-strand bundle
          at meiosis.

        """
        self.m = m
        self.p = p
        self.obligateChiasma = obligateChiasma

        @property
        def m(self):
            return self.m

        @property
        def p(self):
            return self.p

        @property
        def obligateChiasma(self):
            return self.obligateChiasma

    def __calculateLStar(self, L):
        m = self.m
        p = self.p

        if L <= 50.0:
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
                denom = 1 - math.exp(-Lstar / 50.0)
            else:
                lambda1 = Lstar/50.0 * (m + 1) * (1 - p)
                lambda2 = Lstar/50.0 * p

                sm = 0.0
                for i in range(m + 1):
                    sm += dpoisson(i, lambda1) * (m + 1 - i) / (m + 1)

                denom = 1 - sm * math.exp(-lambda2)

            return 2 * L - 2 * Lstar / denom

        return scipy.optimize.brentq(f=funcToZero, a=math.exp(-8),
                                     b=L, args=(L, m, p))

    def simulateCrossover(self, L):
        """Simulate crossing over events

        Parameters
        -----------
        L : positive float
          Length of the chromosome in centiMorgan

        Return
        ------
        out : list
          A list with the positions of crossover events

        """
        m = self.m
        p = self.p
        obligateChiasma = self.obligateChiasma
#        print(m,p,obligateChiasma)

        if obligateChiasma:
            Lstar = self.__calculateLStar(L)
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
            # about twice times faster for small nXO (1-5):
            return sorted([L * random.random() for _ in range(nXO)])
#            return sorted(np.random.uniform(0.0, L, nXO).tolist())

        lambda1 = Lstar/50.0 * (m + 1) * (1.0 - p)  # adjust expectation
        lambda2 = Lstar/50.0 * p

        while True:
            # chiasma and intermediate points
            nPoints = np.random.poisson(lambda1)
            # which point is the first chiasma?
            first = random.randrange(m + 1)  # much faster
#            first = np.random.choice(range(m + 1))
            if first > nPoints:
                nIchi = 0
            else:
                nIchi = nPoints // (m + 1) + int(first < (nPoints % (m + 1)))
                # be careful to use explicit integer division here
            # no. chiasma from no interference process
            if p > 0:
                nNIchi = np.random.poisson(lambda2)
            else:
                nNIchi = 0

            if not obligateChiasma or nIchi + nNIchi > 0:
                break

        # locations of chiasmata and intermediate points for process with
        # interference
        pointLocations = sorted([L * random.random() for _ in range(nPoints)])
#        pointLocations = sorted(np.random.uniform(0.0, L, nPoints).tolist())

        # move every (m+1)st point back to front
        nChi = 0
        for j in range(first, nPoints - 1, m + 1):
            # Here, j > nChi is required, otherwise a mess will happen
            pointLocations[nChi] = pointLocations[j]
            nChi += 1

        # chiasma locations from non-interference process
        NIchiLocations = [L * random.random() for _ in range(nNIchi)]
#        NIchiLocations = np.random.uniform(0.0, L, nNIchi).tolist()

        # combine interference and no interference chiasma locations
        chiLocations = sorted(pointLocations[:nChi] + NIchiLocations)

        # thin by 1/2
        nXO = 0
        XOLocations = list()
        for i in range(len(chiLocations)):
            if random.random() < 0.5:  # flip coin -> chiasma
                nXO += 1
                XOLocations.append(chiLocations[i])

        return XOLocations

# %%




#    for i in range(numChrom):
#        numLoci = numLociPerChrom[i]
#        sizeChrom = sizesChrom[i]
#        positions = np.random.uniform(low=0.0, high=sizeChrom, size=numLoci).tolist()
#        lociNames = [getNewID() for _ in range(numLoci)]
#        loci = [Locus(locusName, position) for locusName, position in zip(lociNames, positions)]
#        chrom = Chromosome(chromName = "chrom" + str(i), headPos = 0.0, tailPos = sizeChrom)
#        chrom.addLoci(*loci)
#        genome.addChrom(chrom)
#
#    return genome

# Test Pedigree
#ID = [1, 2, 3, 4, 5]
#father = [None, None, 1, 2, 1, 3]
#mother = [None, None, 2, 1, 3, 3]
#ped = Pedigree(ID=ID, father=father, mother=mother)













#numChrom = 10
#chromLengths = [100] * numChrom
#numLociPerChrom = [1000] * numChrom
#positions = [sorted(np.random.uniform(0.0, h, n)) for h,n in zip(chromLengths, numLociPerChrom)]
#genome = Genome(chromLengths, positions)
#crossoverSimulator = CrossoverSimulator(m = 3, p = 0.5, obligateChiasma = True)
##
##
##
#pop = Population(crossoverSimulator, genome, 'pop')
#
#
#nF = 5
#fID = ['f' + str(i) for i in range(nF)]
#pop.addFounders(nF, fID)
#pop.founderIndividuals
#for _ in range(3):
##    pop.selfing()
#    pop.randomMating(size = 5)

#for i in pop:
#    print(i)


#lastID = pop.pedigree.ID[-1]
#contr = pop.realGeneContr(lastID, 500)
##time.time() - now
#contr
#
#pop.pairRealGeneContr(IDa='f1', IDp=lastID)
#sum(contr.values())
#
#now = time.time()
#tmp=pop.realGeneContrMatrix(IDa=pop.pedigree.getID(0), IDp=pop.pedigree.getID(5))
#print(time.time() - now)
#print(tmp['matrix'])
#
#now = time.time()
#tmp=pop.pedigree.ExpGeneContrMatrix(IDa=pop.pedigree.getID(0), IDp=pop.pedigree.getID(5))
#print(time.time() - now)
#print(tmp['matrix'])
#
#
#ped = pop.pedigree
#
#now = time.time()
#m = ped.isAncestorMatrix(ped.getID(), ped.getID())['matrix']
#time.time() - now
#
#print(m)


#IDp = 314411950031648894434388644473678155678
#IDa = 175397077608902921401047815579261051393
#ped.isAncestor(IDa, IDp)

##print(pop.pedigree)
#ped = pop.pedigree
#
#now = time.time()
#print(ped.geneContribution(IDa = ped.getID(4), IDp = ped.getID(ped.lastGen())))
#time.time() - now
#
#
#
#
#
#
#import os
#os.getcwd()
#
#
#
#
#
#
#import matplotlib.pyplot as plt
#now = time.time()
#IBD = pop.realIBD(generation = pop.pedigree.lastGen())
#time.time() - now
#print(IBD)
#plt.matshow(IBD)
#
#now = time.time()
#print(pop.pedigree.computeNRM(ID = pop.pedigree.getID(pop.pedigree.lastGen())))
#time.time() - now
#
#
#F1 = Individual(name = 'F1')
#F2 = Individual(name = 'F2')
#I3 = Individual(name = 'I3', father = F1, mother = F2)
#I4 = Individual(name = 'I4', father = F1, mother = F2)
#I5 = Individual(name = 'I5', father = I3, mother = I4)
#I6 = Individual(name = 'I6', father = F1, mother = I5)
#
#newpop = Population(crossoverSimulator)
#newpop.individuals = {'F1' : F1, 'F2' : F2, 'I3' : I3, 'I4' : I4, 'I5' : I5, 'I6' : I6}
#pedigree = newpop.constructPedigree()
#pedigree.pairwiseGeneContribution(ida = 'F1' , idp = 'I6')
#pedigree.geneContribution(IDa, IDp)
#
#self = newpop.constructPedigree()
#newpop.pedigree.toDataFrame()
#self.ID
#
#
#
#
#
#
#
## Generate some data that where each slice has a different range
## (The overall range is from 0 to 2)
#data = np.random.random((4,10,10))
#data *= np.array([0.5, 1.0, 1.5, 2.0])[:,None,None]
#
## Plot each slice as an independent subplot
#fig, axes = plt.subplots(nrows=2, ncols=2)
#for dat, ax in zip(data, axes.flat):
#    # The vmin and vmax arguments specify the color limits
#    im = ax.imshow(dat, vmin=0, vmax=2)
#
## Make an axis for the colorbar on the right side
#cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
#fig.colorbar(im, cax=cax)
#
#plt.show()







