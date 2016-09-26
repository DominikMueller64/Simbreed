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

## math
import numpy as np
import math
from collections import namedtuple, OrderedDict
import collections  # OrderedDict
from collections.abc import Mapping, Sequence  # Population
import pandas as pd
import bisect
from scipy.spatial.distance import correlation, cosine, cdist


from sortedcollections import ValueSortedDict
from sortedcontainers import SortedDict
## tools
from itertools import chain, combinations, count
import itertools  # isSorted
import copy
import warnings
import time
## random numbers
import random
#
from numba import jit

TheDict = (SortedDict, ValueSortedDict, OrderedDict)[1]

#from simbreed.Functions import *



class Genome:

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

        self._numChrom = len(chromLengths)
        self._chromLengths = chromLengths
        self.positions = [np.asarray(a=pos, dtype=np.float64) for pos in positions]
        self.locusNames = locusNames

    @property
    def numChrom(self):
        return self._numChrom

    @property
    def chromLengths(self):
        return self._chromLengths

    @property
    def positions(self):
        return self.__positions

    @positions.setter
    def positions(self, value):
        abs_tol = 1e-6

        if not len(value) == self._numChrom:
            raise ValueError(
                "The length of 'positions' (%d) does not match the number of "
                "chromosomes (%d)." % (len(value), self._numChrom))

        for ix, pos in enumerate(value):
            if not isSorted(pos, abs_tol=0.0):
                raise ValueError("Positions on chromosome %d are not sorted "
                                 % (ix))
            # We can now assume pos to be sorted
            if pos[0] < 0.0:
                raise ValueError("The first position on chromosome %d (%f) "
                                 "must be >= 0.0." % (ix, pos[0]))

            if pos[-1] > self._chromLengths[ix]:
                raise ValueError("The last position on chromosome %d (%f), "
                                 "must be <= %f."
                                 % (ix, pos[0], self._chromLengths[ix]))

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

        if self.__positions is None:
            raise AttributeError("Before setting 'locusNames', 'positions' "
                                 "must be specified.")

        if not len(value) == self._numChrom:
            raise ValueError(
                "The length of 'locusNames' (%d) does not match the number of "
                "chromosomes (%d)." % (len(value), self._numChrom))

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

    def __init__(self, skeleton, ID=None, father=None,
                 mother=None, flesh=None, parent=None,
                 individual=None):
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
        self.individual = individual

    @property
    def ID(self):
        return self._ID

    @ID.setter
    def ID(self, value):
        self._ID = value

    @property
    def father(self):
        return self._father

    @father.setter
    def father(self, value):
        if value is None:
            self._father = value
            return None
        elif not isinstance(value, Gamete):
            raise ValueError("Attribute 'father' must be a 'Gamete'.")
        self._father = value

    @property
    def mother(self):
        return self._mother

    @mother.setter
    def mother(self, value):
        if value is None:
            self._mother = value
            return None
        elif not isinstance(value, Gamete):
            raise ValueError("Attribute 'mother' must be a 'Gamete'.")
        self._mother = value

    @property
    def skeleton(self):
        return self._skeleton

    @skeleton.setter
    def skeleton(self, value):
        for ix, chrom in enumerate(value):
            if not set(chrom.keys()).issuperset(['lineages', 'locations',
                                                 'Xlocations', 'parentStart']):
                raise ValueError("The provided 'skeleton' on chromosome %d "
                                 "does not contain the necessary keys." % (ix))
        self._skeleton = value

    @property
    def flesh(self):
        if self._flesh is None:
            self.enflesh()
        return self._flesh

    @flesh.setter
    def flesh(self, value):
        # Check if flesh fits to the skeletton.
        # This might be a later performance bottleneck and direct access to
        # __flesh should be pursued in such circumstances.
        if value is None:
            self._flesh = value
            return None
        elif not len(value) == len(self.skeleton):
            raise ValueError("The length of 'flesh' (%d) does not match the "
                             "length of 'skeleton' (%d) in gamete %s."
                             % (len(value), len(self.skeleton), str(self.ID)))
        self._flesh = value

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, value):
        if value is None:
            self._parent = value
            return None
        elif not isinstance(value, Individual):
            raise ValueError("Attribute 'parent' must be an 'Individual'.")
        self._parent = value


    def enflesh(self):
        """Add Genotypes. """
        population = self.individual.population
        genome = population.genome
        fleshDict = population.fleshDict
        numChrom = genome.numChrom
        flesh = [None] * numChrom
        for chrom in range(numChrom):
            positions = genome.positions[chrom]
            lineages = self.skeleton[chrom]['lineages']
            locations = self.skeleton[chrom]['locations']
#            tmp = [None] * len(positions)
            tmp = np.empty(len(positions))
            ixl = 0
            lo = 0
            for loc, lin in zip(locations, lineages):
                ixr = bisect.bisect_left(a=positions, x=loc, lo=ixl)
#                ixr = myBisect_left(a=positions, x=loc, lo=ixl)
                tmp[ixl:ixr] = fleshDict[lin][chrom][ixl:ixr]
                ixl = ixr

            flesh[chrom] = tmp

        self._flesh = flesh




# %%


class FounderGamete(Gamete):

    """FounderGamete

    Subclasses ``Gamete``

    Instance variables
    ------------------
    lineage : lineage of the founder gamete

    The user does not need to manually instantiate new FounderGametes.

    """

    def __init__(self, chromLengths, ID=None, flesh=None, individual=None):
        """Instantiate a new FounderGamete

        Parameters
        ----------
        lineage : integer
          The lineage of a founder gamete is an integer identifier labelling
          the gamete. This is used to track the origin of chromosomal
          segments in later generations to founder gametes.

        """
        lineage = getNewID()
#        skeleton = [{'lineages': [lineage],
#                     'locations': [chromLen],
#                     'Xlocations': [None],
#                     'parentStart': None} for chromLen in chromLengths]

        skeleton = [{'lineages': np.array([lineage], dtype=np.int64),
                     'locations': np.array([chromLen], dtype=np.float64),
                     'Xlocations': [None],
                     'parentStart': None} for chromLen in chromLengths]


        Gamete.__init__(self, ID=ID, skeleton=skeleton, flesh=flesh,
                        individual=individual)
#        super().__init__(self, ID=ID, skeleton=skeleton, flesh=flesh)

        self.__lineage = lineage

    @property
    def lineage(self):
        return self.__lineage

# %%









# %%


class Individual:

    def __init__(self, population, parents, generation,
                 DH, F, ID=None, reproGroup=None):

        self.population = population
        self.parents = parents
        self.generation = generation
        self.DH = DH
        self.F = F
        self.ID = ID
        self.reproGroup = reproGroup
        self.isRoot = False


        self.reproGroup = ''.join(str(p) for p in parents)

        self.creationTime = time.time()


#        self.children = set()

#        for parent in self.parents:
#            if not parent == 'NA':
#                self.population[parent].addChild(self)

    @property
    def population(self):
        return self._population

    @population.setter
    def population(self, value):
        self._population = value

    @property
    def ID(self):
        return self._ID

    @ID.setter
    def ID(self, value):
        self._ID = value

    @property
    def parents(self):
        return self._parents

    @parents.setter
    def parents(self, value):
        self._parents = value

#    @property
#    def father(self):
#        return self._father
#
#    @father.setter
#    def father(self, value):
#        if value is None:
#            self._father = value
#            return None
#        elif not isinstance(value, Individual):
#            raise ValueError("Attribute 'father' must be "
#                             "an 'Individual'.")
#        self._father = value
#
#    @property
#    def mother(self):
#        return self._mother
#
#    @mother.setter
#    def mother(self, value):
#        if value is None:
#            self._mother = value
#            return None
#        elif not isinstance(value, Individual):
#            raise ValueError("Attribute 'mother' must be "
#                             "an 'Individual'.")
#        self._mother = value

    @property
    def founder(self):
        return sum(map(lambda x: x == 'NA', self.parents)) == 2

    @property
    def semifounder(self):
        temp = sum(map(lambda x: x == 'NA', self.parents))
        return  0 < temp < len(self.parents)

    @property
    def generation(self):
        return self._generation

    @generation.setter
    def generation(self, value):
        self._generation = value


    @property
    def gametes(self):
        return self._gametes

    @gametes.setter
    def gametes(self, value):
        self._gametes = value

    @property
    def DH(self):
        return self._DH

    @DH.setter
    def DH(self, value):
        self._DH = value

    @property
    def F(self):
        return self._F

    @F.setter
    def F(self, value):
        self._F = value

    def __str__(self):
        return str(self.ID)

    def __getitem__(self, item):
        return self.gametes[item]

    def __gt__(self, I2):
        try:
            sg = self.generation
            I2g = I2.generation
        except (TypeError, AttributeError):  # If either is None or absent.
            return True
        else:
            if sg == I2g:
                return self.creationTime > I2.creationTime
            else:
                return sg > I2g

#    def addChild(self, childID):
#        """ Add a child. """
#        self.children.add(childID)
#
#    def hasChild(self, childID=None):
#        """ Has child/children?
#
#        Parameters
#        ----------
#        child : Entry (default=None)
#
#        If child is None, it is tested if there is any child at all.
#
#        """
#        if childID is not None:
#            return childID in self.children
#        else:
#            return bool(self.children)

    def enskelet(self):
        chromLengths = self.population.chromLengths
        crossoverSimulator = self.population.crossoverSimulator
        if self.DH:  # Only one gamete is produced.
            if not checkEqual(self.parents):
                raise RuntimeError("'father' is not identical to "
                                   "mother for a doubled haploid.")
            if self.parents[0] == 'NA':
                gametes = [FounderGamete(chromLengths=chromLengths,
                                         individual=self)] * 2
            else:
                parent = self.population[self.parents[0]]
                skeleton = simulateMeiosis(parent=parent,
                                           crossoverSimulator=crossoverSimulator)
                gametes = [Gamete(skeleton=skeleton, ID=getNewID(),
                                  father=parent.gametes[0],
                                  mother=parent.gametes[1],
                                  individual=self)] * 2
        else:
            gametes = [None] * 2
            for i, parent in enumerate(self.parents):
                if parent == 'NA':
                    gametes[i] = FounderGamete(chromLengths=chromLengths,
                                               individual=self)
                else:
                    parent = self.population[parent]
                    skeleton = simulateMeiosis(parent=parent,
                                               crossoverSimulator=crossoverSimulator)
                    gametes[i] = Gamete(skeleton=skeleton, ID=getNewID(),
                                        father=parent.gametes[0],
                                        mother=parent.gametes[1],
                                        individual=self)

        self.gametes = gametes

#    @abstractmethod
#    def _update(self):
#        pass

class Root:

    def __init__(self):
        self.ID = 'NA'
        self.parents = ()
        self.generation = -1
        self.isRoot = True
        self.founder = False
        self.semifounder = False
        self.children = set()


    def __gt__(self, I2):
        try:
            return self.generation > I2.generation
        except (TypeError, AttributeError):  # If either is None or absent.
            return True



# %%

#class populationIndividual(Individual):
#
#    def __init__(self, population, father, mother, generation, DH,
#                 ID=None):
#        super().__init__(population=population,
#                         father=father, mother=mother,
#                         generation=generation, DH=DH, F=None, ID=ID)
#
#    def enskelet(self):
#        if self.DH:  # Only one gamete is produced.
#            if self.father is not self.mother:
#                raise RuntimeError("'father' is not identical to "
#                                   "mother for a doubled haploid.")
#            gametes = [simulateMeiosis(self.father, self.crossoverSimulator)] * 2
#        else:
#            gametes = [None] * 2
#            for i, parent in enumerate((self.father, self.mother)):
#                gametes[i] = simulateMeiosis(parent, self.crossoverSimulator)
#
#        self.gametes = gametes
#
#
#class populationFounderIndividual(Individual):
#
#    def __init__(self, population, DH, F, ID=None):
#
#        super().__init__(population=population,
#                         father=None, mother=None,
#                         generation=0, DH=DH, F=F, ID=ID)
#    def enskelet(self):
#
#        if self.DH:  # Only one gamete is produced.
#            gametes = [FounderGamete(chromLengths=self.population.chromLengths)] * 2
#        else:
#            gametes = [None] * 2
#            for i in range(2):
#                gametes[i] = FounderGamete(chromLengths=self.population.chromLengths)
#        self.gametes = gametes



#    """Individual
#
#    Instance variables
#    ------------------
#    ID : Identifier
#
#    father : father individual
#
#    motehr : mother individual
#
#    gametes : gametes carried by the individual
#
#
#    """
#
#    def __init__(self, ID=None, father=None, mother=None):
#        """Instantiate an Individual
#
#        Parameters
#        ---------
#        ID : integer or string
#          A unique identifier
#        father : Individual
#          The father of the individual
#        mother : Individual
#          The mother of the individual
#        gametes : list
#          A list with the two Gametes carried by the individual
# %%










# %%


class pedPop(TheDict):

    @classmethod
    def fromPedigree(cls, ID, father, mother, DH=None, F=None):
        """Construct a 'Pedigree' object from a known pedigree.

        Parameters
        ----------
        ID : list
          List with unique identifiers.
        fatherID : list
          List with IDs of the male parents. Must be 'None' if unknown.
        motherID: list
          List with IDs of the female parents. Must be 'None' if unknown.
        DH: list (default: there are no doubled haploids)
          List with booleans indicating the if the respective individual is
          a doubled haploid. Ignored for founders.
        F: dictionary (default: all founders are outbred)
          Dictionary specifing inbreeding coefficients for founders

        This function should be employed by the user to buil a 'Pedigree'
        object for simulation from a previously known pedigree.

        If founders are inbreed, this can be specified in the
        'F' argument. Founders can per definition be no
        doubledHaploids, but can be completely inbreed, which is equivalent
        for founders.

        Generation information will automatically be inferred by the given
        relationships between individuals. Individuals having both parents
        unknown are treated as founder individuals.

        Having one parent known and another one unknown is untested and might
        lead to undefined behaviour. However, this case is rather irrelevant
        for simulation studies.

        Relationships between founder individuals cannot be accounted for.

        The pedigree will automatically be extended and sorted if necessary.
        """

        def checkIDNone(ID):
            """ Check if any ID is None. """
            if any(i is None for i in ID):
                raise ValueError("IDs cannot be None.")

        def checkIDUnique(ID):
            """ Check if all IDs are unique. """
            if not allUnique(ID):
                raise ValueError("Not all IDs are unique")

        def checkParentNone(parent):
            """ Check if parent is None. """
            if any(i is None for i in parent):
                raise ValueError("None is not allowed as parent ID, use 'NA'.")

        def checkSameLength(*args):
            """ Check if all given arguments have the same length. """
            if not len(set(map(len, args))) == 1:
                raise ValueError("Arguments do not have the same length.")

        def checkDHSelfing(DH, fatherID, motherID):
            """Check if all DH are selfing progeny."""
            if any(not f == m for d, f, m in
                   zip(DH, fatherID, motherID) if d):
                raise ValueError("Not all doubled haploids are "
                                 "selfing progeny.")
#        def checkF():
#            """Check if inbreeding coefficients of all founders are provided."""
#            founderID = [i for i in self.ID if self.isFounder(i)]
#            diff = set(founderID) - self.F.keys()
#            if diff:
#                raise ValueError(
#                        "If 'F' is provided, inbreeding "
#                        "coefficients for all founder individuals must "
#                        "be specified. They are missing for the following "
#                        "individuals:\n\n " + '{}'.format(diff))

       # Construct Pedgiree
        checkIDNone(ID=ID)
        checkIDUnique(ID=ID)
        checkSameLength(ID, father, mother)
        checkParentNone(parent=father)
        checkParentNone(parent=mother)
        if DH is not None:
            checkSameLength(ID, DH)
            checkDHSelfing(DH=DH, fatherID=father, motherID=mother)
        else:
            DH = [False] * len(ID)
#        root = Root()
        graph = dict()
#        graph.append(root)
        for i, f, m, d in zip(ID, father, mother, DH):

            if not f == m == 'NA':  # not founder
                temp = None
            else:
                try:
                    temp = F.get(i, 0.0)  # Try to retrieve.
                except AttributeError:
                    temp = 0.0  # If F is not specified.

            ind = Individual(population=graph, parents=(f, m),
                             generation=None, DH=d, F=temp, ID=i)
            graph.update({i: ind})

        if isCyclic(graph):  # Check for loops.
            raise ValueError("The Pedigree contains loops.")
        sort(graph)  # Sort the graph. Only possible for an ordinary dict.
#        graph.sort()
        return cls(graph)


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.append(Root())

    def __copy__(self):
        newone = type(self)()
        newone.update(self)
        newone.__dict__.update(self.__dict__)
        return newone

    def append(self, value, reproGroup=None):
        """ Append an entry. """
        self.update({value.ID: value})

    def extend(self, values):
        """ Extend by appending several entries. """
        for value in values:
            self.append(value)

    def allID(self):
        """ Return all IDs, except for the root of the population. """
        temp = list(self.keys())
        try:
            temp.remove('NA')
        except ValueError:
            pass
        return temp

    def getID(self, generation=None):
        """ Return all IDs of a specific generation. """
        if generation is None:
            return self.allID()
        return tuple(i for i, v in self.items() if
                getattr(v, 'generation', None) == generation)

    def getFounderID(self):
        """ Return all IDs of founder individuals. """
        return self.getID(generation=0)

    # TODO: Change to getInd.
    def returnInd(self, ID):
        """ Return all Indivdiuals with given IDs. """
        if isSingle(ID):
            ID = (ID,)
        return tuple(self[i] for i in ID)

    def lastGen(self):
        """ Return the number of the last generation. """
        return max(e.generation for e in self.values())
    
    def getLastID(self):
        """ Return all IDs of the last generation. """
        return self.getID(generation=self.lastGen())

    def _subset(self, Ind):
        """ For internal use only. """
        if isSingle(Ind):
            Ind = (Ind,)
        seen = set([self['NA']])
        for ind in Ind:
            if ind not in seen:
                stack = [ind]
                while stack:
                    v = stack.pop()
                    seen.add(v)
                    stack.extend(set(self[p] for p in v.parents) - seen)
#        seen.discard(None)
        pop = pedPop()
        for v in self.values():  # Preserve the order.
            if v in seen:
                pop.append(v)
        return pop

    def subset(self, ID):
        """ Subset the population and return a copy.

        This leaves the original population unchanged.
        """
        if isSingle(ID):
            ID = (ID,)
        seen = set(['NA'])
        for i in ID:
            if i not in seen:
                stack = [i]
                while stack:
                    v = self[stack.pop()]
                    seen.add(v.ID)
                    stack.extend(set(v.parents) - seen)
#        seen.discard('NA')
        pop = pedPop()
        for v in self.values():  # Preserve the order.
            if v.ID in seen:
                pop.append(v)
        return pop

#    def _scramble(self):
#        """ Scramble the populaiton. Only for testing purposes. """
#        temp = [x.generation for x in self.values()]
#        random.shuffle(temp)
#        for ix, ind in enumerate(self.values()):
#            ind.generation = temp[ix]
#        ID = list(self.keys())
#        random.shuffle(ID)
#        pop = type(self)()
#        for i in range(len(self)):
#            pop.append(self[ID[i]])
#        self.clear()
#        self.update(pop)

#    def isSorted(self):
#        """Check if the pedigree is sorted."""
#        seen = set(['NA'])
#        return not any(any(p not in seen for p in v.parents) or
#                       seen.add(v.ID) for v in self.values())

    def _createIndividual(self, ID=None, father='default', mother='default',
                         generation=None, DH=False, F=None):

        if father == 'default':
            father = self['NA']
        if mother == 'default':
            mother = self['NA']

        if ID is None:
            ID = getNewID()

        if generation is not None:
            if any(generation <= g for g in (father.generation, mother.generation)):
                raise ValueError("The provided 'generation' value is "
                                 "not possible.")
        else:
            generation = max(father.generation, mother.generation) + 1

        self.append(Individual(population=self,
                               parents=(father.ID, mother.ID),
                               generation=generation, DH=DH, F=F, ID=ID))



    def createFounders(self, n, ID=None, F=None, returnID=False):
        if ID is None:
            ID = getNewID(n=n)
        for i in range(n):
            self._createIndividual(ID=ID[i], generation=0, DH=None,
                                  F=F[i] if F is not None else 0.0)
        if returnID:
            return ID


    def _cross(self, ID, father, mother, generation=None):
        """ Internal use only. """
        self._createIndividual(ID=ID, father=father, mother=mother,
                              generation=generation, DH=False)

    def cross(self, fatherID, motherID, generation=None,
              ID=None, returnID=False):
        """ Cross two individuals. """
        if ID is None:
            ID = getNewID()
        self._cross(ID=ID, father=self[fatherID], mother=self[motherID],
                    generation=generation)
        if returnID:
            return ID


    def _selfing(self, ID, parent, generation):
        self._cross(ID=ID, father=parent, mother=parent,
                    generation=generation)

    def _DH(self, ID, parent, generation):
        self._createIndividual(ID=ID, father=parent, mother=parent,
                              generation=generation, DH=True)

    def _getReproGroups(self, Ind, reproGroupGen):
        if reproGroupGen < 0:
            raise ValueError("'reproGroupGen' must be >= 0.")

        groupDict = dict()
        for ind in Ind:
            stack = [ind]
            while stack:
                vertex = stack.pop()
                if vertex.generation < reproGroupGen:
                    next
                elif vertex.generation == reproGroupGen:
                    try:
                        groupDict[vertex.reproGroup].add(ind.ID)
                    except KeyError:
                        groupDict[vertex.reproGroup] = set([ind.ID])
                else:
                    stack.extend(list(self[v] for v in set(vertex.parents)))
        return groupDict


    def getReproGroups(self, reproGroupGen, ID=None, generation=None):
        ID = IDGenerationHandler(idGetter=self.getID, ID=ID,
                         generation=generation)
        Ind = tuple(self[i] for i in ID)
        return self._getReproGroups(Ind=Ind, reproGroupGen=reproGroupGen)

    def randomMating(self, size, allowSelfing=True, ID=None,
                     returnID=False, reproGroupGen=None):
        """Perform random mating."""
        lastgen = self.lastGen()
        if ID is None:
            ID = self.getID(lastgen)
        else:
            if not set(self.allID()).issuperset(set(ID)):
                raise ValueError("'ID' is not a subset of the IDs present.")

        if not allowSelfing and len(ID) < 2:
            raise ValueError(
                "It is not possible to have 'allowSelfing = False' if there "
                "are less than 2 individuals in the population.")

        Ind = tuple(self[i] for i in ID)
        if reproGroupGen is not None:
            groupDict = self._getReproGroups(Ind=Ind, reproGroupGen=reproGroupGen)
            for group in groupDict.values():
                self.randomMating(size=size, allowSelfing=allowSelfing,
                                  ID=group)
            return  # to break out

        newID = getNewID(n=size, sequence=True)
        for i in newID:
            father = random.choice(Ind)
            mother = random.choice(Ind)
            if not allowSelfing:
                while father is mother:
                    mother = random.choice(Ind)
            self._cross(ID=i, father=father, mother=mother)

        if returnID:
            return newID
    # TODO: methods inbreeding and dH.. are so similar.... code duplication!!!
    def inbreeding(self, size=1, ID=None, returnID=False):
        """Perform self-fertilization.

        This naturally does not need to account for reproductive groups.
        """
        lastgen = self.lastGen()
        if ID is None:
            ID = self.getID(lastgen)
        else:
            if not set(self.allID()).issuperset(set(ID)):
                raise ValueError("'ID' is not a subset of the IDs present.")
        Ind = tuple(self[i] for i in ID)
        newID = getNewID(n=size * len(Ind), sequence=True)
        ct = itertools.count()
        for parent, _ in itertools.product(Ind, range(size)):
            self._selfing(ID=newID[next(ct)], parent=parent,
                          generation=lastgen + 1)
#            for parent in Ind:
#                for _ in range(size):
        if returnID:
            return newID

    def doubledHaploidization(self, size=1, ID=None, returnID=False):
        """Produce doubled haploids.

        This naturally does not need to account for reproductive groups.
        """

        lastgen = self.lastGen()
        if ID is None:
            ID = self.getID(lastgen)
        else:
            if not set(self.allID()).issuperset(set(ID)):
                raise ValueError("'ID' is not a subset of the IDs present.")
        Ind = tuple(self[i] for i in ID)
        newID = getNewID(n=size * len(Ind), sequence=True)
        ct = itertools.count()
        for parent, _ in itertools.product(Ind, range(size)):
#        for parent in Ind:
#            for _ in range(size):
            self._DH(ID=newID[next(ct)], parent=parent,
                     generation=lastgen + 1)
        if returnID:
            return newID

    def roundRobin(self, size=1, ID=None, returnID=False,
                   reproGroupGen=None):
        """Perform mating according to a round robin design."""
        lastgen = self.lastGen()
        if ID is None:
            ID = self.getID(lastgen)
        else:
            if not set(self.allID()).issuperset(set(ID)):
                raise ValueError("'ID' is not a subset of the IDs present.")
        if len(ID) < 3:
            raise ValueError("RoundRobin with less than 3 individuals "
                             "is not possible.")
        Ind = tuple(self[i] for i in ID)

        if reproGroupGen is not None:
            groupDict = self._getReproGroups(Ind=Ind, reproGroupGen=reproGroupGen)
            for group in groupDict.values():
                self.roundRobin(size=size, ID=group)
            return  # to break out

        newID = getNewID(n=size * len(Ind), sequence=True)  # Should be the correct number.
        ct = itertools.count()
        for (father, mother), _ in itertools.product(zip(Ind, Ind[1:] + Ind[:1]), range(size)):
            self._createIndividual(ID=newID[next(ct)], father=father,
                               mother=mother)
#        for id1 in curID:
#            for id2 in curID[1:] + curID[:1]:
#                for i in range(size):
        if returnID:
            return newID


    def synthetic(self, size=1, ID=None, reciprocals=False, returnID=False,
                  reproGroupGen=None):
        """Perform mating according to a synthetic design.

        Note: Making reciprocals does not really make sense in a model where
        traits are only influenced by nuclear DNA.
        """
        lastgen = self.lastGen()
        if ID is None:
            ID = self.getID(lastgen)
        else:
            if not set(self.allID()).issuperset(set(ID)):
                raise ValueError("'ID' is not a subset of the IDs present.")
        n = len(ID)
        if n < 2:
            raise ValueError("A Synthetic with less than 2 individuals "
                             "is not possible.")

        nmat = int(n * (n - 1) / (2 - reciprocals))  # Number of pairs.
        Ind = tuple(self[i] for i in ID)

        if reproGroupGen is not None:
            groupDict = self._getReproGroups(Ind=Ind, reproGroupGen=reproGroupGen)
            for group in groupDict.values():
                self.synthetic(size=size, reciprocals=reciprocals)
            return  # to break out

        newID = getNewID(n=size * nmat, sequence=True)
        ct = itertools.count()
        for ix1, parent1 in enumerate(Ind):
            for (ix2, parent2), _ in itertools.product(enumerate(Ind[ix1 + 1:],
                    start = ix1 + 1), range(size)):
                self._createIndividual(ID=newID[next(ct)], father=parent1,
                                       mother=parent2)
                if reciprocals:
                    self._createIndividual(ID=newID[next(ct)], father=parent2,
                                          mother=parent1)
#            for ix2, id2 in enumerate(curID[ix1 + 1:], start = ix1 + 1):
#                self.createEntry(getNewID(), id1, id2, curgen + 1, False, None)
#                if reciprocals:
#                    self.createEntry(getNewID(), id2, id1, curgen + 1, False, None)
        if returnID:
            return newID


    def _tabularMethod(self, IDa, IDb):
        # Subset pedigree for minimal computation time
        pop = self.subset(ID=set().union(*(IDa, IDb)))
        ID = pop.allID()
        x = dict(zip(ID, itertools.count()))
        indicesA = tuple(x[i] for i in IDa)
        if IDb is IDa or IDb == IDa:
            indicesB = indicesA
        else:
            indicesB = tuple(x[i] for i in IDb)
        mat = np.identity(len(ID), dtype=np.float_)
        for i, v in enumerate(ID):
            v = self[v]
            f = v.parents[0]
            m = v.parents[1]
            # Diagonal
            if v.founder:
                if v.DH:
                    mat[i, i] = 2.0
                else:
                    mat[i, i] = 1.0 + v.F
#            elif not any(x == 'NA' for x in v.parents):
            elif f != 'NA' and m != 'NA':
                if v.DH:
                    mat[i,i] = 2.0
                else:
                    mat[i, i] = 1.0 + 0.5 * mat[x[f], x[m]]
            # Off-diagonal
            temp = 0.0
            for p in (f, m):
                if p != 'NA':
                    temp += 0.5 * mat[x[p], 0:i]
            mat[i, 0:i] = mat[0:i, i] = temp

        return mat[np.ix_(indicesA, indicesB)]

    def _LDLMethod(self, IDa, IDb):
        pop = self.subset(ID=set().union(*(IDa, IDb)))
        pop.pop('NA')
        ID = pop.allID()
        x = dict(zip(ID, itertools.count()))
        n = len(pop)
        L = np.identity(n=n, dtype=np.float_)
        for j, vj in enumerate(pop.values()):
            for i, vi in enumerate(list(pop.values())[j + 1:n], start=j + 1):
                fi = vi.parents[0]
                mi = vi.parents[1]
                if fi != 'NA':
                    if mi != 'NA':
                        L[i, j] = 0.5 * (L[x[fi], j] + L[x[mi], j])
                    else:
                        L[i, j] = 0.5 * L[x[fi], j]
                else:
                    if mi != 'NA':
                        L[i, j] = 0.5 * L[x[mi], j]
                    else:
                        L[i, j] = 0

        # Relies on the first individual beeing a founder.
        D = np.zeros(n)
        F = np.zeros(n)
        for j, vj in enumerate(pop.values()):
            fj = vj.parents[0]
            mj = vj.parents[1]
            if fj != 'NA':
                if mj != 'NA':
                    if vj.DH:
                        D[j] = 1.0
                    else:
                        D[j] = 0.5 - 0.25 * (F[x[fj]] + F[x[mj]])
                else:
                    D[j] = 0.75 - 0.25 * F[x[fj]]
            else:
                if mj != 'NA':
                    D[j] = 0.75 - 0.25 * F[x[mj]]
                else:  # The individual must be a founder!
                    D[j] = 1.0 + vj.F
#                    D[j] = 1 # normal code
            F[j] = np.dot(np.square(L[j,]), D) - 1.0

        rtuple = collections.namedtuple('rtuple', ['L', 'D', 'F', 'ID'])
        return rtuple(L=L, D=D, F=F, ID=ID)
#        return {'L': L, 'D': D, 'F': F, 'names' = tuple(i.ID for i in Ind)}

    def makeA(self, ID=None, generation=None,
              IDb=None, generationB=None, method='tabular'):
        ID, IDb =  dualIDGenerationHandler(self.getID, ID, generation, IDb, generationB)
        if method == 'tabular':
            mat = self._tabularMethod(IDa=ID, IDb=IDb)
        elif method == 'LDL':
            LDL = self._LDLMethod(IDa=ID, IDb=IDb)
            L = LDL.L
            mat = (L * LDL.D).dot(L.T)
            x = dict(zip(LDL.ID, itertools.count()))
            mat = mat[np.ix_(tuple(x[i] for i in ID), tuple(x[i] for i in IDb))]
        else:
            raise ValueError("'method' must be either 'tabular' "
                             "or 'LDL'.")
        rtuple = collections.namedtuple('rtuple', ['matrix', 'rowID', 'colID'])
        return rtuple(matrix=mat, rowID=ID, colID=IDb)

    def makeAinv(self, ID=None, generation=None,
                 IDb=None, generationB=None, method='tabular'):
        ID, IDb =  dualIDGenerationHandler(self.getID, ID, generation,
                                           IDb, generationB)

        if method == 'tabular':
            mat = np.linalg.inv(self._tabularMethod(IDa=ID, IDb=IDb))
        elif method == 'LDL':
            LDL = self._LDLMethod(IDa=ID, IDb=IDb)
            Linv = np.linalg.inv(LDL.L)
            mat = (Linv.T * (1.0 / LDL.D)).dot(Linv)
            x = dict(zip(LDL.ID, itertools.count()))
            mat = mat[np.ix_(tuple(x[i] for i in ID), tuple(x[i] for i in IDb))]
            # Cholesky approach not faster even for large (5000) matrices.
#            n = A.shape[0]
#            L = np.linalg.cholesky(A)
#            Y = scipy.linalg.solve_triangular(L, np.eye(n), lower=True, check_finite=False)
#            X = scipy.linalg.solve_triangular(L, Y, trans=1, lower=True, overwrite_b=True, check_finite=False)
        else:
            raise ValueError("'method' must be either 'tabular' "
                             "or 'LDL'.")
        rtuple = collections.namedtuple('rtuple', ['matrix', 'rowID', 'colID'])
        return rtuple(matrix=mat, rowID=ID, colID=IDb)

    def inbreedingCoef(self, ID=None, generation=None, method='tabular'):
        ID = IDGenerationHandler(self.getID, ID, generation)
        if method == 'tabular':
            A = self._tabularMethod(IDa=ID, IDb=ID)
            F = (np.diag(A) - 1.0).tolist()
        elif method == 'LDL':
            LDL = self._LDLMethod(IDa=ID, IDb=ID)
            x = dict(zip(LDL.ID, itertools.count()))
            F = LDL.F.tolist()
            F = [F[x[i]] for i in ID]
        else:
            raise ValueError("'method' must be either 'tabular' "
                             "or 'LDL'.")
        rtuple = collections.namedtuple('rtuple', ['F', 'ID'])
        return rtuple(F=F, ID=ID)

    def mendelianSamplingVariance(self, ID=None, generation=None):
        ID = IDGenerationHandler(self.getID, ID, generation)
        LDL = self._LDLMethod(IDa=ID, IDb=ID)
        x = dict(zip(LDL.ID, itertools.count()))
        D = LDL.D.tolist()
        D = [D[x[i]] for i in ID]
        rtuple = collections.namedtuple('rtuple', ['D', 'ID'])
        return rtuple(D=D, ID=ID)


    # TODO: Add support for founder individuals
    # TODO: update
    def makeD(self, ID=None, generation=None,
              IDb=None, generationB=None, method='tabular'):
        raise NotImplementedError("makeD not yet implemented.")

        ID = IDGenerationHandler(idGetter=self.getID, ID=ID, generation=generation)
        A = self.makeA(ID = self.allID(), method = method)
        x = dict(zip(self, itertools.count()))
        if IDb is None and generationB is None:
            mat = np.zeros(shape=(len(ID), len(ID)), dtype=np.float_)
            for ixa, ida in enumerate(ID):
                va = self[IDa]; fa = va.fatherID; ma = va.motherID
                if va.founder:
                    raise ValueError("Individual %s is a founder an cannot be "
                                     "included." % (IDa))
                for ixb, idb in enumerate(ID[ixa + 1:], start=ixa + 1):
                    vb = self[IDb]; fb = vb.fatherID; mb = vb.motherID
                    if vb.founder:
                            raise ValueError("Individual %s is a founder an cannot be "
                                             "included." % (IDb))
                    dom = 0.25 * (A[x[fa], x[fb]] + A[x[ma], x[mb]] +
                                  A[x[fa], x[mb]] + A[x[ma], x[fb]])
                    mat[ixa, ixb] = mat[ixb, ixa] = dom

                mat[ixa, ixa] = 0.25 * (A[x[fa], x[fa]] + A[x[ma], x[ma]] +
                                        2.0 * A[x[fa], x[ma]])
            pkg = {'matrix': mat, 'rownames': ID, 'colnames': ID}
        else:
            IDb = IDGenerationHandler(self.getID, IDb, generationB)
            mat = np.zeros(shape=(len(ID), len(IDb)), dtype=np.float_)
            for ixa, ida in enumerate(ID):
                va = self[IDa]; fa = va.fatherID; ma = va.motherID
                if va.founder:
                    raise ValueError("Individual %s is a founder an cannot be "
                                     "included." % (IDa))
                for ixb, idb in enumerate(IDb):
                    vb = self[IDb]; fb = vb.fatherID; mb = vb.motherID
                    if vb.founder:
                            raise ValueError("Individual %s is a founder an cannot be "
                                             "included." % (IDb))
                    mat[ixa, ixb] = 0.25 * (A[x[fa], x[fb]] + A[x[ma], x[mb]] +
                                            A[x[fa], x[mb]] + A[x[ma], x[fb]])
            pkg = {'matrix': mat, 'rownames': ID, 'colnames': IDb}

        return pkg

#        x = dict(zip(self, itertools.count()))  # Gives index in A, not D!
        for i, idi in enumerate(ID):
            vi = self[idi]; fi = vi.fatherID; mi = vi.motherID
            if vi.founder:
                raise ValueError("Individual %s is a founder an cannot be "
                                 "included." % (idi))
            for j, idj in enumerate(ID[i + 1:n], start=i + 1):
                vj = self[idj]; fj = vj.fatherID; mj = vj.motherID
                dom = 0.25 * (A[x[fi], x[fj]] + A[x[mi], x[mj]] +
                              A[x[fi], x[mj]] + A[x[mi], x[fj]])
                D[i, j] = D[j, i] = dom
            D[i, i] = 0.25 * (A[x[fi], x[fi]] + A[x[mi], x[mi]] +
                                 2.0 * A[x[fi], x[mj]])  # Check if necessary
        return D


    def _pairExpGeneContr(self, IDa, IDb = None):
        if IDb is None:
            IDb = IDa
        if IDb is IDa or IDb == IDa:
            return 0.0
        paths = []
        def myDFS(graph, start, end, path=[]):
            path = path + [start]
            if start == end:
                paths.append(path)
            for node in graph[start].parents:
                if node not in path:
                    myDFS(graph, node, end, path)

        myDFS(graph=self, start=IDb, end=IDa)
        return sum(map(lambda x: 0.5**(len(x) - 1), paths))

#    def _pairExpGeneContr(self, IDa, IDb=None):
#        if IDb is None:
#            IDb = IDa
#        ga = self[IDa].generation
##        if not ga < self[IDp].generation:
##            raise ValueError("%s is not older than %s" % (IDa, IDp))
#        qdict = collections.Counter([IDb])
#        contr = 0.0  # Contribution of IDa to IDp.
#        c = 0.5  # Weighting factor, halved for each step back.
#        while qdict:
#            nqdict = collections.Counter()
#            for i in qdict.keys():
#                vi = self[i]
#                if vi.generation > ga:  # Otherwise, parents cannot be target.
#                    for p in (vi.father, vi.mother):
#                        if p is None:
#                            continue
#                        if p.ID == IDa:  # Contribute c times number of paths.
#                            contr += c * qdict[i]
#                        elif p is not None:  # Step back.
#                            nqdict.update({p.ID: qdict[i]})
#            c *= 0.5
#            qdict = nqdict
#        return contr

    def expGeneContr(self, ID=None, generation=None,
                           IDb=None, generationB=None):
        if ID is not None and isSingle(ID) and isSingle(IDb):
            return self._pairExpGeneContr(IDa=ID, IDb=IDb)
        return matrixComputer(idGetter = self.getID,
                              function = self._pairExpGeneContr,
                              ID = ID, generation = generation,
                              IDb = IDb, generationB = generationB,
                              symmetric=False)


    def _pairIsAncestor(self, IDa, IDb=None):
        if IDb is None:
            return False
        # Depth-first search algorithm for an arbitrary graph.
        visited, stack = set(['NA']), [IDb]
        while stack:
            vertex = stack.pop()
            if vertex == IDa:
                return True
            if vertex not in visited:
                visited.add(vertex)
                stack.extend(set(self[vertex].parents) - visited)
        return False

#    def isRelated(self, ID1, ID2):
#        """ Test if two individuals are related by pedigree. """
#        return isAncestor(ID1, ID2) or isAncestor(ID2, ID1)

    def isAncestor(self, ID=None, generation=None,
                           IDb=None, generationB=None):
        if ID is not None and isSingle(ID) and isSingle(IDb):
            return self._pairIsAncestor(IDa=ID, IDb=IDb)

        return matrixComputer(idGetter=self.getID,
                              function=self._pairIsAncestor,
                              ID=ID, generation=generation,
                              IDb=IDb, generationB=generationB,
                              symmetric=False, dtype=np.bool)

    def recode(self):
        # Keys are immutable, so entries must be eventually deleted.
        keys = self.keys() - set(['NA'])
        x = dict(zip(keys, itertools.count(start=1)))
        x.update({'NA':'NA'})
        for i in keys:
            ni = x[i]
            ind = self.pop(i)
            ind.ID = ni
            ind.parents = tuple(x[p] for p in ind.parents)
            self.update({ni: ind})

    def numInd(self, generation):
        """ Return the number of individuals in a certain generation. """
        ct = itertools.count()
        for e in self.values():
            if getattr(e, 'generation', None) == generation:
                next(ct)
        return next(ct)
    
    
    def getDoubledHaploids(self, generation):
        return [i for i, v in self.items() if
                getattr(v, 'generation', None) == generation and
                getattr(v, 'DH', None)]

    def offspring(self, ID=None, generation=None):
        ID = IDGenerationHandler(idGetter=self.getID, ID=ID, generation=generation)
        if not set(self).issuperset(set(ID)):
            raise ValueError("Not all individuals in ID are contained in "
                             "the pedigree")
        ret = {i:set() for i in ID}
        keys = self.keys() - set(['NA'])
        for i in ret:
            for j in keys - set([i]) :
                if i in self[j].parents:
                    ret[i].add(j)

        return ret

    def countOffspring(self, ID=None, generation=None):
        off = offspring(ID=ID, generation=generation)
        return {i: len(v) for i, v in off.items()}


    def toDataFrame(self, showGeneration=True, showDH=True, showF=True):
        if not self:
            return pd.DataFrame()

        columns = ['ID', 'parents']
        colnames = ['ID', 'parentA', 'parentB']

        for temp in (columns, colnames):
            if showGeneration:
                temp += ['generation']
            if showDH:
                temp += ['DH']
            if showF:
                temp += ['F']

        temp = []
        for v in self.values():
            t = []
            for a in columns:
                if a == 'parents':
                    if v.isRoot:
                        for _ in range(2):
                            t.append(str(None))
                    else:
                        parents = getattr(v, a)
                        for parent in parents:
                            t.append(str(parent))
                else:
                    t.append(str(getattr(v, a, None)))
            temp.append(t)



#        temp = ([str(getattr(v, a, None)) for a in columns] for v in self.values())
        return pd.DataFrame(temp, columns=colnames)

    def __repr__(self):
        return str(self.toDataFrame())

    def write(self, file, sep=' ', na_rep='NA'):
        temp = self.toDataFrame()
        temp.to_csv(file, sep=sep, na_rep=na_rep, index=False)



#    def tabulateInd(self):
#        return tuple(map(self.numInd,
#                        sorted(set(e.generation for e in self.values()))))

#    def numGen(self):
#        """ Return the number of generations in the pedigree."""
#        return len(set(e.generation for e in self.values))



#    def isFounder(self, ID):
#        flag = True
#        if isSingle(ID):
#            ID = (ID,)
#            flag = False
#        if not set(self).issuperset(set(ID)):
#            raise ValueError("Not all individuals in ID are contained in "
#                             "the pedigree")
#        f = [getattr(self[i], 'founder', None) for i in ID]
#        return f if flag else f[0]






class skPop(pedPop):

    @classmethod
    def fromPedigree(cls, ID, father, mother, chomLengths,
                     crossoverSimulator, DH=None, F=None):
        pP = pedPop.fromPedigree(ID=ID, father=father,
                                 mother=mother, DH=DH, F=F)
        return skPop.fromPedPop(
                   population=pP,
                   chomLengths=chromLengths,
                   crossoverSimulator=crossoverSimulator)


    @classmethod
    def fromPedPop(cls, population, chromLengths,
                               crossoverSimulator):
        sP = cls(chromLengths=chromLengths, crossoverSimulator=crossoverSimulator)
        sP.update(population)
        for ind in population.values():
            if not ind.isRoot:
                ind.population = sP
                ind.enskelet()

        return sP

    def __init__(self, chromLengths, crossoverSimulator):
        super().__init__()
        self.chromLengths = chromLengths
        self.crossoverSimulator = crossoverSimulator

    def __copy__(self):
        newone = type(self)(chromLengths=self.chromLengths,
                            crossoverSimulator=self.crossoverSimulator)
        newone.update(self)
        newone.__dict__.update(self.__dict__)
        return newone
#    def fleshDecorator(oldMethod):
#        def newMethod(self, ID, DH, F):
#            oldMethod(ID=ID, DH=DH, F=F)
#            self[ID].enskelet()
#
#    createFounder = fleshDecorator(self.createFounder)

    def _createIndividual(self, ID=None, father='default', mother='default',
                         generation=0, DH=False, F=None):

        if father == 'default':
            father = self['NA']
        if mother == 'default':
            mother = self['NA']

        if ID is None:
            ID = getNewID()

        if generation is not None:
            if any(generation <= g for g in (father.generation, mother.generation)):
                raise ValueError("The provided 'generation' value is "
                                 "not possible.")
        else:
            generation = max(father.generation, mother.generation) + 1

        ind = Individual(population=self,
                         parents=(father.ID, mother.ID),
                         generation=generation, DH=DH, F=F, ID=ID)
        self.append(ind)
        ind.enskelet()


    def _realIBDGam(self, gam1, gam2):
        """Compute the realized IBD coefficient between two gametes."""
        gam1, gam2 = gam1.skeleton, gam2.skeleton
        chromLengths = self.chromLengths
        numChrom = len(chromLengths)
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


    def _pairRealIBD(self, IDa, IDb=None):
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
        ind1 = self[IDa]
        if IDb is None:
            ind2 = ind1
        else:
            ind2 = self[IDb]
        IBD = 0.25 * (self._realIBDGam(ind1.gametes[0], ind2.gametes[0]) +
                      self._realIBDGam(ind1.gametes[0], ind2.gametes[1]) +
                      self._realIBDGam(ind1.gametes[1], ind2.gametes[0]) +
                      self._realIBDGam(ind1.gametes[1], ind2.gametes[1]))
        return IBD

    def realIBD(self, ID=None, generation=None,
                      IDb=None, generationB=None):
        """Compute the realized IBD coefficient between all individuals in a
        generation

        Parameters
        ----------
        generation : integer (defaults to the last generation)
          The generation for which realized IBD coefficients should be
          calculated.

        """
        if ID is not None and isSingle(ID) and isSingle(IDb):
            return self._pairRealIBD(IDa=ID, IDb=IDb)
        return matrixComputer(idGetter=self.getID,
                              function=self._pairRealIBD,
                              ID=ID, generation=generation,
                              IDb=IDb, generationB=generationB)


    def realFounderGeneContrMatrix(self, ID=None, generation=None,
                                   includeUnknown=False):
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
        ID = IDGenerationHandler(self.getID, ID, generation)
        lookup = {}
        fix = collections.OrderedDict()
        ct = itertools.count()
        for i, v in self.items():
            if v.founder or (v.semifounder if includeUnknown else False):
                fix.update({i: next(ct)})
                for gamete in v.gametes:
                    try:
                        lin = getattr(gamete, 'lineage')
                    except AttributeError:
                        pass
                    else:
                        lookup[lin] = i
        mat = np.zeros(shape=(len(ID), len(fix)), dtype=np.float_)

        for ix, i in enumerate(ID):
            contr = self._realFounderGeneContr(i)
            for lin, v in contr.items():
                try:  # Possibility of founder alleles not present in founders.
                    f = lookup[lin]
                except KeyError:
                    pass
                else:
                    mat[ix, fix[f]] += v
        rtuple = collections.namedtuple('rtuple', ['matrix', 'rownames', 'colnames'])
        return rtuple(matrix=mat, rownames=tuple(ID), colnames=tuple(fix.keys()))

    def _realFounderGeneContr(self, ID):
        """Calculate realized gene contribution of founders to a single
        individual."""
        ind = self[ID]
        d = {}
        totalLen = 2.0 * sum(self.chromLengths)
        for gamete in ind.gametes:
            for ix, chrom in enumerate(gamete.skeleton):
                lineages = chrom['lineages']
                locations = [loc / totalLen for loc in chrom['locations']]
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


    def _realGeneContr(self, ID, maxGenBack):
        raise NotImplementedError("Not correctly implemented due to a "
                                  "subtle bug.")
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

        ind = self[ID]
        maxGenBack = min(maxGenBack, self[ID].generation)
        chromLengths = self.chromLengths
        allChunks = set()
        

        for chrom, chromLen in enumerate(chromLengths):
            chunkSet = set([Chunk(ind, ind.gametes[0], chrom, ind.parents[0],
                                  0.0, chromLen),
                            Chunk(ind, ind.gametes[1], chrom, ind.parents[1],
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
                    if parent == 'NA':  # reached end of pedigree, place back
                        newChunkSet.add(chunk)
                        continue

                    tmp = gamete.skeleton[chunk.chrom]  # info stored here
                    Xlocations = tmp['Xlocations'].copy()
                    grandparents = self[parent].parents
                    parentAtStart = tmp['parentStart']

                    parentalGametes = (gamete.father, gamete.mother)
                    cur_par = parentalGametes.index(parentAtStart)
#                    IX = iter(zip([0.0] + Xlocations,
#                                  Xlocations + [chromLen]))
#                    xol, xor = next(IX)
                    for xol, xor in list(zip(np.hstack(([0.0], Xlocations)),
                                        np.hstack((Xlocations, [chromLen])))):
                        if xor <= start:
                            continue
                        elif xol < start:
                            if start < xor <= stop:
                                nstart = start
                                nstop = xor
                            elif stop < xor:
                                nstart = start
                                nstop = stop
                        elif start <= xol < stop:
                            if start < xor <= stop:
                                nstart = xol
                                nstop = xor
                            elif stop < xor:
                                nstart = xol
                                nstop = stop
                        elif stop < xol:
                            break

                        print("start: %.2f \t stop: %.2f \t xol: %.2f \t "
                              "xor: %.2f" % (start, stop, xol, xor))
                        print("%.5s" % gamete.ID + "\tLevel " +
                              str(counter + 1) + "\tnew Chunk " +
                              "\tnstart:%.2f " % (nstart),
                              "\tnstop:%.2f" % (nstop))
                        print("\n")
                        newChunkSet.add(Chunk(gamete=parentalGametes[cur_par],
                                              individual=self[parent],
                                              chrom=chrom,
                                              parent=grandparents[cur_par],
                                              start=nstart,
                                              stop=nstop))

                        cur_par = 1 - cur_par

                chunkSet = newChunkSet

            allChunks.update(chunkSet)
        # analyze allChunks for contributions
        for chunk in allChunks:
            print(chunk.stop-chunk.start)
        
        
        
        
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

    def _pairRealGeneContr(self, IDa, IDb=None):
        """Calculate realized gene contribution of one individual to another.

        Parameters
        ----------
        IDa : string/integer
          Identifier of parent individual
        IDb : string/integer
          Identifier of progeny individual

        Return
        ------
        out : float
          Realized gene contribution of IDa to IDb.

        """
        if not self._pairIsAncestor(IDa, IDb):
            return 0.0  # Caters also for the case IDb = None.
        ga = self[IDa].generation
        gb = self[IDb].generation
        contr = self._realGeneContr(ID=IDb, maxGenBack=gb - ga)
#        try:
#            return(contr[IDa])
#        except KeyError:
#            return 0.0
        return contr.get(IDa, 0.0)  # If not found, contribution is 0.0.


    def realGeneContr(self, ID=None, generation=None,
                      IDb=None, generationB=None):
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

          Note: This function repeatedly calls 'pairRealGeneContr', which is
          highly inefficient. For a more efficient version, directly work
          with the function 'realGeneContr' or devise a new algorithm.

        """
        if ID is not None and isSingle(ID) and isSingle(IDb):
            return self._pairRealGeneContr(IDa=ID, IDb=IDb)
        return matrixComputer(idGetter=self.getID,
                              function=self._pairRealGeneContr,
                              ID=ID, generation=generation,
                              IDb=IDb, generationB=generationB,
                              symmetric=False)





    def recombinationNumber(self, ID=None, generation=None):
        """ Return the total number of recombinations that occured, starting
        from the founder genomes. """
        ID = IDGenerationHandler(idGetter=self.getID, ID=ID, generation=generation)
        ret = dict()
        for i in ID:
            ind = self[i]
            ct = 0
            for g in ind.gametes:
                for chrom in g.skeleton:
                    ct += len(chrom['lineages'])
            ret[i] = ct
        return ret

    def getLineage(self, pos, chrom, ID=None, generation=None):
        if not (0.0 <= pos <= self.chromLengths[chrom]):
            raise ValueError("'pos' must be a location within the chromosome.")
        ID = IDGenerationHandler(idGetter=self.getID, ID=ID, generation=generation)
        ret = dict()
        lin = namedtuple('lin', ['paternal', 'maternal'])
        for i in ID:
            ind = self[i]
            temp = list()
            for g in ind.gametes:
                sk = g.skeleton[chrom]
                temp.append(sk['lineages'][bisect_left(a=sk['locations'], x=pos)])

            ret[i] = lin(paternal=temp[0], maternal=temp[1])
        return ret

    def getSegmentBorders(self, ID=None, generation=None):
        ID = IDGenerationHandler(idGetter=self.getID, ID=ID, generation=generation)
        ret = list()   
        for cx in range(len(self.chromLengths)):
            tmp = set()
            for i in ID:
                for g in self[i].gametes:
                   tmp.update(g.skeleton[cx]['locations'])
                   tmp.discard(self.chromLengths[cx])
            ret.append(sorted(list(tmp)))
        return ret
                
        





class genoPop(skPop):


    @classmethod
    def fromSkPop(cls, population, genome, founderFlesh):

        fP = cls(chromLengths=population.chromLengths,
                 crossoverSimulator=population.crossoverSimulator)
        fP.genome = genome
        for ind in population.values():
            ind.population = fP
        fP.update(population)

        fleshDict = dict()
        for fid in fP.getFounderID():
            try:
                founderFlesh.get(fid)
            except KeyError:
                raise ValueError("'founderFlesh' must be a dictionary with "
                                 "genotypes of all founder individuals.")
            else:
#                for i in range(2):
#                    fP[fid].gametes[i].flesh = founderFlesh[fid][i]
                for gam, flesh in zip(fP[fid].gametes, founderFlesh[fid]):
                    gam.flesh = flesh
                    fleshDict.update({gam.lineage: flesh})
        fP.fleshDict = fleshDict
        return fP


#    def _addFleshGamete(self, gamete):
#        """Add Genotypes to a single gamete"""
#        numChrom = self.genome.numChrom
#        flesh = [None] * numChrom
#        for ichrom in range(numChrom):
#            positions = self.genome.chromosomes[ichrom].positions
#            lin = gamete.skeleton[ichrom]['lineages']
#            loc = gamete.skeleton[ichrom]['locations']
#            tmp = [None] * len(positions)
#            ixl = 0
#            for lc, li in zip(loc, lin):
#                ixr = bisect.bisect(positions, lc)
#                # find founder
#                tmp[ixl:ixr] = self.founderGametes[li].flesh[ichrom][ixl:ixr]
#                ixl = ixr
#            flesh[ichrom] = tmp
#        gamete.flesh = flesh

    def addFlesh(self, ID=None, generation=None):
        """Add genotypes to individuals

        Parameters
        ----------
        generation : integer (default = all generations)
          The generation for which genotypes should be computed.

        """
        ID = IDGenerationHandler(self.getID, ID, generation)
        for id_ in ID:
            ind = self[id_]
            ind.gametes[0].enflesh()
            ind.gametes[1].enflesh()


    def genoIndividuals(self, ID=None, generation=None):
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
        ID = IDGenerationHandler(idGetter=self.getID, ID=ID,
                                 generation=generation)
        ret = OrderedDict()
        for i in ID:
            ind = self[i]
            ret[i] = np.vstack((np.hstack(ind.gametes[0].flesh),
                                np.hstack(ind.gametes[1].flesh)))
        return ret

    def genoMatrix(self, ID=None, generation=None):
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
        dic = self.genoIndividuals(ID=ID, generation=generation)
        ret = namedtuple('ret', ('matrix','rownames', 'colnames'))

        if self.genome.locusNames is None:
            colnames = None
        else:
            colnames = tuple(chain(*self.genome.locusNames))

        return ret(matrix=np.vstack(dic.values()),
                   rownames=tuple(dic.keys()),
                   colnames=colnames)

    def GRM(self, ID=None, generation=None, method="vanRaden1", p=None,
            assureBiallelic=False):
        mA = ('vanRaden1', 'vanRaden2', 'SMC', 'euclidean')
        if method not in mA:
            raise ValueError("'method' must be one of {0}.".format(mA))
#        ID = IDGenerationHandler(idGetter=self.getID, ID=ID,
#                                 generation=generation)
        X = self.genoMatrix(ID=ID, generation=generation)
        ID = X.rownames
        X = X.matrix

        if method == 'SMC':
            if not np.all(np.apply_along_axis(lambda x: len(set(x)), 0, X)):
                raise NotImplementedError("SMC is currently only implemented "
                                          "for biallelic loci.")
            GRM = SMC_routine(X)

        if method == 'euclidean':
            Z = X[::2,:] + X[1::2,:]
            GRM = cdist(XA=Z, XB=Z, metric='euclidean')

        if method in ('vanRaden1', 'vanRaden2'):
            if not assureBiallelic:
                if not np.all(np.apply_along_axis(lambda x: len(set(x)), 0, X)):
                    raise ValueError("Not all loci are biallelic.")

            if p is None:
                p = X.mean(axis=0)

            Z = X[::2,:] + X[1::2,:]
            m = Z.shape[0]

            # Drop monomorphic loci.
            poly = np.abs(p - 0.5) < 0.5
            p = p[poly]
            Z = Z[:, poly]
            n = Z.shape[1]

            Z = Z - 2.0 * p   # center
            if method == 'vanRaden2':
                Z /= np.sqrt(n * 2.0 * p * (1.0 - p))
            elif method == 'vanRaden1':
                Z /= np.sqrt(2.0 * np.sum(p * (1.0 - p)))

            GRM = np.dot(Z, Z.T)

        rtuple = collections.namedtuple('rtuple', ['matrix', 'ID'])
        return rtuple(matrix=GRM, ID=ID)

        return GRM

    def LPS(self, distMin, distMax, ID=None, IDb=None, generation=None,
            generationB=None, method='correlation'):

        mA = ('correlation', 'cosine', 'direction')
        if method not in mA:
            raise ValueError("'method' must be one of {0}.".format(mA))

        ID, IDb =  dualIDGenerationHandler(idGetter=self.getID, ID=ID,
                                           generation=generation, IDb=IDb,
                                           generationB=generationB)
        positions = self.genome.positions
        lstA = list()
        lstB = list()
        for cx, pos in enumerate(positions):
            chromA = np.vstack(g.flesh[cx] for i in ID for g in self[i].gametes)
            chromB = np.vstack(g.flesh[cx] for i in IDb for g in self[i].gametes)

            corA, corB = LPS_routine(chromA=chromA, chromB=chromB, pos=pos,
                                     distMin=distMin, distMax=distMax)
            #  %timeit -n1 -r1
            lstA.append(corA)
            lstB.append(corB)

        corA = np.hstack(lstA)
        corB = np.hstack(lstB)
        sub = np.logical_and(np.logical_not(np.isnan(corA)),
                             np.logical_not(np.isnan(corB)))
        corA = corA[sub]
        corB = corB[sub]

        if method == 'correlation':
            return 1.0 - correlation(corA, corB)
        elif method == 'cosine':
            return 1.0 - cosine(corA, corB)
        else:
            return np.mean(np.sign(corA) == np.sign(corB))

    # TODO: implement generation handling
    def LD(self, ID, distMin, distMax, generation=None,
           metric='r2'):
        mA = ('r2', 'Dprime', 'Dcor2')
        if metric not in mA:
            raise ValueError("'method' must be one of {0}.".format(mA))

        ID = IDGenerationHandler(idGetter=self.getID, ID=ID,
                                 generation=generation)
        positions = self.genome.positions
        lst = list()
        for cx, pos in enumerate(positions):
            chrom = np.vstack(g.flesh[cx] for i in ID for g in self[i].gametes)
            if metric == 'r2':
                ret = LD_r2_routine(chrom=chrom, pos=pos, distMin=distMin,
                                    distMax=distMax)
                ret = np.square(ret)

            elif metric == 'Dcor2':
                p = chrom.mean(0)
                ret = LD_routine(chrom=chrom, p=p, pos=pos, distMin=distMin,
                                   distMax=distMax, metric=2)
                ret = np.square(ret)

            elif metric == 'Dprime':
                p = chrom.mean(0)
                ret = LD_routine(chrom=chrom, p=p, pos=pos, distMin=distMin,
                                   distMax=distMax, metric=3)
            lst.append(ret)
        ret = np.hstack(lst)
        return np.nanmean(ret)

    # TOD0: Revise, test, implement D, Dprime
    def pairwiseLD(self, ID, generation=None, metric='r2', asDataFrame=True):
        ID = IDGenerationHandler(idGetter=self.getID, ID=ID,
                                 generation=generation)
        mA = ('r','r2', 'D', 'Dprime', 'Dcor', 'Dcor2')
        if metric not in mA:
            raise ValueError("'method' must be one of {0}.".format(mA))

        ID = IDGenerationHandler(idGetter=self.getID, ID=ID,
                                 generation=generation)
        positions = self.genome.positions
        lst = list()
        for cx, pos in enumerate(positions):
            chrom = np.vstack(g.flesh[cx] for i in ID for g in self[i].gametes)
            if metric in ('r', 'r2'):
                ret = np.corrcoef(chrom, rowvar=0)
                ret = ret[np.triu_indices(n=ret.shape[1], k=1)]
                if metric == 'r2':
                    ret = np.square(ret)
            elif metric == 'D':
                p = chrom.mean(0)
                ret = LD_routine(chrom=chrom, p=p, pos=pos, distMin=0.0,
                                   distMax=pos[-1] - pos[0], metric=1)
            elif metric in ('Dcor', 'Dcor2'):
                p = chrom.mean(0)
                ret = LD_routine(chrom=chrom, p=p, pos=pos, distMin=0.0,
                                   distMax=pos[-1] - pos[0], metric=2)
                if metric == 'Dcor2':
                    ret = np.square(ret)

            elif metric == 'Dprime':
                p = chrom.mean(0)
                ret = LD_routine(chrom=chrom, p=p, pos=pos, distMin=0.0,
                                   distMax=pos[-1] - pos[0], metric=3)

            N = ret.size
            out = np.recarray(N, dtype=[('LocusA', object), ('LocusB', object),
                                        ('chrom', int),
                                        ('dist', float), ('LD', float)])

            if hasattr(self.genome, 'locusNames'):
                temp = self.genome.locusNames[cx]
            else:
                temp = range(len(pos))

            for i, (la, lb), d, LD in zip(count(),
                                          combinations(temp, r=2),
                                          combinations(pos, r=2),
                                          ret):

                out[i, ] = (la, lb, cx + 1, d[1] - d[0], LD)

            lst.append(out)

        out = np.hstack(lst)
        if asDataFrame:
            return pd.DataFrame(out)
        return out

#class Population(Mapping):

#    """Population
#
#    Instance variables
#    ------------------
#    name : population name
#
#    founderGametes : gametes that founded the population
#
#    individuals : individuals of the population
#
#    founderIndividuals: individuals that founded the population
#
#    genome : Genome object
#
#    pedigree : Pedigree object
#
#    crossoverSimulator : CrossoverSimulator object
#
#    Instance methods
#    ----------------
#
#    """
#
#    def __init__(self, crossoverSimulator, genome, ID=None):
#        """Instantiate a new Population
#
#        Parameters
#        ----------
#        crossoverSimulator : CrossoverSimulator
#          Used to simulate crossover events.
#
#        genome: Genome
#          A Genome that has at least chromosome lengths specified.
#
#        name : integer or string
#          Name or identifier of the population
#
#
#        """
#        self._ID = ID
#        self.crossoverSimulator = crossoverSimulator
#        self.genome = genome
#        self._pedigree = Pedigree(population=self)  # Empty pedigree.
#        self._individuals = collections.OrderedDict()
#        self._founderIndividuals = collections.OrderedDict()
#        self._gametes = collections.OrderedDict()
#        self._founderGametes = collections.OrderedDict()

#    # PROPERTIES
#
#    # ID
#    @property
#    def ID(self):
#        return self._ID
#
#    # crossoverSimulator
#    @property
#    def crossoverSimulator(self):
#        return self.__crossoverSimulator
#
#    @crossoverSimulator.setter
#    def crossoverSimulator(self, value):
#        if not isinstance(value, CrossoverSimulator):
#            raise ValueError("A valid CrossoverSimulator must be provided.")
#        self.__crossoverSimulator = value
#
#    # genome
#    @property
#    def genome(self):
#        return self._genome
#
#    @genome.setter
#    def genome(self, value):
#        if not isinstance(value, Genome):
#            raise ValueError("A valid Genome must be provided.")
#
#        if not value.chromLengths:
#            raise ValueError("The genome must at minimum have chromLengths.")
#
#        self._genome = value
#
#    # pedigree
#    @property
#    def pedigree(self):
#        return self._pedigree
#
##    # individuals
##    @property
##    def individuals(self):
##        return self._individuals
##
##    # founderIndividuals
##    @property
##    def founderIndividuals(self):
##        return self._founderIndividuals
#
#    # gametes
#    @property
#    def gametes(self):
#        return self._gametes
#
#    # founderGametes
#    @property
#    def founderGametes(self):
#        return self._founderGametes
#
#    def __getitem__(self, item):
#        return self.pedigree[item]
#
#    def __len__(self, item):
#        return self.pedigree.__len__
#
#    def __iter__(self):
#        return self.pedigree.items().__iter__()

#    def checkInTune(self):
#        """ Checks if the population and its pedigree are in tune. """
#        IDpop = list(self.individuals.keys())
#        IDped = list(self.pedigree.keys())
#        if not set(IDpop) == set(IDped):
#            raise RuntimeError("Not in tune: Sets of IDs are unequal.")
#        for i in IDped:
#            vped = self.pedigree[i]
#            fIDped, mIDped = vped.fatherID, vped.motherID
#            vpop = self.individuals[i]
#            fIDpop, mIDpop = vpop.father.ID, vpop.mother.ID
#            if not fIDped == fIDpop or not mIDped == mIDpop:
#                raise RuntimeError("Not in tune: Different parents occured.")

#    def _add(self, ind):
#
#        if not isinstance(ind, Individual):
#            raise ValueError("Must be of class 'Individual'.")
#
#        if ind.ID in self.individuals:
#            raise ValueError("An individual with the same ID is "
#                             "already present.")
#
#        fatherGamete = ind.gametes[0]
#        motherGamete = ind.gametes[1]
#        for gamete in (fatherGamete, motherGamete):
#            if gamete.ID in self.gametes:
#                raise ValueError("A gamete with the same ID is "
#                                 "already present.")
#
#        isFatherFounder = isinstance(fatherGamete, FounderGamete)
#        isMotherFounder = isinstance(motherGamete, FounderGamete)
#        if isFatherFounder or isMotherFounder:
#            self.founderIndividuals[ind.ID] = ind
#
#        if isFatherFounder:
#            self.founderGametes[fatherGamete.ID] = fatherGamete
#
#        if isMotherFounder:
#            self.founderGametes[motherGamete.ID] = motherGamete
#
#        self.individuals[ind.ID] = ind
#        self.gametes[fatherGamete.ID] = fatherGamete
#        self.gametes[motherGamete.ID] = motherGamete

    # NEW INDIVIDUALS

#    def createFounder(self, ID, F, DH):
#        self.pedigree.createFounder(ID=ID, DH=DH, F=F)
#
#    def createFounders(self, n, ID=None, F=None, DH=None):
#        """ Create founders to the population and the pedigree
#
#        Parameters
#        ----------
#        nFounder : integer
#          The number of founder individuals that should be added
#
#        ID : list (default = None)
#          Unique identifiers for the founders (integer or strings)
#
#        F : list (default = None)
#          Inbreeding coefficients of founders.
#
#        Note: By adding Founders, the population and its pedigree are
#        overwritten. You cannot add founders except at the very beginning of
#        the population.
#
#        """
##        self._pedigree = Pedigree()
##        self._individuals.clear()
##        self._founderIndividuals.clear()
##        self._gametes.clear()
##        self._founderGametes.clear()
#        for i in range(n):
#            self.createFounder(ID=ID[i] if ID is not None else getNewID(),
#                               F=F[i] if F is not None else 0.0,
#                               DH=DH[i] if DH is not None else False)
#
#

#
#    def cross(self, ID, fatherID, motherID):
#        """Cross two individuals
#
#        Parameters
#        ----------
#        father : Individual
#          The father of the cross
#        mother : Individual
#          The mother of the cross
#        crossoverSimulator : CrossoverSimulator
#          An object of type 'CrossoverSimulator'
#        ID : string/integer (default = None)
#          A new identifier for the child. If 'None', a new UUID is obtained.
#
#        Return
#        ------
#        out : Individual
#          The child of the cross.
#
#        """
#        self.pedigree.cross(ID, fatherID, motherID, generation)
#
#    def selfing(self, ID, parentID, generation):
#        """Self-fertilize an individual."""
#        self.pedigree.selfing(ID, parentID, parentID, generation)
#
#    def DH(self, ID, parentID, generation):
#        """Produce a doubled haploid."""
#        self.pedigree.DH(ID, parentID, parentID, generation, True)
#
#    def randomMating(self, size, allowSelfing=True, ID=None, returnID=True):
#        """Perform random mating
#
#        Parameters
#        ----------
#        size : integer
#          Number of progeny to be produced.
#        allowSelfing : boolean
#          Should self-fertilization be allowed?
#
#        """
#        self.pedigree.randomMating(size, size, allowSelfing, ID, returnID)
#
#    def inbreeding(self, size):
#        """Perform self-fertilization
#
#        Parameters
#        ----------
#        size : integer
#
#        If size is a single integer, a number of 'size' selfing progeny are
#        produced from each individual. If it is a list, a variable number of
#        progeny is produced as specified in the list.
#
#        """
#        self.pedigree.selfing(size)
#
#    def roundRobin(self, size):
#        self.pedigree.roundRobin(size)

#    def synthetic(self, size, reciprocals):
#        self.pedigree.synthetic(size, reciprocals)
#
#    def allID(self):
#        """ Return all IDs currently present in the population """
#        return list(self.individuals.keys())
#
#    def getID(self, generation):
#        """ Return IDs of individuals in a given generation
#
#        Parameters
#        ----------
#        generation : integer
#          The generation for which IDs should be returned.
#
#        Return
#        ------
#        out : list
#          A list with IDs.
#
#        """
#        return self.pedigree.getID(generation)
#
#    def getInd(self, ID):
#        """ Get an Individual. """
#        return self[ID]




#    def _createIndividual(self, ID=None, fatherID=None, motherID=None, DH=None):
#        founder = fatherID is None and motherID is None
#        gametes = [None, None]
#        if founder:
#
##            populationFounderIndividual(pedigree=self.pedigree, population, DH, F, alleles, ID=None)
#
#            if DH:  # Only one gamete is produced.
#                gametes = [FounderGamete(chromLengths=self.genome.chromLengths,
#                                         lineage=getNewID(), ID=getNewID())] * 2
#            else:
#                for i in range(2):
#                    gametes[i] = FounderGamete(
#                                     chromLengths=self.genome.chromLengths,
#                                     lineage=getNewID(),
#                                     ID=getNewID())
#        else:
#            father, mother = self.getInd(fatherID), self.getInd(motherID)
#
#
#        self._add(Individual(ID=ID, father=father, mother=mother,
#                             gametes=gametes))

#    @classmethod
#    def fromPedigree(cls, crossoverSimulator, genome,
#                     ID, fatherID, motherID, DH=None, F=None):
#        population = cls(crossoverSimulator, genome, ID)
#        population.pedigree.fromPedigree(ID, fatherID, motherID, DH, F)
#        return population


    #@classmethod
#    def fromPedigree(cls, ID, fatherID, motherID, generation=None):

#        # Assistant functions
#        def isFounder(ID, fatherID, motherID):
#            return fatherID is None and motherID is None
#
#        def checkIDNone(ID):
#            """ Check if any ID is None. """
#            if any(i is None for i in ID):
#                raise ValueError("IDs cannot be 'None'.")
#
#        def checkIDUnique(ID):
#            """ Check if all IDs are unique. """
#            if not allUnique(ID):
#                raise ValueError("Not all IDs are unique")
#
#        def checkSameLength(*args):
#            """ Check if all given arguments have the same length. """
#            if not len(set(map(len, args))) == 1:
#                raise ValueError("Arguments do not have the same length.")
#
#        def checkOwnAncestor(ID, fatherID, motherID):
#            """Check if all Individuals are not their own ancestor."""
#            if any((i == f or i == m for i, f, m in
#                    zip(ID, fatherID, motherID))):
#                raise ValueError("An individual has itself as ancestor.")
#
#        def checkDHSelfing(DH, fatherID, motherID):
#            """Check if all DH are selfing progeny."""
#            if any(not f == m for d, f, m in
#                   zip(DH, fatherID, motherID) if d):
#                raise ValueError("Not all doubled haploids are "
#                                 "selfing progeny.")
#        def checkF():
#            """Check if inbreeding coefficients of all founders are provided."""
#            founderID = [i for i in self.ID if self.isFounder(i)]
#            diff = set(founderID) - self.F.keys()
#            if diff:
#                raise ValueError(
#                        "If 'F' is provided, inbreeding "
#                        "coefficients for all founder individuals must "
#                        "be specified. They are missing for the following "
#                        "individuals:\n\n " + '{}'.format(diff))
#
#       # Construct Pedgiree
#        checkIDNone(ID); checkIDUnique(ID)
#        checkSameLength(ID, fatherID, motherID)
#        checkOwnAncestor(ID, fatherID, motherID)
#        if DH is not None:
#            checkSameLength(ID, DH)
#            checkDHSelfing(DH, fatherID, motherID)
#
#        # Fill Pedigree
#        for ix, i, f, m in zip(itertools.count(), ID, fatherID, motherID):
#            try:
#                tmp = F.get(i, 0.0)
#            except (AttributeError, KeyError):
#                tmp = 0.0
#
#            self.createEntry(ID=i, father=f, mother=m, generation=None,
#                             DH=False if DH is None else DH[ix], F=tmp)






