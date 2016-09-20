#from abc import ABCMeta, abstractmethod, abstractproperty
import collections  # OrderedDict


class Pedigree(collections.OrderedDict):  # Inherits from OrderedDict
    
    @classmethod
    def fromPedigree(cls, ID, father, mother, generation, DH, F):

        pedigree = cls()
        allID = set().union(*(ID, father, mother)) - {None}

        for i in allID:
            pedigree.createEntry(ID=i, father=None, mother=None)

        for ix, i, f, m in zip(itertools.count(), ID, fatherID, motherID):
            ind = pedigree[i]
            for c, p in zip(('father', 'mother'), (f, m)):
                if p is not None:
                    setattr(ind, c, pedigree[p])

        pedigree.sort()
        return(pedigree)
    
    def __init__(self, population=None):
        super().__init__()
        self.population = population
    
    @property
    def population(self):
        return self._population

    @population.setter
    def population(self, value):
        self._population = value

#    def append(self, value):
#        self.update({value.ID: value})

    def append(self, value):
        """ Append an entry. """
        # If child of a parent, parent gets a child.
        for parent in (value.father, value.mother):
            if parent is not None:
                if parent.ID in self:
                    parent.addChild(value)
        # If parent, gets child.
        for ind in self.values():
            if value in (ind.father, ind.mother):
                value.addChild(ind)
        self.update({value.ID: value})

    def extend(self, values):
        for value in values:
            self.append(value)
    
    
    ## The pedigree should not have methods createFounder and createProgeny??
    def createFounder(self, ID, DH, F):
        self.append(populationFounderIndividual(pedigree=self,
                                population=self.population,
                                DH=DH, F=F, ID=ID))
    
    def createProgeny(self, ID, father, mother, generation, DH):
        self.append(populationIndividual(pedigree=self,
                                       population=self.population,
                                       father=father,
                                       mother=mother,
                                       generation=generation,
                                       DH=DH, ID=ID))
        
#    def createEntry(self, ID, father, mother, generation, DH, F):
#
#        if father is mother is None:
#            ind = populationFounderIndividual(pedigree=self,
#                                        population=self.population,
#                                        DH=DH, F=F, ID=ID)
#        else:
#            ind = populationIndividual(pedigree=self,
#                                       population=self.population,
#                                       father=father,
#                                       mother=mother,
#                                       generation=generation,
#                                       DH=DH, ID=ID)
#        self.append(ind)


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
        ct = itertools.count()
        for e in self.values():
            if getattr(e, 'generation', None) == generation:
                next(ct)
        return next(ct)

    def allID(self):
        return list(self.keys())

    def getID(self, generation=None):
        """ Return IDs of individuals in a given generation

        Parameters
        ----------
        generation : integer
          The generation for which IDs should be returned.

        Return
        ------
        out : list
          A list with IDs.

        """
        if generation is None:
            return self.allID()
        return tuple(i for i, v in self.items() if
                getattr(v, 'generation', None) == generation)


    def returnInd(self, ID):
        if isinstance(ID, (int, str)):
            ID = (ID,)
        return tuple(self[i] for i in ID)


    def tabulateInd(self):
        return tuple(map(self.numInd,
                        sorted(set(e.generation for e in self.values()))))

    def lastGen(self):
        """ Return the number of the last generation. """
        return max(e.generation for e in self.values())

    def numGen(self):
        """ Return the number of generations in the pedigree."""
        return len(set(e.generation for e in self.values))


    def _subset(self, Ind):
        """ For internal use only. """
        nPed = Pedigree()
        for ind in Ind:
            deque = collections.deque((ind,))
            while deque:
                v = deque.popleft()
                nPed.append(v)
                for parent in (v.father, v.mother):
                    if parent is not None:
                        deque.append(parent)
        nPed.sort()
        return nPed


    #TODO: Make this public method, use _subset for this
#    def subset(self, ID=None, generation=None):
#        """ Subset a pedigree
#
#        Parameters
#        ----------
#        ID : string/integer or list
#          Identifiers for individuals
#
#        Return
#        ------
#        out : Pedigree
#          A new pedigree instance containing only individuals in ID and their
#          ancestors.
#
#        """
#        if not self.isSorted():
#            raise ValueError("The pedigree must be sorted.")
#        ID = IDGenerationHandler(self.getID, ID, generation)
#
#        nPed = Pedigree()
#        for i in ID:
#            deque = collections.deque([i])
#            while deque:
#                v = self[deque.popleft()]
#                nPed.append(v)
#                for pID in (v.fatherID, v.motherID):
#                    if pID is not None:
#                        deque.append(pID)
#        nPed.sort()
#        return nPed

    def isPadded(self):
        """Check if the pedigree is extended.

        Return
        ------
        out : boolean
          True, if the pedigree is extended, otherwise, False.

        A pedigree is extended if all parents are themselfs in the pedigree
        with unknown origin.

        """
        for ind in self.values():
            for parent in (e.father, e.mother):
                if not (parent is None or parent.ID in self):
                    return False
        return True

    # TODO: update this function
    def pad(self):
        """ Extend the pedigree.

        Inbreeding of founders is assumed to be 0.0.

        """
        for ind in self.values():
            for parent in (e.father, e.mother):
                if not (parent is None or parent.ID in self):
                    createEntry(ID=parent.ID, father=None, mother=None,
                                generation=0)
                    self.move_to_end(parent.ID, False)

    def isSorted(self):
        """Check if the pedigree is sorted."""
        seen = set([None])
        return not any(v.father not in seen or v.mother not in seen or
                       seen.add(v) for v in self.values())
    # TODO: Check if individual are always put in smalles/greatest
    # possible generation.
    def sort(self, placeFounderFirst=True):
        """Sort the pedigree.

        If no data on generations are available, the algorithm provided in
        `Zhang et al. (2009) <http://www.medwelljournals.com/fulltext/
        ?doi=javaa.2009.177.182>`_ is applied to sort from scratch.

        """
        n = len(self)  # We are sure that all IDs are unique.
        for ind in self.values():  # Set all generation data to zero.
            ind.generation = 0
        nPed = Pedigree()  # A new Pedigree.
        while self:
            shc = self.copy()
            for ind1 in shc.values():
                for ind2 in shc.values():
                    if not ind1 is ind2:
                        if ind1 in (ind2.father, ind2.mother):
                            ind1.generation += 1
                            break
                else:
                     nPed.append(self.pop(ind1.ID))

        maxgen = next(reversed(nPed.values())).generation
        for ind in nPed.values():
            if placeFounderFirst and ind.founder:
                ind.generation = 0
            else:
                ind.generation = abs(ind.generation - maxgen)

        self.update(sorted(nPed.items(), key = lambda x: x[1].generation))


#    def lastGen(self):
#        """ Return the number of the last generation. """
#        return next(reversed(self.values())).generation


#    def sort(self):
#        """ In place sort not possible. """
#        temp = self.copy()  # Shallow copy.
#        self.clear()
#        self.update(sorted(temp.items(), key = lambda x: x[1].generation))


    def _tabularMethod(self, IndA, IndB):
        """ Compute the numberator relationship (IBD-) matrix.

        Parameters
        ----------
        IDa : string/integer or list (default = all individuals)
          Identifiers for the individuals.

        Return
        ------
        out : dictionary
          A dictionary with the IBD matrix and the column/row names.

        """
        # Subset pedigree for minimal computation time
        pedigree = self._subset(Ind=set().union(*(IndA, IndB)))
        Ind = pedigree.values()
        x = dict(zip(Ind, itertools.count()))
        indicesA = tuple(x[i] for i in IndA)
        if IndB is IndA or IndB == IndA:
            indicesB = indicesA
        else:
            indicesB = tuple(x[i] for i in IndB)
        mat = np.identity(len(Ind), dtype=np.float_)
        for i, v in enumerate(Ind):
            f = v.father
            m = v.mother
            # Diagonal
            if v.founder:
                if v.DH:
                    mat[i, i] = 2.0
                else:
                    mat[i, i] = 1.0 + v.F
            elif f is not None and m is not None:
                if v.DH:
                    mat[i,i] = 2.0
                else:
                    mat[i, i] = 1.0 + 0.5 * mat[x[f], x[m]]
            # Off-diagonal
            temp = 0.0
            for p in (f, m):
                if p is not None:
                    temp += 0.5 * mat[x[p], 0:i]  # takes 0.9s per 1000 reps
            mat[i, 0:i] = mat[0:i, i] = temp  # takes 0.5s per 1000 reps

        return mat[np.ix_(indicesA, indicesB)]  # takes 0.1s per 1000 reps

#    mat = numpy.zeros([len(self)] * 2, dtype=np.float_, order = 'F')
#
#    now = time.time()
#    for _ in range(1000):
#
#        for i, v in enumerate(self.values()):
#            fID = v.fatherID
#            mID = v.motherID
#            # Diagonal
##            if v.founder:
##                pass #mat[i, i] = 1.0 + v.F
##            elif fID is not None and mID is not None:
##                pass #mat[i, i] = 1.0 + (1.0 if v.DH else 0.5 * mat[ix[fID], ix[mID]])
#            # Off-diagonal
#            j = range(i)  # working with lower triangle
#            temp = 0.0
#            for pID in (fID, mID):
#                if pID is not None:
#                    temp += 0.5 * mat[ix[pID], 0:i]
#            mat[i, 0:i] = mat[0:i, i] = temp
#        indices = [ix[id_] for id_ in ID]
#        mat[np.ix_(indices, indices)]
#    time.time() - now

    def _LDLMethod(self, IndA, IndB):
        pedigree = self._subset(Ind=set().union(*(IndA, IndB)))
        Ind = pedigree.values()
        x = dict(zip(Ind, itertools.count()))
        n = len(pedigree)
        L = np.identity(n=n, dtype=np.float_)
        for j, vj in enumerate(self.values()):
            for i, vi in enumerate(list(self.values())[j + 1:n], start=j + 1):
                fi = vi.father
                mi = vi.mother
                if fi is not None:
                    if mi is not None:
                        L[i, j] = 0.5 * (L[x[fi], j] + L[x[mi], j])
                    else:
                        L[i, j] = 0.5 * L[x[fi], j]
                else:
                    if mi is not None:
                        L[i, j] = 0.5 * L[x[mi], j]
                    else:
                        L[i, j] = 0

        # Relies on the first individual beeing a founder.
        D = np.zeros(n)
        F = np.zeros(n)
        for j, vj in enumerate(self.values()):
            fj = vj.father
            mj = vj.mother
            if fj is not None:
                if mj is not None:
                    if vj.DH:
                        D[j] = 1.0
                    else:
                        D[j] = 0.5 - 0.25 * (F[x[fj]] + F[x[mj]])
                else:
                    D[j] = 0.75 - 0.25 * F[x[fj]]
            else:
                if mj is not None:
                    D[j] = 0.75 - 0.25 * F[x[mj]]
                else:  # The individual must be a founder!
                    D[j] = 1.0 + vj.F
#                    D[j] = 1 # normal code
            F[j] = np.dot(np.square(L[j,]), D) - 1.0

        rtuple = collections.namedtuple('rtuple', ['L', 'D', 'F', 'Ind'])
        return rtuple(L=L, D=D, F=F, Ind=Ind)
#        return {'L': L, 'D': D, 'F': F, 'names' = tuple(i.ID for i in Ind)}

    def makeA(self, ID=None, generation=None,
              IDb=None, generationB=None, method='tabular'):
        ID, IDb =  dualIDGenerationHandler(self.getID, ID, generation,
                                           IDb, generationB)
        Ind, IndB = map(self.returnInd, (ID, IDb))
        
        if method == 'tabular':
            mat = self._tabularMethod(IndA=Ind, IndB=IndB)
        elif method == 'LDL':
            LDL = self._LDLMethod(IndA=Ind, IndB=IndB)
            L = LDL.L
            mat = (L * LDL.D).dot(L.T)
            x = dict(zip(LDL.Ind, itertools.count()))
            mat = mat[np.ix_(tuple(x[i] for i in Ind), tuple(x[i] for i in IndB))]
        else:
            raise ValueError("'method' must be either 'tabular' "
                             "or 'LDL'.")
        rtuple = collections.namedtuple('rtuple', ['matrix', 'rowID', 'colID'])
        return rtuple(matrix=mat, rowID=ID, colID=IDb)
    
    def makeAinv(self, ID=None, generation=None,
                 IDb=None, generationB=None, method='tabular'):
        ID, IDb =  dualIDGenerationHandler(self.getID, ID, generation,
                                           IDb, generationB)
        Ind, IndB = map(self.returnInd, (ID, IDb))

        if method == 'tabular':
            mat = np.linalg.inv(self._tabularMethod(IndA=Ind, IndB=IndB))
        elif method == 'LDL':
            LDL = self._LDLMethod(IndA=Ind, IndB=IndB)
            Linv = np.linalg.inv(LDL.L)
            mat = (Linv.T * (1.0 / LDL.D)).dot(Linv)
            x = dict(zip(LDL.Ind, itertools.count()))
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
        Ind = self.returnInd(ID)
        if method == 'tabular':
            A = self._tabularMethod(IndA=Ind, IndB=Ind)
            F = (np.diag(A) - 1.0).tolist()
        elif method == 'LDL':
            x = dict(zip(LDL.Ind, itertools.count()))
            F = self._LDLMethod(IndA=Ind, IndB=Ind).F.tolist()
            F = [F[x[i]] for i in ID]
        else:
            raise ValueError("'method' must be either 'tabular' "
                             "or 'LDL'.")
        rtuple = collections.namedtuple('rtuple', ['F', 'ID'])
        return rtuple(F=F, ID=ID)

#    def _dominanceRelationship(self, A, x, IDa, IDb):
#        va = self[IDa]; fa = va.fatherID; ma = va.motherID
#        if va.founder:
#                raise ValueError("Individual %s is a founder an cannot be "
#                                 "included." % (IDa))
#        if IDa is IDb or IDa == IDb:
#            return 0.25 * (A[x[fa], x[fa]] + A[x[ma], x[ma]] +
#                           2.0 * A[x[fa], x[ma]])
#        vb = self[IDb]; fb = vb.fatherID; mb = vb.motherID
#        if vb.founder:
#                raise ValueError("Individual %s is a founder an cannot be "
#                                 "included." % (IDb))
#        return 0.25 * (A[x[fa], x[fb]] + A[x[ma], x[mb]] +
#                       A[x[fa], x[mb]] + A[x[ma], x[fb]])
    # TODO: Add support for founder individuals
    # TODO: update
    def makeD(self, ID=None, generation=None,
              IDb=None, generationB=None, method='tabular'):

        ID = IDGenerationHandler(self.getID, ID, generation)
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

    # TODO: Check if matrix L in A = LDL^T contains this information.
    def pairExpGeneContr(self, IDa, IDb=None):
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
        if IDb is None:
            IDb = IDa
        ga = self[IDa].generation
#        if not ga < self[IDp].generation:
#            raise ValueError("%s is not older than %s" % (IDa, IDp))
        qdict = collections.Counter([IDb])
        contr = 0.0  # Contribution of IDa to IDp.
        c = 0.5  # Weighting factor, halved for each step back.
        while qdict:
            nqdict = collections.Counter()
            for i in qdict.keys():
                vi = self[i]
                if vi.generation > ga:  # Otherwise, parents cannot be target.
                    for p in (vi.fatherID, vi.motherID):
                        if p == IDa:  # Contribute c times number of paths.
                            contr += c * qdict[i]
                        elif p is not None:  # Step back.
                            nqdict.update({p: qdict[i]})
            c *= 0.5
            qdict = nqdict
        return contr

    def expGeneContrMatrix(self, ID=None, generation=None,
                           IDb=None, generationB=None):
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
        return matrixComputer(idGetter = self.getID,
                              function = self.pairExpGeneContr,
                              ID = ID, generation = generation,
                              IDb = IDb, generationB = generationB)


    def isAncestor(self, IDa, IDb=None):
        """Check if an individual as an ancestor of another one

        A depth-first search algorithm is used.

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
        if IDb is None:
            return False
        # Depth-first search algorithm for an arbitrary graph.
        visited, stack = set([None]), [IDb]
        while stack:
            vertex = stack.pop()
            if vertex == IDa:
                return True
            if vertex not in visited:
                visited.add(vertex)
                v = self[vertex]
                stack.extend(set([v.fatherID, v.motherID]) - visited)
        return False
        # old implementation
#        ga = self.getGeneration(IDa)
#        compset = set([IDp])
#        while compset:
#            new_compset = set()
#            for cand in compset:
#                if cand == IDa:
#                    return True
#                gcand = self.getGeneration(cand)
#                if ga < gcand:
#                    dt = self.getDetails(cand)
#                    new_compset.update((dt['father'], dt['mother']))
#                    new_compset.discard(None)
#
#            compset = new_compset
#        return False

    def isRelated(self, ID1, ID2):
        return isAncestor(ID1, ID2) or isAncestor(ID2, ID1)

    def isAncestorMatrix(self, ID=None, generation=None,
                           IDb=None, generationB=None):
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

        Note: This function repeatedly calls 'isAncestor', which is highly
        suboptimal to solve the problem. For an efficient version, the
        depth-first search algorithm in 'isAncestor' could be easily
        minimally modified.

        """
        return matrixComputer(idGetter = self.getID,
                              function = self.isAncestor,
                              ID = ID, generation = generation,
                              IDb = IDb, generationB = generationB)


    def toDataFrame(self):
        """Convert pedigree to a pandas DataFrame

        Return
        ------
        out : DataFrame
          A pandas DataFrame holding the pedigree.

        """
        if not self:
            return pd.DataFrame()

        columns = ('ID', 'father', 'mother', 'generation', 'DH')
        tmp = ([str(getattr(v, a, None)) for a in columns] for v in self.values())
#        for v in self.values():
#            for a in columns:
#                print(str(getattr(v, a, None)))
        return pd.DataFrame(tmp, columns=columns)

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
        temp = self.toDataFrame()
        temp.to_csv(file, sep=sep, na_rep=na_rep, index=False)


    def recode(self):
        # Keys are immutable, so entries must be eventually deleted.
        x = dict(zip(self, itertools.count(start=1)))
        for i in tuple(self):
            ni = x[i]
            ind = self.pop(i)
            ind.ID = ni
            self.update({ni: ind})


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
        return [i for i, v in self.items() if
                getattr(v, 'generation', None) == generation and
                getattr(v, 'DH', None)]

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
        flag = True
        if isinstance(ID, (str, int)):
            ID = [ID]
            flag = False
        if not set(self).issuperset(set(ID)):
            raise ValueError("Not all individuals in ID are contained in "
                             "the pedigree")
        noff = tuple(len(self[i].children) for i in ID)
        return f if flag else f[0]

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
        flag = True
        if isinstance(ID, (str, int)):
            ID = [ID]; flag = False
        if not set(self).issuperset(set(ID)):
            raise ValueError("Not all individuals in ID are contained in "
                             "the pedigree")
        f = [getattr(self[id_], 'founder', None) for id_ in ID]
        return f if flag else f[0]



#    def createEntry(self, ID, father, mother, generation=None, DH=False, F=None):
#        """ Adds an entry to the pedigree and deligates the creation of a
#        new individual if the pedigree is attached to a population. """
#
#        ind = Individual(ID, father, mother)
#        ind.pedigree = self
#        ind.generation = generation
#        ind.DH = DH
#        ind.F = F
#        self.append(ind)
##        if not self.isPadded():
##            self.pad()
##        if not self.isSorted() or generation is None:
##            self.sort()
#        # TODO: create gamete
#        if self.population is not None:
#            pass



#class rawPedigree(Pedigree):
#    # TODO: adapt
#    def append(self, value):
#        self.update({value.ID: value})
#
#    def isPadded(self):
#        """Check if the pedigree is extended.
#
#        Return
#        ------
#        out : boolean
#          True, if the pedigree is extended, otherwise, False.
#
#        A pedigree is extended if all parents are themselfs in the pedigree
#        with unknown origin.
#
#        """
#        for ind in self.values():
#            for parent in (e.father, e.mother):
#                if not (parent is None or parent.ID in self):
#                    return False
#        return True
#
#    # TODO: update this function
#    def pad(self):
#        """ Extend the pedigree.
#
#        Inbreeding of founders is assumed to be 0.0.
#
#        """
#        for ind in self.values():
#            for parent in (e.father, e.mother):
#                if not (parent is None or parent.ID in self):
#                    createEntry(ID=parent.ID, father=None, mother=None,
#                                generation=0)
#                    self.move_to_end(parent.ID, False)
#
#    def isSorted(self):
#        """Check if the pedigree is sorted."""
#        seen = set([None])
#        return not any(v.father not in seen or v.mother not in seen or
#                       seen.add(v) for v in self.values())
#    # TODO: Check if individual are always put in smalles/greatest
#    # possible generation.
#    def sort(self, placeFounderFirst=True):
#        """Sort the pedigree.
#
#        If no data on generations are available, the algorithm provided in
#        `Zhang et al. (2009) <http://www.medwelljournals.com/fulltext/
#        ?doi=javaa.2009.177.182>`_ is applied to sort from scratch.
#
#        """
#        n = len(self)  # We are sure that all IDs are unique.
#        for ind in self.values():  # Set all generation data to zero.
#            ind.generation = 0
#        nPed = rawPedigree()  # A new Pedigree.
#        while self:
#            shc = self.copy()
#            for ind1 in shc.values():
#                for ind2 in shc.values():
#                    if not ind1 is ind2:
#                        if ind1 in (ind2.father, ind2.mother):
#                            ind1.generation += 1
#                            break
#                else:
#                     nPed.append(self.pop(ind1.ID))
#
#        maxgen = next(reversed(nPed.values())).generation
#        for ind in nPed.values():
#            if placeFounderFirst and ind.founder:
#                ind.generation = 0
#            else:
#                ind.generation = abs(ind.generation - maxgen)
#
#        self.update(sorted(nPed.items(), key = lambda x: x[1].generation))
#
#
#    def createEntry(self, ID, father, mother, generation=None):
#        """ Adds an entry to the pedigree and deligates the creation of a
#        new individual if the pedigree is attached to a population. """
#
#        ind = Individual(pedigree=self, ID=ID, father=father,
#                                 mother=mother, generation=generation)
#        self.append(ind)
#
#    @classmethod
#    def fromPedigree(cls, ID, fatherID, motherID):
#
#        pedigree = cls()
#        allID = set().union(*(ID, fatherID, motherID)) - {None}
#
#        for i in allID:
#            pedigree.createEntry(ID=i, father=None, mother=None)
#
#        for ix, i, f, m in zip(itertools.count(), ID, fatherID, motherID):
#            ind = pedigree[i]
#            for c, p in zip(('father', 'mother'), (f, m)):
#                if p is not None:
#                    setattr(ind, c, pedigree[p])
#
#        pedigree.sort()
#        return(pedigree)


