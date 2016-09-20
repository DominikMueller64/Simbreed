# math
import numpy as np
import math
import scipy.optimize  # brentq
## random numbers
import random

from numba import jit

## %%
#class CrossoverSimulator:
#
#    """Simulate crossover locations on a single meiotic product using
#    the Stahl model.
#
#    Instance Variables
#    ------------------
#    m : Interference parameter
#
#    p : Proportion of chiasmata from non-interference mechanism
#
#    obligateChiasma : Boolean indicating if obligate chiasma is required
#
#    Chiasma locations are a superposition of two processes: a proportion p
#    exhibiting no interference, and a proportion (1-p) following the
#    chi-square model with interference parameter m. Crossover locations are
#    derived by thinning the chiasma locations with probability 1/2.
#
#    Simulations are under the Stahl model with the interference parameter
#    being an integer. This is an extension of the chi-square model, but with
#    chiasmata being the superposition of two processes, one following the
#    chi-square model and the other exhibiting no interference.
#
#    The source code is a translation from
#    `R/simcross <https://github.com/kbroman/simcross/>`_ package
#
#    References
#    ----------
#
#    Copenhaver, G. P., Housworth, E. A. and Stahl, F. W. (2002) Crossover
#    interference in arabidopsis. Genetics 160, 1631–1639.
#
#    Foss, E., Lande, R., Stahl, F. W. and Steinberg, C. M. (1993) Chiasma
#    interference as a function of genetic distance. Genetics 133, 681–691.
#
#    Zhao, H., Speed, T. P. and McPeek, M. S. (1995) Statistical analysis of
#    crossover interference using the chi-square model. Genetics 139, 1045–1056.
#
#
#    """
#
#    def __init__(self, m=10, p=0, obligateChiasma=False):
#        """ Instantiate a CrossoverSimulator
#
#        Parameters
#        ----------
#        m : nonnegative integer (default = 10)
#          Interference parameter (m = 0 corresponds to no interference)
#        p : float between 0 and 1 (default = 0)
#          Proportion of chiasmata from no-interference mechanism
#          (p = 0 gives a pure chi-square model)
#        obligateChiasma : boolean
#          If TRUE, require an obligate chiasma on the 4-strand bundle
#          at meiosis.
#
#        """
#        self.m = m
#        self.p = p
#        self.obligateChiasma = obligateChiasma
#
#        @property
#        def m(self):
#            return self.m
#
#        @property
#        def p(self):
#            return self.p
#
#        @property
#        def obligateChiasma(self):
#            return self.obligateChiasma
#
#    def __calculateLStar(self, L):
#        m = self.m; p = self.p
#
#        if L <= 50.0:
#            raise ValueError("Chromosome must have length > 50 cM")
#        if m < 0 or not isinstance(m, int):
#            raise ValueError("Parameter 'm' must be a non-negative integer.")
#
#        if p < 0 or p > 1:
#            raise ValueError("Parameter 'p' must be in [0, 1]")
#        if p == 1:
#            m = p = 0  # If p == 1, might take m = 0, p = 0
#
#        def funcToZero(Lstar, L, m, p):
#            def dpoisson(k, lambda_):
#                """Density function of the Poisson distribution.
#        
#                The Poisson distribution is the limit of the binomial distribution
#                for large N.
#        
#                Parameters
#                ----------
#                k : non-negative integer
#                  The non-negative integer for which the discrete density function
#                  should be evaluated
#                lambda_ : positive real number
#                  The expectation/variance.
#        
#                Returns
#                -------
#                out : float
#                  Value of the Possion density function.
#        
#                """
#                return lambda_**k * math.exp(-lambda_) / math.factorial(k)
#            if m == 0:
#                denom = 1 - math.exp(-Lstar / 50.0)
#            else:
#                lambda1 = Lstar/50.0 * (m + 1) * (1 - p)
#                lambda2 = Lstar/50.0 * p
#                sm = 0.0
#                for i in range(m + 1):
#                    sm += dpoisson(i, lambda1) * (m + 1 - i) / (m + 1)
#                denom = 1 - sm * math.exp(-lambda2)
#            return 2 * L - 2 * Lstar / denom
#
#        return scipy.optimize.brentq(f=funcToZero, a=math.exp(-8),
#                                     b=L, args=(L, m, p))
#
#    def simulateCrossover(self, L):
#        """Simulate crossing over events
#
#        Parameters
#        -----------
#        L : positive float
#          Length of the chromosome in centiMorgan
#
#        Return
#        ------
#        out : list
#          A list with the positions of crossover events
#
#        """
#        m = self.m; p = self.p
#        obligateChiasma = self.obligateChiasma
#        if obligateChiasma:
#            Lstar = self.__calculateLStar(L)
#        else:
#            Lstar = L
#        # no-interference model is a lot easier
#        if m == 0:
#            if obligateChiasma:
#                # rejection sampling to get at least one chiasma
#                while True:
#                    nXO = np.random.poisson(Lstar / 50.0)
#                    if nXO > 0:
#                        break
#                nXO = np.random.binomial(nXO, 0.5)
#            else:
#                nXO = np.random.poisson(L / 100.0)
#            # about twice times faster for small nXO (1-5):
#            return sorted([L * random.random() for _ in range(nXO)])
##            return sorted(np.random.uniform(0.0, L, nXO).tolist())
#        lambda1 = Lstar / 50.0 * (m + 1) * (1.0 - p)  # adjust expectation
#        lambda2 = Lstar / 50.0 * p
#        while True:
#            # chiasma and intermediate points
#            nPoints = np.random.poisson(lambda1)
#            # which point is the first chiasma?
#            first = random.randrange(m + 1)  # much faster
##            first = np.random.choice(range(m + 1))
#            if first > nPoints:
#                nIchi = 0
#            else:
#                nIchi = nPoints // (m + 1) + int(first < (nPoints % (m + 1)))
#                # be careful to use explicit integer division (//) here
#            # no. chiasma from no interference process
#            if p > 0:
#                nNIchi = np.random.poisson(lambda2)
#            else:
#                nNIchi = 0
#            if not obligateChiasma or nIchi + nNIchi > 0:
#                break
#        # locations of chiasmata and intermediate points for process with
#        # interference
#        pointLocations = sorted([L * random.random() for _ in range(nPoints)])
#        # move every (m+1)st point back to front
#        nChi = 0
#        for j in range(first, nPoints - 1, m + 1):
#            # Here, j > nChi is required, otherwise a mess will happen
#            pointLocations[nChi] = pointLocations[j]
#            nChi += 1
#
#        # chiasma locations from non-interference process
#        NIchiLocations = [L * random.random() for _ in range(nNIchi)]
#        # combine interference and no interference chiasma locations
#        chiLocations = sorted(pointLocations[:nChi] + NIchiLocations)
#        # thin by 1/2
#        nXO = 0
#        XOLocations = list()
#        for i in range(len(chiLocations)):
#            if random.random() < 0.5:  # flip coin -> chiasma
#                nXO += 1
#                XOLocations.append(chiLocations[i])
#
#        return XOLocations
# %%

#import numpy as np
#from numba import jit
#import random
#import math
#from numba import jit

class CrossoverSimulator:

    """Simulate crossover locations on a single meiotic product using
    the Stahl model.

    Instance Variables
    ------------------
    m : Interference parameter

    p : Proportion of chiasmata from non-interference mechanism

    obligateChiasma : Boolean indicating if obligate chiasma is required
    
    This is a Python translation of the crossover simulator used in the
    R-package 'simcross' of Karl. W. Broman [Add link to github here].

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
        if self.obligateChiasma:
            Lstar = _calculateLStar(L=L, m=self.m, p=self.p)
        else:
            Lstar = L
        return _simulateCrossover(L=L, m=self.m, p=self.p,
                                  obligateChiasma=self.obligateChiasma,
                                  Lstar=Lstar)



@jit(nopython=True)
def factorial(k):
    fk = 1.0
    for i in range(1, k + 1):
        fk *= i
    return fk
    
#%timeit factorial(5)
#%timeit math.factorial(5)
#%timeit math.gamma(6)

@jit(nopython=True)
def dpoisson(k, lambd):
    return lambd**k * math.exp(-lambd) / factorial(k)
#    return math.exp(k * math.log(lambd) - lambd) / factorial(k)  # maybe more stable?

#%timeit dpoisson(1, 4)
#%timeit scipy.stats.poisson.pmf(1,4)

@jit(nopython=True)
def funcToZero(Lstar, L, m, p):

    if m == 0:
        denom = 1.0 - math.exp(-Lstar / 50.0)
    else:
        lambda1 = Lstar / 50.0 * (m + 1) * (1.0 - p)
        lambda2 = Lstar / 50.0 * p
        sm = 0.0
        for i in range(m + 1):
            sm += dpoisson(i, lambda1) * (m + 1 - i) / (m + 1)
        denom = 1.0 - sm * math.exp(-lambda2)
    return 2.0 * L - 2.0 * Lstar / denom


#%timeit funcToZero(Lstar=200, L=300, m=10, p=0.5)

 
def _calculateLStar(L, m, p):
    if L <= 50.0:
        raise ValueError("Chromosome must have length > 50 cM")
    if m < 0 or not isinstance(m, int):
        raise ValueError("Parameter 'm' must be a non-negative integer.")
    if p < 0 or p > 1:
        raise ValueError("Parameter 'p' must be in [0, 1]")
    if p == 1:
        m = p = 0  # If p == 1, might take m = 0, p = 0

    return scipy.optimize.brentq(f=funcToZero, a=math.exp(-8.0),
                                 b=L, args=(L, m, p))



   
@jit(nopython=True)
#@jit
def _simulateCrossover(L, m, p, obligateChiasma, Lstar): 
    # no-interference model is a lot easier
    if m == 0:
        if obligateChiasma:
            # rejection sampling to get at least one chiasma
            while True:
                nXO = np.random.poisson(Lstar / 50.0)
                if nXO > 0:
                    break
            nXO = np.random.binomial(nXO, 0.5)
        else:
            nXO = np.random.poisson(L / 100.0)
        # about twice times faster for small nXO (1-5):
#            return sorted([L * random.random() for _ in range(nXO)])
        temp = np.random.uniform(0.0, L, nXO)
        temp.sort()
        return temp

    lambda1 = Lstar / 50.0 * (m + 1) * (1.0 - p)  # adjust expectation
    lambda2 = Lstar / 50.0 * p
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
            # be careful to use explicit integer division (//) here
        # no. chiasma from no interference process
        if p > 0:
            nNIchi = np.random.poisson(lambda2)
        else:
            nNIchi = 0
        if not obligateChiasma or (nIchi + nNIchi) > 0:
            break
    # locations of chiasmata and intermediate points for process with
    # interference
#        pointLocations = sorted([L * random.random() for _ in range(nPoints)])
    pointLocations = np.random.uniform(0.0, L, nPoints)
    pointLocations.sort()
    # move every (m+1)st point back to front
    nChi = 0
    for j in range(first, nPoints - 1, m + 1):
        # Here, j > nChi is required, otherwise a mess will happen
        pointLocations[nChi] = pointLocations[j]
        nChi += 1

    # chiasma locations from non-interference process
#        NIchiLocations = [L * random.random() for _ in range(nNIchi)]
    NIchiLocations = np.random.uniform(0.0, L, nNIchi)
    # combine interference and no interference chiasma locations
#        chiLocations = sorted(pointLocations[:nChi] + NIchiLocations)
    chiLocations = np.empty(nChi + nNIchi)
    ct = 0
    for i in range(chiLocations.size):
        if i < nChi:
            chiLocations[i] = pointLocations[i]
        else:
            chiLocations[nChi + ct] = NIchiLocations[ct]
            ct += 1
    
#    chiLocations = np.hstack((pointLocations[:nChi], NIchiLocations))
#    chiLocations.sort()  # XXX: The sort is not necessary I guess, only
    # thin by 1/2        # XOLocations has to be finally sorted!

    nXO = 0
    XOLocations = np.empty_like(chiLocations)
    # TODO: Check out numpy.random.choice() here.
    for i in range(chiLocations.size):
        if random.random() < 0.5:  # flip coin -> chiasma
            XOLocations[nXO] = chiLocations[i]
            nXO += 1
    ret = XOLocations[:nXO]  
    
    ret.sort()
    return ret


#%timeit -n10000 -r10 _simulateCrossover(L=100, m=-0.0, p=0.0, obligateChiasma=False, Lstar=100)
#
#cS = CrossoverSimulator(m=0, p=0, obligateChiasma=False)
#%timeit -n1000 -r10 cS.simulateCrossover(L=300.0) # simcross: 3.4 micro
#
#cS = CrossoverSimulator(m=10, p=0, obligateChiasma=False)
#%timeit -n1000 -r10 cS.simulateCrossover(L=300.0)  # simcross: 6.2 micro
#
#cS = CrossoverSimulator(m=10, p=0.5, obligateChiasma=True)
#%timeit -n1000 -r10 cS.simulateCrossover(L=300.0)  # simcross: 91 micro
#
#
#np.hstack((pointLocations[:nChi], NIchiLocations))
#np.stack((pointLocations[:nChi], NIchiLocations), 2)