#!/usr/bin/python

from time import sleep



@profile
def its_time_for_the_calculator(foo):
    """ It's time for the calculator. """
    if not isinstance(foo, int):
        return None

    a = []
    for i in range(foo):
        a.append(i)
    
    
    def so_slow(bar):
        """ Simulate a slow function. """
        sleep(5)
        return bar
    
    

    b = so_slow(a)

    c = 0
    for i in range(foo):
        c += i

    return None

#def main():
#    print(its_time_for_the_calculator(100000))
#
#if __name__ == "__main__":
#    main()




import math
import scipy.optimize
import numpy as np
import random

class CrossoverSimulator:

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

    @profile
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
            # about twice times faster for small nXO (1-5):
            return sorted([L * random.random() for _ in range(nXO)])
#            return sorted(np.random.uniform(0.0, L, nXO).tolist())


        lambda1 = Lstar/50.0 * (m + 1) * (1.0 - p)
        lambda2 = Lstar/50.0 * p

        while True:
            # chiasma and intermediate points
            nPoints = np.random.poisson(lambda1)

            # which point is the first chiasma?
            first = random.randrange(m + 1)
#            first = np.random.choice(range(m + 1))
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
            if random.random() < 0.5: # flip coin -> chiasma
                nXO += 1
                XOLocations.append(chiLocations[i])

        return XOLocations


if __name__ == "__main__":
    cS = CrossoverSimulator(2, 0.5, False)
    cS.simulateCrossover(100.0)