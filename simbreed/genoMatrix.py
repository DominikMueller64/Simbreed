from collections import Iterable, Counter
from itertools import count



class namedMatrix:
    
    def __init__(self, matrix, rownames=None, colnames=None):
        self.matrix=matrix
        self._rownames=dict((str(i), ix) for i, ix in zip(rownames, count()))
        self._colnames=dict((str(i), ix) for i, ix in zip(colnames, count()))
        
        
    @property
    def rownames(self):
        return self._rownames
    
    @rownames.setter
    def rownames(self, value):
        if not nrow(self.matrix) == len(rownames):
            raise ValueError("The number of rownames must correspond "
                             "to the number of rows of the matrix.")
           
    @property
    def colnames(self):
        return self._colnames
    
    @colnames.setter
    def colnames(self, value):
        if not ncol(self.matrix) == len(colnames):
            raise ValueError("The number of colnames must correspond "
                             "to the number of columns of the matrix.")

        
    def __getitem__(self, value):
        try:
            return self.matrix[value]
        except IndexError:
            if not len(value) == 2:
                raise IndexError("Indixes must be separated by a comma.")
            rows, cols = value
            
            if not isinstance(rows, (int, slice)):
                rows = self.getIndices(value=rows, 
                                       ref=self.rownames)
            if not isinstance(cols, (int, slice)):
                cols = self.getIndices(value=cols, 
                                       ref=self.colnames)
                
            return self.matrix[rows, cols]
                
    def getIndices(self, value, ref):
        if isinstance(value, str):
            try:
                ix = ref[value]
            except TypeError:
                raise ValueError("Names are required.")
            except KeyError:
                raise ValueError("{0} is not present. "
                                 "".format(value))
            else:
                return ix
                
        elif isinstance(value, Iterable):
            return tuple(self.getIndices(elem, ref) for elem in value)






class genoMatrix:
    
#    def __init__(self, matrix, rownames=None, colnames=None):
#        self.matrix=matrix
#        self._rownames=dict((str(i), ix) for i, ix in zip(rownames, count()))
#        self._colnames=dict((str(i), ix) for i, ix in zip(colnames, count()))
#        
#        
#    @property
#    def rownames(self):
#        return self._rownames
#    
#    @rownames.setter
#    def rownames(self, value):
#        if not nrow(self.matrix) == len(rownames):
#            raise ValueError("The number of rownames must correspond "
#                             "to the number of rows of the matrix.")
#           
#    @property
#    def colnames(self):
#        return self._colnames
#    
#    @colnames.setter
#    def colnames(self, value):
#        if not ncol(self.matrix) == len(colnames):
#            raise ValueError("The number of colnames must correspond "
#                             "to the number of columns of the matrix.")
#
#        
#    def __getitem__(self, value):
#        try:
#            return self.matrix[value]
#        except IndexError:
#            if not len(value) == 2:
#                raise IndexError("Indixes must be separated by a comma.")
#            rows, cols = value
#            
#            if not isinstance(rows, (int, slice)):
#                rows = self.getIndices(value=rows, 
#                                       ref=self.rownames)
#            if not isinstance(cols, (int, slice)):
#                cols = self.getIndices(value=cols, 
#                                       ref=self.colnames)
#                
#            return self.matrix[rows, cols]
#                
#    def getIndices(self, value, ref):
#        if isinstance(value, str):
#            try:
#                ix = ref[value]
#            except TypeError:
#                raise ValueError("Names are required.")
#            except KeyError:
#                raise ValueError("{0} is not present. "
#                                 "".format(value))
#            else:
#                return ix
#                
#        elif isinstance(value, Iterable):
#            return tuple(self.getIndices(elem, ref) for elem in value)


    def freq(self):
        return {col: dict(Counter(self[:, col])) for col in self.colnames}

    def allBiallelic(self):
        return all(1 <= len(x) <= 2 for x in self.freq().values())

    # FIXME: Is this correct?
    def GRM(self, method='vanRaden2', freq=None,
            assureBiallelic=False):
        
        if not assureBiallelic:
            if not self.allBiallelic():
                raise ValueError("All loci must be biallelic for the "
                                 " genomic relationship matrix.")
        if freq is None:
            freq = self.freq()
        
        mat = self.matrix[::2,:] + self.matrix[1::2,:]
        m = nrow(self.matrix)
        
        p = self.matrix.mean(axis=0)
#        p = np.fromiter((freq[col][1] / m for col in self.colnames),
#                                   dtype=np.float64)
        
        # Drop monomorphic loci.      
        poly = np.abs(p - 0.5) < 0.5
        p = p[poly]
        mat = mat[:, poly]
        n = ncol(mat)
        
        mat = mat - 2.0 * p   # center
        if method == 'vanRaden2':
            mat /= np.sqrt(n * 2.0 * p * (1.0 - p))
        elif method == 'vanRaden1':
            mat /= np.sqrt(2.0 * np.sum(p * (1.0 - p)))

        GRM = np.dot(mat, mat.T)
        return GRM






































