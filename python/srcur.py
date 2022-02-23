from lusol import lusol
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

class srcur:
    def __init__(self,A, rank):
        self.A = A
        (m,n) = A.shape
        self.m = m
        self.n = n

        trlu = lusol(A, rank)
        rank = trlu.rank
        self.rank = rank
        self.ap = trlu.p-1
        self.aq = trlu.q-1

        L = trlu.getL0()
        L = L[np.ix_(self.ap,self.ap[:rank])]
        U = trlu.getU()
        U = U[np.ix_(self.ap[:rank],self.aq)]

        L21 = L[np.ix_(range(rank,m),range(rank))]
        U12 = U[np.ix_(range(rank),range(rank,n))]

        self.A11 = A[np.ix_(self.ap[:rank],self.aq[:rank])]
        self.A12 = A[np.ix_(self.ap[:rank],self.aq[rank:])]
        self.A21 = A[np.ix_(self.ap[rank:],self.aq[:rank])]
        self.A22 = A[np.ix_(self.ap[rank:],self.aq[rank:])]
        self.O = np.random.rand(20,m-rank)

        OL = csc_matrix.dot(self.O,L21)
        OLU = csc_matrix.dot(OL,U12)
        self.S = csc_matrix.dot(self.O,self.A22) - OLU
        
        self.corelu = lusol(self.A11, rank)

    def refactor(self):
        rank = self.rank
        self.corelu = lusol(A, rank)

        L11 = self.corelu.getL0()
        U11 = self.corelu.getU()

        U12 = spsolve(L11,self.A12)
        L21 = spsolve(U11.transpose(),self.A21.transpose())

        OL = csc_matrix.dot(self.O,L21)
        OLU = csc_matrix.dot(OL,U12)
        self.S = csc_matrix.dot(self.O,self.A22) - OLU

    def maxS(self):
        s_c = np.argmax(np.sum(self.S**2,0))
        nrank = self.rank
        u = A[self.ap[nrank:],self.aq[nrank+s_c]].todense()
        u = np.ndarray.flatten(np.array(u))
        v1 = self.A12[:,s_c].todense()
        v2 = self.corelu.solveAv(v1)
        v3 = csc_matrix.dot(self.A21,v2);
        max_col = u - v3
        s_r = np.argmax(np.abs(max_col))
        alpha = max_col[s_r]
        return (alpha,s_r,s_c)
    
    def  maxA11tinv(self,alpha,s_r,s_c):
        nrank = self.rank
        Omega = np.random.rand(20,nrank+1)
        a21 = np.ndarray.flatten(np.array(self.A21[s_r,:].transpose().todense()))
        a12 = np.ndarray.flatten(np.array(self.A12[:,s_c].todense()))

        w21 = np.ndarray.flatten(self.corelu.solveAtv(a21))
        w12 = np.ndarray.flatten(self.corelu.solveAv(a12))
        Omega[:,nrank] = Omega[:,nrank] - Omega[:,:nrank]@w21
        Omega[:,:nrank] = self.corelu.solveAm(Omega[:,:nrank].transpose()).transpose()
        Omega[:,nrank] = Omega[:,nrank]/alpha
        Omega[:,:nrank] = Omega[:,:nrank]- np.outer(Omega[:,nrank],w12);

        a_c = np.argmax(np.sum(Omega**2,0))
        e = np.zeros(nrank+1)
        e[a_c] = 1

        e[nrank] = e[nrank] - np.dot(w12,e[:nrank])
        e[:nrank] = self.corelu.solveAtv(e[:nrank])
        e[nrank] = e[nrank]/alpha
        e[:nrank] = e[:nrank] - w21*e[nrank]

        a_r = np.argmax(np.abs(e))
        beta = e[a_r]
        return (beta, s_r, s_c)
A = csc_matrix(np.random.rand(10,10))
mycur = srcur(A,4)
(alpha,s_r,s_c) = mycur.maxS()
mycur.maxA11tinv(alpha,s_r,s_c)
