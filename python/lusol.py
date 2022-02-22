#from ctypes import c_longlong, c_double, cdll, byref, c_void_p, pointer
from ctypes import *
import numpy as np
from scipy.sparse import csc_matrix 
import os.path
class lusol:
    libclusol = 0
    @classmethod
    def loadlibrary(cls):
        libdir = os.path.dirname(__file__)
        libpath = os.path.abspath(os.path.join(libdir,'libclusol.so'))
        if cls.libclusol == 0:
            cls.libclusol = cdll.LoadLibrary(libpath)
        return
    
    def __init__(self, A, rank ):
        self.rank = rank
        self.A = A
        self.m = A.shape[0]
        self.n = A.shape[1]
        # LUSOL input parameters
        self.maxcol = 5
        self.pivot = 'TRP'
        self.keepLU = 1
        self.Ltol1 = 10
        self.Ltol2 = 10
        self.small = np.finfo(float).eps**.8
        self.Utol1 = np.finfo(float).eps**67
        self.Utol2 = 0
        self.Uspace = 3.0
        self.dens1 = 0.3
        self.dens2 = 0.5

        self.loadlibrary()
        self.set_options()
        self.allocate_and_copy(A)
        self.factorize()
    def size(self):
        return (self.m, self.n)
    
    def set_options(self):
        A = self.A
        m = A.shape[0]
        n = A.shape[1]
        self.luparm = np.zeros(30,dtype=np.int64)
        self.parmlu = np.zeros(30,dtype=np.float64)
        self.luparm[1] = -1
        self.luparm[2] = self.maxcol

        if self.pivot == 'TRP':
            self.luparm[5] = 1
        elif self.pivot == 'TCP':
            self.luparm[5] = 2
        self.luparm[5] = 2
        self.luparm[7] = self.keepLU

        self.parmlu[0] = self.Ltol1
        self.parmlu[1] = self.Ltol2
        self.parmlu[2] = self.small
        self.parmlu[3] = self.Utol1
        self.parmlu[4] = self.Utol2
        self.parmlu[5] = self.Uspace
        self.parmlu[6] = self.dens1
        self.parmlu[7] = self.dens2
        self.luparm_ptr = c_void_p(self.luparm.ctypes.data)
        self.parmlu_ptr = c_void_p(self.parmlu.ctypes.data)

    def allocate_and_copy(self, A):
        m = A.shape[0]
        n = A.shape[1]
        self.m = m
        self.n = n
        # set storage sizes
        self.nelem = csc_matrix.count_nonzero(A)
        self.lena = max(2*self.nelem,10*self.m,10*self.n,10000)

        self.a = np.zeros(self.lena, dtype=np.float64)
        self.indc = np.zeros(self.lena,dtype=np.int64)
        self.indr = np.zeros(self.lena,dtype=np.int64)

        # extract data from A for use in LUSOL
        (indc_tmp, indr_tmp) = np.nonzero(A)
        a_tmp = A[indc_tmp,indr_tmp]
        self.a[0:self.nelem] = a_tmp
        self.indc[0:self.nelem] = indc_tmp+1
        self.indr[0:self.nelem] = indr_tmp+1

        # vectors of length m
        self.p = np.zeros(m,dtype=np.int64)
        self.lenr = np.zeros(m,dtype=np.int64)
        self.locr = np.zeros(m,dtype=np.int64)
        
        self.iqloc = np.zeros(m,dtype=np.int64)
        self.ipinv = np.zeros(m,dtype=np.int64)

        # vectors of length n
        self.q = np.zeros(n,np.int64)
        self.lenc = np.zeros(n,np.int64)
        self.locc = np.zeros(n,np.int64)

        self.iploc = np.zeros(n,np.int64)
        self.iqinv = np.zeros(n,np.int64)
        

        # Integer Scalars
        self.m_ptr = pointer(c_longlong(self.m))
        self.n_ptr = pointer(c_longlong(self.n))
        self.nelem_ptr = pointer(c_longlong(self.nelem))
        self.lena_ptr = pointer(c_longlong(self.lena))
        self.rank_ptr = pointer(c_longlong(self.rank))

        # Vectors of length alen
        self.a_ptr = c_void_p(self.a.ctypes.data)
        self.indc_ptr = c_void_p(self.indc.ctypes.data)
        self.indr_ptr = c_void_p(self.indr.ctypes.data)

        # Vectors of length m 
        self.p_ptr = c_void_p(self.p.ctypes.data)
        self.lenr_ptr = c_void_p(self.lenr.ctypes.data)
        self.locr_ptr = c_void_p(self.locr.ctypes.data)
        self.iqloc_ptr = c_void_p(self.iqloc.ctypes.data)
        self.ipinv_ptr = c_void_p(self.ipinv.ctypes.data)

        # Vectors of length n
        self.q_ptr = c_void_p(self.q.ctypes.data)
        self.lenc_ptr = c_void_p(self.lenc.ctypes.data)
        self.locc_ptr = c_void_p(self.locc.ctypes.data)
        self.iploc_ptr = c_void_p(self.iploc.ctypes.data)
        self.iqinv_ptr = c_void_p(self.iqinv.ctypes.data)
        


        
    def factorize(self):
        w = np.zeros(self.n,dtype=c_double)
        w_ptr = c_void_p(w.ctypes.data)
        
        inform = c_longlong(0)
        inform_ptr = pointer(inform)
        
        self.libclusol.clu1fac(
            self.m_ptr,
            self.n_ptr,
            self.nelem_ptr,
            self.lena_ptr,
            self.rank_ptr,
            self.luparm_ptr,
            self.parmlu_ptr,
            self.a_ptr,
            self.indc_ptr,
            self.indr_ptr,
            self.p_ptr,
            self.q_ptr,
            self.lenc_ptr,
            self.lenr_ptr,
            self.locc_ptr,
            self.locr_ptr,
            self.iploc_ptr,
            self.iqloc_ptr,
            self.ipinv_ptr,
            self.iqinv_ptr,
            w_ptr,
            inform_ptr)
        # Correct the rank 
        self.rank = self.luparm[15]
        return

    def getL0(self):
        m = self.m
        n = self.n
        lenL0 = self.luparm[20]
        lena = self.lena

        li = np.zeros(lenL0+m,dtype=np.int64)
        lj = np.zeros(lenL0+m,dtype=np.int64)
        la = np.zeros(lenL0+m)

        li[0:lenL0] = self.indc[lena-lenL0:lena]-1
        lj[0:lenL0] = self.indr[lena-lenL0:lena]-1
        la[0:lenL0] = -self.a[lena-lenL0:lena]

        li[lenL0:] = range(m)
        lj[lenL0:] = range(m)
        la[lenL0:] = np.ones(m)
        L = csc_matrix((la, (li, lj)))
        return L

    def getU(self):
        m = self.m
        n = self.n 
        lenU = np.sum(self.lenr)

        ui = np.zeros(lenU,dtype=np.int64)
        uj = np.zeros(lenU,dtype=np.int64)
        ua = np.zeros(lenU)

        k1 = 0;
        k2 = 0;
        p = self.p
        for k in range(self.rank):
            i = p[k]-1
            rlen = self.lenr[i]
            rloc = self.locr[i]-1

            k2 = k1+rlen

            ui[k1:k2] = i*np.ones(rlen)
            uj[k1:k2] = self.indr[rloc:rloc+rlen]-1
            ua[k1:k2] = self.a[rloc:rloc+rlen]

            k1 = k1 + rlen
            
        U = csc_matrix((ua, (ui,uj)),shape = (m,n))
        return U


    def clu6sol(self, b, mode = 5):
        # Modes
        # 1 x solves L x = b
        # 2 x solves L'X = b
        # 3 x solves U x = b
        # 4 x solves U'x = b
        # 5 x solves A x = b
        # 6 x solves A'x = b
        lenb = len(b)
        m,n = self.size()
        v = np.zeros(m,dtype=np.float64)
        w = np.zeros(n,dtype=np.float64)
        if mode == 1:
            if lenb != m:
                raise(ValueError('b has incorrect size'))
            v = b
        elif mode == 2:
            v = b
        elif mode == 3:
            v = b
        elif mode == 4:
            w = b
        elif mode == 5:
            v = b
        elif mode == 6:
            w = b
        else:
            raise(ValueError('Unrecognized Mode'))

        v_ptr = c_void_p(v.ctypes.data)
        w_ptr = c_void_p(w.ctypes.data)
        mode_ptr = pointer(c_longlong(mode))
        ret_inform_ptr = pointer(c_longlong(0))

        self.libclusol.clu6sol(
            mode_ptr,
            self.m_ptr,
            self.n_ptr,
            v_ptr,
            w_ptr,
            self.lena_ptr,
            self.luparm_ptr,
            self.parmlu_ptr,
            self.a_ptr,
            self.indc_ptr,
            self.indr_ptr,
            self.p_ptr,
            self.q_ptr,
            self.lenc_ptr,
            self.lenr_ptr,
            self.locc_ptr,
            self.locr_ptr,
            ret_inform_ptr)

        if mode in [1,2,3]:
            x = v
        else:
            x = w

        inform = self.luparm[9]
        resid = self.parmlu[19]
        return (x,inform, resid)

    def solveA(self,b):
        return self.clu6sol(b,5)

    def solveAt(self,b):
        return self.clu6osl(b,6)
    


