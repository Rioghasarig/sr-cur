#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 12:18:32 2022

@author: oekenta
"""
from ctypes import c_ulonglong, c_double, cdll, byref
import numpy as np

class lusol:
    liblusol = 0
    @classmethod
    def loadlibrary(cls):
        cls.liblusol = cdll.LoadLibrary('/home/grad/oekenta/sr-cur/src/libclusol.so')
    def __init__(self, A : np.array ):
        # LUSOL input parameters
        self.rank = 0
        self.maxcol = 0
        self.pivot = 0
        self.keepLU = 0
        self.Ltol1 = 0
        self.Ltol2 = 0
        self.small = 0
        self.Utol1 = 0
        self.Utol2 = 0
        self.Uspace = 0
        self.dens1 = 0
        self.dens2 = 0
        
        #LU1FAC Inputs
        self.m = c_ulonglong(A.shape[0])
        self.n = c_ulonglong(A.shape[1])
        self.nelem = c_ulonglong(np.count_nonzero(A))
        self.lena = c_ulonglong(10000)
        self.ap = c_ulonglong*A.shape[0]
        self.aq = c_ulonglong*A.shape[1] 

    def factorize():
        
        
A = np.array([[1,2],[3,4]])
l = lusol(A) 
l.loadlibrary()
