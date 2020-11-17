#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 11:29:38 2020

@author: tdraper
"""

import numpy as np
from numpy.random import randint
from math import pi
from PRNG import middleSquare, logisticMap, linearCongruential, rCarry, fibonacci, s24r55
from cmath import exp

I = complex(0,1)

#%%

class pointSet:
    def __init__(self, dim=1, size = 1000, modulus = 10**5, method = 'rCarry', normalize = True):
        methods = ['MiddleSquare', 'LogisticMap', 'LinearCongruential', 'rCarry', 'Fibonacci', 'S24R55']
        assert method in methods
        self.dim    = dim
        self.size   = size
        self.random = self._initiateRandomNumbers(dim, size, modulus, method, normalize)
    
    def _pickRandomizer(self, size, modulus, method, normalize):
        picker = {'MiddleSquare'       : middleSquare(seed    = randint(0, high = modulus-1),
                                                      modulus = modulus,
                                                      size    = size),
                  'LogisticMap'        : logisticMap(seed    = randint(0, high = modulus-1),
                                                     modulus = modulus,
                                                     size    = size),
                  'LinearCongruential' : linearCongruential(seed    = randint(0, high = modulus-1),
                                                            modulus = modulus,
                                                            size    = size),
                  'rCarry'             : rCarry(seed      = randint(0, high = modulus-1, size = 24),
                                                modulus   = modulus,
                                                size      = size,
                                                s         = 10,
                                                carry_bit = 0),
                  'Fibonacci'          : fibonacci(seed = randint(0, high = modulus-1, size = 2),
                                                   modulus = modulus,
                                                   size = size),
                  'S24R55'             : s24r55(seed = randint(0, high=modulus-1, size = 55),
                                                modulus = modulus,
                                                size = size)}
        alg = picker[method]
        return alg.random/modulus if normalize else alg.random
   
    def _initiateRandomNumbers(self, dim, size, modulus, method, normalize):
        random = np.zeros((dim, size))
        for idx, _ in enumerate(random):
            random[idx] = self._pickRandomizer(size, modulus, method, normalize)
        return random.T
                


#%%


class Diaphony:
    def __init__(self, ps):
        self.pointset = ps
        self.dim      = ps.dim
        self.size     = ps.size
        self.strength = self._initStrength()

    def __str__(self):
        return 'Sum of Strengths: {}'.format(self._sumOfStrengths())
        
    def _initStrength(self):
        pass
    
    def _sumOfStrengths(self):
        return 2**self.dim*np.sum(self.strength)
    
    def betaApprox(self, x):
        pass
    
    def betaExact(x):
        pass

#%%
        

class Euler1D(Diaphony):
    def __init__(self, ps):
        assert ps.dim == 1
        super().__init__(ps)

    def _initStrength(self):
        strength = np.arange(start = 1, stop = self.size + 1, step = 1)**2
        strength = np.divide(np.ones(self.size),
                             strength)
        return 3/pi**2*strength
    
    def betaApprox(self, x):
        res = 0
        for n in range(self.size):
            res += 2*self.strength[n]*np.cos(2*pi*(n+1)*x)
        return res
        
    def betaExact(self, x):
        return 1-6*np.abs(x)*(1-np.abs(x))
    
    def diaphony(self):
        t = self.size
        row, col = np.triu_indices(self.size, k=1)
        for n in range(len(row)):
            t += self.betaExact(self.pointset.random[row[n]] - self.pointset.random[col[n]])
        return t/self.size
            
    
class Euler(Diaphony):
    def __init__(self, ps):
        assert ps.dim > 1
        super().__init__(ps)
        # due to symmetry, only calculate the strength of vectors whose components are all positive

    def _initStrength(self):
        strength = np.zeros([self.size for i in range(self.dim)])
        for vec, _ in np.ndenumerate(strength): # the index of 'strength' is the vector n
            res = np.ones(self.dim)
            for idx, comp in enumerate(vec):
                if comp != 0:
                    res[idx] = 3/(pi*comp)**2
            strength[vec] = np.prod(res)/(2**self.dim-1)
        return strength
    
    def betaApprox(self, x):
        res = 0
        for vec, str_n in np.ndenumerate(self.strength):
            ndotz =  np.matmul(permuteSignMatrix(self.dim)*vec, x)
            res   += str_n*np.sum(2*np.cos(2*pi*ndotz))
        return res            
            
    def betaExact(self, x):
        res = np.ones(self.dim)
        for i, xmu in enumerate(x):
            res[i] += Euler1D.beta(xmu)
        return (np.prod(res)-1)/(2**self.dim-1)

    def diaphony(self):
        pass
   
    
#%%


class Gulliver(Diaphony):
    def __init__(self, dim, size, s = 0.1):
        assert 0<s and s<1
        super().__init__(dim, size)
        self.s        = s
        self.strength = self._initStrength()
        
    def _initStrength(self):
        strength = np.zeros([self.size for i in range(self.dim)])
        prefac   = 1/(((1+self.s)/(1-self.s))**self.dim - 1)
        for vec, _ in np.ndenumerate(strength):
            exp = -1*np.sum(np.absolute(vec))
            strength[vec] = prefac*self.s**exp
        return strength
        
    def _beta(self, x):
        res = np.zeros(self.dim)
        for i, xmu in enumerate(x):
            res[i] = (1-self.s**2)/(1-2*self.s*np.cos(2*pi*xmu)+self.s**2)
        res    = np.prod(res) - 1
        prefac = 1/(((1+self.s)/(1-self.s))**self.dim - 1)
        return prefac*res

#%%

ps = pointSet(dim = 2, size = 1000, modulus = 10**5, method = 'LinearCongruential')
eu = Euler1D(ps)
print(eu.diaphony())

#%%
            
dia_eul = Euler(dim = 4, size = 10)
print(dia_eul.strength)
print('sum of strengths: {}'.format(np.sum(dia_eul.strength)))
        
#%%

dia_gul = Gulliver(dim = 5, size = 2, s = 0.5)
print(dia_gul.strength)        
print('sum of strengths: {}'.format(np.sum(dia_gul.strength))) 
        
#%%

def permuteSignMatrix(dim = 3):
    if dim == 1:
        return np.array([1])
    elif dim == 2:
        return np.array([[1,1],[1,-1]])
    else:
        mat = np.zeros((2**(dim-1), dim))
        mat[:,0] = 1
        mat[:2**(dim-2),1:] = permuteSignMatrix(dim-1)
        mat[2**(dim-2):,1:] = -1*permuteSignMatrix(dim-1)
        return mat

mat = permuteSignMatrix(dim = 5)
print(mat)
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        