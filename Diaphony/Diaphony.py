#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 11:29:38 2020

@author: tdraper
"""

import numpy as np
from numpy.random import randint
from math import pi
from PRNG import middleSquare, logisticMap, linearCongruential, rCarry


#%%

class pointSet:
    def __init__(self, dim=1, size = 1000, modulus = 10**5, method = 'rCarry'):
        methods = ['MiddleSquare', 'LogisticMap', 'LinearCongruential', 'rCarry']
        assert method in methods
        self.dim    = dim
        self.random = self._initiateRandomNumbers(dim, size, modulus, method)
    
    def _pickRandomizer(self, size, modulus, method):
        picker = {'MiddleSquare'       : middleSquare(seed    = randint(0, high = modulus-1),
                                                      modulus = modulus,
                                                      size    = size,),
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
                                                carry_bit = 0)}
        alg = picker[method]
        return alg.random
   
    def _initiateRandomNumbers(self, dim, size, modulus, method):
        random = np.zeros((dim, size))
        for idx, _ in enumerate(random):
            random[idx] = self._pickRandomizer(size, modulus, method)
        return random.T
                


#%%


class Diaphony:
    def __init__(self, dim, size):
        self.dim      = dim
        self.size     = size
        self.strength
        
    def _initStrength(self):
        pass
    
    def beta(self, x):
        pass


#%%
        

class Euler1D(Diaphony):
    def __init__(self, size):
        super().__init__(1, size)
        self.strength = self._initStrength()

    def _initStrength(self):
        strength = np.arange(start = 1, stop = self.size + 1, step = 1)**2
        strength = np.divide(np.ones(self.size),
                             strength)
        return 3/pi**2*strength
    
    def beta(x):
        return 1-6*np.abs(x)*(1-np.abs(x))


class Euler(Diaphony):
    def __init__(self, dim, size):
        super().__init__(dim, size)
        self.strength = self._initStrength()

    def _initStrength(self):
        strength = np.zeros([self.size for i in range(self.dim)])
        for vec, _ in np.ndenumerate(strength): # the index of 'strength' is the vector n
            res = np.ones(self.dim)
            for idx, comp in enumerate(vec):
                if comp != 0:
                    res[idx] = 3/(pi*comp)**2
            strength[vec] = np.prod(res)/(2**self.dim-1)
        return strength
    
    def beta(self, x):
        res = np.ones(self.dim)
        for i, xmu in enumerate(x):
            res[i] += Euler1D.beta(xmu)
        return (np.prod(res)-1)/(2**self.dim-1)


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
            
dia_eul = Euler(dim = 4, size = 10)
print(dia_eul.strength)
print('sum of strengths: {}'.format(np.sum(dia_eul.strength)))
        
#%%

dia_gul = Gulliver(dim = 5, size = 2, s = 0.5)
print(dia_gul.strength)        
print('sum of strengths: {}'.format(np.sum(dia_gul.strength))) 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        