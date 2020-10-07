#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 11:20:06 2020

@author: tdraper
"""

import numpy as np
from math import factorial
import matplotlib.pyplot as plt
from PRNG import rCarry, linearCongruential, fibonacci, s24r55

#%%

class chisq:
    def __init__(self, prng, binsize = 100, normalize = False):
        # The binsize
        self.binsize  = binsize                    # number of bins
        # Random numbers to evaluate
        self.size     = len(prng.random)           # number of random numbers
        self.random   = prng.random/prng.modulus if normalize else prng.random
        # For the chi2-test
        self.expected = len(prng.random)/binsize*np.ones(binsize) # sums to size
        # The bin edges and observed values
        self.observed, self.bins = np.histogram(a     = self.random, 
                                                bins  = self.binsize, 
                                                range = (0,1))    
        
    # How to print the class  
    def __str__(self):
        return r"chi2-test: {}".format(self.test())
        
    # Perform the chi2 test          
    def test(self):
        res = np.divide((self.observed - self.expected)**2, 
                         self.expected)
        return np.sum(res)
    
    # Show a histogram of the observed data
    def showBins(self):
        plt.hist(x        = self.random, 
                 bins     = self.bins,
                 density  = True,
                 histtype = 'step',
                 alpha = 0.5)
        plt.plot(np.arange(0, 1, 1/self.binsize),
                 np.ones(self.binsize),
                 alpha = 0.5)
        plt.show()

#%%

class permutationTest:
    def __init__(self, prng, tuplesize = 3):
        self.random      = prng.random
        self.tuplesize   = tuplesize
        self.tuplelength = len(prng.random)//tuplesize
        self.tuples      = self._createTuples()
        self.bins        = self._fillBins()
    
    def __str__(self):
        return r"chi2-test: {}".format(self.chi2())
    
    def _createTuples(self):
        tuples  = np.zeros((self.tuplelength, self.tuplesize))
        stop    = self.tuplesize*self.tuplelength
        step    = self.tuplesize
        for start in range(self.tuplesize):
            tuples[:,start] = self.random[start:stop:step]
        return tuples
    
    def _tupleToBin(self, ntuple, n, binnumber = 0):
        assert len(ntuple)-1 == n
        if n!=0:
            index         = np.argmax(ntuple)
            newbinnumber  = binnumber + index*factorial(n)
            newcoordinate = np.delete(ntuple,index)
            return self._tupleToBin(newcoordinate, n-1, newbinnumber)
        return binnumber
    
    def _fillBins(self):
        bins = np.zeros(factorial(self.tuplesize))
        for i in range(self.tuplelength):
            binnumber  = self._tupleToBin(self.tuples[i], 
                                          self.tuplesize-1)
            bins[binnumber]+=1
        return bins
            
    def showBins(self):
        plt.hist(x          = np.arange(factorial(self.tuplesize)),
                 weights    = self.bins,
                 bins       = factorial(self.tuplesize),
                 histtype   = 'bar',
                 density    = True,
                 rwidth     = 0.9,
                 color      = 'orange',
                 alpha      = 0.5)
        plt.plot()
    
    def chisq(self):
        observed = self.bins
        expected = self.tuplelength/factorial(self.tuplesize)*np.ones(factorial(self.tuplesize))
        res = np.divide((observed-expected)**2,
                        expected)
        return np.sum(res)
        
        
        
        
#%% Initiate RCARRY
        
modulus   = 2**20
# Random numbers of O(modulus) of # length 24.
seeds     = np.array([999999, 123456, 535322, 874372, 1000000, 875767, 232001, 600201, 43421, 688765, 213813, 564123, 
                      34237, 142342, 123321, 613549, 97876, 927496, 827248, 248382, 324294, 952349, 23424, 993494])                                       
s         = 10
carry_bit = 0
size      = 10**5
# Generate random numbers for RCARRY

rcarry = rCarry(seeds, modulus, size, s, carry_bit)

#%% Initiate Fibonacci

#%% Initiate s24r55

#%% Chi^2 test RCARRY

chisq_rcarry = chisq(rcarry,
                     binsize = 100,
                     normalize = True)
chisq_rcarry.showBins()
chisq_rcarry.__str__()

#%% Permutation test RCARRY

pt_rcarry = permutationTest(rcarry)
pt_rcarry.showBins()
pt_rcarry.chisq()

#%% Permutation test Fibonacci

#%% Permutation test s24r55























