#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 11:20:06 2020

@author: tdraper
"""

import numpy as np
import matplotlib.pyplot as plt
from PRNG import RCARRY, LinearCongruential

#%%

class chi2:
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
    def show(self):
        plt.hist(x        = self.random, 
                 bins     = self.bins,
                 density  = True,
                 histtype = "step",
                 alpha = 0.5)
        plt.plot(np.arange(0, 1, 1/self.binsize),
                 np.ones(self.binsize),
                 alpha = 0.5)
        plt.show()
       
#%% Test RCARRY
        
modulus   = 2**20
seeds     = np.array([999999, 123456, 535322, 874372, 1000000, 875767, 232001, 
                      600201, 43421, 688765, 213813, 564123, 34237, 142342, 
                      123321, 613549, 97876, 927496, 827248, 248382, 324294, 
                      952349, 23424, 993494]) # Random numbers of O(modulus) of
                                              # length 24.
s         = 10
carry_bit = 0

# Generate random numbers for RCARRY

rcarry = RCARRY(seeds, modulus, s, carry_bit)
rcarry.generateNumbers()

# Test and show

test_rcarry = chi2(rcarry, 
                   normalize = True)
print(test_rcarry)
test_rcarry.show()

#%% Test Linear Congruential Generator

lincon = LinearCongruential(seed    = 0, 
                            modulus = 10000,
                            a = 65617,
                            c = 23432)
lincon.generateNumbers()

test_lincon = chi2(lincon, 
                   normalize = True)
print(test_lincon)
test_lincon.show()


















