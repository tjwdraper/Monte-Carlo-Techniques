#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 11:20:06 2020

@author: tdraper
"""

from numpy.random import randint
from PRNG import rCarry, fibonacci, s24r55
from Test import chisq, permutationTest
       

#%% Initiate RCARRY
        

modulus   = 2**20
# Random seed of O(modulus) of length 24.
seed_rc   = randint(0, high = modulus-1, size = 24) 
# And define s, the initial carry bit and number size.                                   
s         = 10
carry_bit = 0
size      = 10**5

# Generate random numbers for RCARRY.
rcarry = rCarry(seed      = seed_rc, 
                modulus   = modulus, 
                size      = size, 
                s         = s, 
                carry_bit = carry_bit)


#%% Initiate Fibonacci


modulus = 10**6
# Random seed of O(modulus) of length 2.
seed_f  = randint(0, high = modulus-1, size = 2)
size    = 10**5

# Generate random numbers from Fibonacci
fibo = fibonacci(seed    = seed_f,
                 modulus = modulus,
                 size    = size)



#%% Initiate s24r55


modulus = 10**6
# Random seed of O(modulus) of length 55.
seed_sr = randint(0, high = modulus-1, size = 55)
size    = 10**5

# Generate random numbers for s24r55.
sr = s24r55(seed    = seed_sr,
            modulus = modulus,
            size    = size)


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


pt_fibo = permutationTest(fibo)
pt_fibo.showBins()
pt_fibo.chisq()


#%% Permutation test s24r55


pt_s24r55 = permutationTest(sr)
pt_s24r55.showBins()
pt_s24r55.chisq()




















