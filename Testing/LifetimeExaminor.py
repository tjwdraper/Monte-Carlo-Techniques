#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 10:12:12 2020

In this file, we will examine several (bad) examples of pseudo-random number
generators (PRNG). We will examine their lifetime, i.e. how long does can we 
produce a stream of random numbers before we run in a duplicate. The analysis
will show that the average lifetime is way below the expected value.

Secondly, we will look at several shift register PRNG's. These are more popular, 
as the life time increases significantly with the register size, for a 
fixed modulus. We will look how the pairs (x_{2n+1}, x_{2n}) are distributed
to see if the RCARRY and Mersene Twister mimic true randomness.

@author: Tom Draper (s4468201)
"""

import numpy as np
from PRNG import LinearCongruential, RCARRY
from Lifetimes import Lifetimes

#%% Set the parameters. Reproduce the lifetime figure in the lecture notes.
        
size = 10**4              # Size of the domain
method = 'Middle Square'  # Middle Square, Logistic Map or Linear Congruential.      
        
lifetimes = Lifetimes(size, method)
lifetimes.calcLifetimes() 
lifetimes.truncate()
lifetimes.showLifetimes()
  
#%% Show the distribution of points of the Linear Congruential method.

lincon = LinearCongruential(seed    = 0, 
                            modulus = 10000,
                            a = 65617,
                            c = 23432)
lincon.generateNumbers()
lincon.showPairs()

# Completely different behaviour for different moduli, a and c.

#%% Analyse the distribution of the shift register PRNG's

# Initiate parameters, the modulus, the seed (with r=24), s and the carry bit.
modulus   = 2**20
seeds     = np.array([999999, 123456, 535322, 874372, 1000000, 875767, 232001, 
                      600201, 43421, 688765, 213813, 564123, 34237, 142342, 
                      123321, 613549, 97876, 927496, 827248, 248382, 324294, 
                      952349, 23424, 993494]) # Random numbers of O(modulus) of
                                              # length 24.
s         = 10
carry_bit = 0

alg = RCARRY(seeds, modulus, s, carry_bit)
alg.generateNumbers()
alg.showPairs('histogram') # can also do 'scatter' but I like this style more.




