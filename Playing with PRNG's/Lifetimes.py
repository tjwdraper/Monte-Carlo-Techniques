#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 10:48:58 2020

Creation of the Lifetimes class, which can be used to analyse the lifetimes of 
several PRNG's. Methods include storing the lifetimes of an algorithms from 
different starting seeds, calculating the average lifetime, and showing the
lifetimes in a histogram.

@author: Tom Draper (s4468201)
"""
import numpy as np
from PRNG import MiddleSquare, LogisticMap, LinearCongruential
import matplotlib.pyplot as plt

#%% Define functions to calculate the average lifetime and show a histogram
#   of the lifetimes.

class Lifetimes:
    def __init__(self, size, method):
        self.lifetimes = np.zeros(size)
        self.weights   = np.arange(start = 1,
                                   stop = size + 1)
        self.size      = size
        self.method    = method
    
    def calcLifetimes(self, show = False):
        for seed in range(1, self.size):
            # Initiate a PRNG.
            if self.method == 'Middle Square':
                algorithm = MiddleSquare(seed, self.size)
            elif self.method == 'Logistic Map':
                algorithm = LogisticMap(seed, self.size)
            elif self.method == 'Linear Congruential':
                algorithm = LinearCongruential(seed, self.size)
            else:
                print('Pick a valid method (logistic map, middle square or linear congruential)')
                break         
            # Calculate the lifetime.
            algorithm.generateNumbers()
            # Increment the corresponding element of lifetimes.
            self.lifetimes[algorithm.lifetime] += 1
            # If true, show the lifetime of this PRNG with seed.
            if show:
                print('Seed: {} \t Lifetime: {}'.format(seed, algorithm.lifetime))
    
    def avgLifetime(self):
        return np.dot(self.lifetimes, self.weights) / np.sum(self.lifetimes)    
    
    # Truncate lifetimes, omit final zeros.
    def truncate(self):
        final_index = np.nonzero(self.lifetimes)[0][-1] + 1
        self.lifetimes = self.lifetimes[:final_index]
    
    # Show of histogram of the distribution of lifetimes.
    def showLifetimes(self):
        plt.hist(np.arange(len(self.lifetimes)), # An array of possible lifetimes (1,2,...).
                 bins = 1000,
                 weights = self.lifetimes,       # How many times lifetime n has appeared.
                 density = True,
                 histtype = 'step')
        plt.xlabel('Lifetime')
        plt.ylabel('Density')
        plt.title(self.method)
        plt.show()