# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 17:10:28 2020

Implementation of the PRNG's. Some are simple PRNG's, which use only one 
argument to generate new random numbers. In recent years, these algorithms have 
become absolete and have been replaced by shift register PRNG's. I will 
implement one the the first srPRNG's, the rcarry algorithm.

I will create one (abstract) parent class, from which all the algorithms
will be derived.

@author: Tom Draper (s4468201)
"""

import numpy as np
from abc import abstractmethod, ABCMeta
import matplotlib.pyplot as plt

#%% Create a general class for Pseudo Random Number Generators

class prng(metaclass = ABCMeta):
    def __init__(self, seed, modulus):
        self.seed         = seed
        self.register     = seed
        self.random       = np.zeros(modulus) # The list of random numbers
        self.modulus      = modulus           # Size of the domain
        self.lifetime     = 0                 # Lifetime of algorithm
        self.no_generated = 0
    
    def __str__(self):
        print('List of random numbers: \n {}'.format(self.random))
    
    # How to generate new numbers depends on the method. Leave it abstract in
    # the parent class
    @abstractmethod
    def randomGen(self):
        pass
    
    # Update the register, random list and no_generated. This is different for
    # the a shift register and a normal PRNG.
    def update(self, new_num):
        if isinstance(self.register, int):
            self.register = new_num
        elif isinstance(self.register, np.ndarray):
            self.register     = np.roll(self.register, -1)
            self.register[-1] = new_num
        else:
            print('something wrong, check the registers')
    
    # Check if the latest random number (stored in the register) appears in the
    # list of random numbers.
    def appeared(self):
        # For now, only test the normal PRNG's. Shift registers have in general 
        # such a long lifetime that it it not necessary (and implementation is 
        # a bit tedious)
        if isinstance(self.register, int):
            return self.register in self.random[:self.no_generated]
                    
    # Calculate the lifetime of this PRNG.
    def generateNumbers(self):
        while self.no_generated < self.modulus:
            # Generate a new random number
            newElement = self.randomGen()
            # Check if it already appeared in our list of random numbers
            if self.appeared():
                # Truncate the list of random numbers (omit the final 0's)
                self.random   = self.random[:self.no_generated]
                try:
                    self.lifetime = self.no_generated
                except AttributeError:
                    print('The number generator has reached its lifetime, but'
                          + 'the class has no attribute for it to be stored.')
                break
            else:
                # Add new random number to the list
                self.random[self.no_generated] = newElement
                self.no_generated += 1
 
    # Plots pairs (x_{2n+1}, x_{2n}) in a scatter plot
    def showPairs(self, scatterhist = 'scatter'):
        x = self.random[1:self.no_generated:2]/self.modulus
        y = self.random[0:self.no_generated:2]/self.modulus
        # If the list of random numbers is uneven, omit the last entry.
        if len(x) < len(y):
            y = y[:-1]
        # Generate the plot, labels, title etc.
        if scatterhist == 'scatter':
            plt.scatter(x,y)
        elif scatterhist == 'histogram':
            plt.hist2d(x,y, 
                       bins = (100, 100),
                       cmap = plt.cm.jet)
            plt.colorbar()
        else:
            print('Please enter a valid method (scatter or histogram)')
        plt.xlabel('$x_{2n+1}$')
        plt.ylabel('$x_{2n}$')
        plt.title('Distribution of $(x_{2n+1}, x_{2n})$')
        plt.show()
    
#%% Implement several examples, including the midsquare method, the logistic 
#   map and the linear congruential method

# The Middle-Square method    
class MiddleSquare(prng):
    def __init__(self, seed, modulus = 10**4):
        super().__init__(seed, modulus)
    
    def randomGen(self):    
        rand = int(np.floor(self.register**2 / self.modulus**(0.5)) % self.modulus)
        self.update(rand)
        return rand

# The Logistic Map method         
class LogisticMap(prng):
    def __init__(self, seed, modulus = 10**4, a = 4):
        super().__init__(seed, modulus)
        self.a = a
    
    def randomGen(self):
        x_norm = self.register/self.modulus
        rand   = int(np.round(self.a*x_norm*(1-x_norm),
                              decimals = 4)*self.modulus)
        self.update(rand)
        return rand

# The Linear-Congruential method        
class LinearCongruential(prng):
    def __init__(self, seed, modulus = 10**4, a = 4001, c = 945):
        super().__init__(seed, modulus)
        self.a = a
        self.c = c
    
    def randomGen(self):
        rand = (self.a*self.register + self.c) % self.modulus
        self.update(rand)
        return rand
    
#%% Implement an examples of a shift register PRNG: RCARRY

class RCARRY(prng):
    def __init__(self, seeds, modulus, s, carry_bit):
        super().__init__(seeds, modulus)
        self.r     = len(seeds)
        self.s     = s
        self.carry = carry_bit
     
    def randomGen(self):
        rand = self.register[-self.s] - self.register[-self.r] - self.carry
        if rand < 0:
            rand += self.modulus
            self.carry = 0
        else:
            self.carry = 1
        self.update(rand)
        return rand      
        
