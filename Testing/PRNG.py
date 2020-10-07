# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 17:10:28 2020

Implementation of the PRNG's. Some are simple PRNG's, which use only one 
argument to generate new random numbers. In recent years, these algorithms have 
become absolete and have been replaced by shift register PRNG's. I will 
implement simple shift register PRNG's, including the rcarry algorithm.

I will create one (abstract) parent class, from which all the algorithms
will be derived.

@author: Tom Draper (s4468201)
"""

import numpy as np
from abc import abstractmethod, ABCMeta
import matplotlib.pyplot as plt

#%% Create a general class for Pseudo Random Number Generators.

class prng(metaclass = ABCMeta):
    def __init__(self, seed, modulus, size):
        self.seed         = seed
        self.register     = seed
        self.modulus      = modulus                     # Size of the domain.
        self.random       = self._generateNumbers(size)  # The list of random numbers.
        self.lifetime     = 0                           # Lifetime of algorithm.

    def __str__(self):
        print('List of random numbers: \n {}'.format(self.random))

    # How to generate new numbers depends on the method. Leave it abstract in
    # the parent class.
    @abstractmethod
    def _randomGen(self):
        pass

    # Update the register, random list and no_generated. This is different for
    # the a shift register and a normal PRNG.
    def _update(self, new_num):
        if isinstance(self.register, int):
            self.register = new_num
        elif isinstance(self.register, np.ndarray):
            self.register     = np.roll(self.register, -1)
            self.register[-1] = new_num
        else:
            print('Error: the register should be either an int or an np.array')

    # Check if the latest random number (stored in the register) appears in the
    # list of random numbers.
    def _appeared(self, random, no_generated):
        # For now, only test the normal PRNG's. Shift registers have in general 
        # such a long lifetime that it it not necessary (and implementation is 
        # a bit tedious).
        if isinstance(self.register, int):
            return self.register in random[:no_generated]
        elif isinstance(self.register, np.ndarray):
            return False

    # Calculate the lifetime of this PRNG.
    def _generateNumbers(self, size):
        random       = np.zeros(size)
        no_generated = 0
        while no_generated < size:
            # Generate a new random number.
            newElement = self._randomGen()
            # Check if it already appeared in our list of random numbers.
            if self._appeared(random, no_generated):
                print(self._appeared(random, no_generated))
                # Truncate the list of random numbers (omit the final 0's).
                random = random[:no_generated]
                self.lifetime = no_generated
                break
            else:
                # Add new random number to the list
                random[no_generated] = newElement
                no_generated += 1
        return random

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
#   map and the linear congruential method.

# The Middle-Square method    
class middleSquare(prng):
    def __init__(self, seed, modulus, size):
        assert isinstance(seed, int)
        super().__init__(seed, modulus, size)

    def _randomGen(self):
        rand = (self.register**2 // self.modulus**(0.5)) % self.modulus
        #rand = int(np.floor(self.register**2 / self.modulus**(0.5)) % self.modulus)
        self._update(rand)
        return rand

# The Logistic Map method         
class logisticMap(prng):
    def __init__(self, seed, modulus, size, a = 4):
        assert isinstance(seed, int)
        self.a = a
        super().__init__(seed, modulus, size)

    def _randomGen(self):
        x_norm = self.register/self.modulus
        rand   = int(np.round(self.a*x_norm*(1-x_norm),
                              decimals = 4)*self.modulus)
        self._update(rand)
        return rand

# The Linear-Congruential method        
class linearCongruential(prng):
    def __init__(self, seed, modulus, size, a = 4001, c = 945):
        self.a = a
        self.c = c
        super().__init__(seed, modulus, size)
        
    def _randomGen(self):
        rand = (self.a*self.register + self.c) % self.modulus
        self._update(rand)
        return rand
    
#%% Implement an examples of shift register PRNG's.

# The RCARRY algorithm.
class rCarry(prng):
    def __init__(self, seed, modulus, size, s, carry_bit):
        assert isinstance(seed, np.ndarray)
        self.s = s
        self.r = len(seed)
        self.carry = carry_bit
        super().__init__(seed, modulus, size)
     
    def _randomGen(self):
        rand = self.register[-self.s] - self.register[-self.r] - self.carry
        if rand < 0:
            rand += self.modulus
            self.carry = 0
        else:
            self.carry = 1
        self._update(rand)
        return rand      

# A fibonacci shift register with s=1, r=2.
class fibonacci(prng):
    def __init__(self, seed, modulus, size):
        assert isinstance(seed, np.ndarray) and len(seed) == 2        
        self.s = 1
        self.r = 2
        super().__init__(seed, modulus, size)
    
    def _randomGen(self):
        rand = (self.register[-self.s] + self.register[-self.r]) % self.modulus
        self._update(rand)
        return rand

# A shift register with s=24 and r=55        
class s24r55(prng):
    def __init__(self, seed, modulus, size):
        assert isinstance(seed, np.ndarray) and len(seed) == 55        
        self.s = 24
        self.r = 55
        super().__init__(seed, modulus, size)
    
    def _randomGen(self):
        rand = (self.register[-self.s] + self.register[-self.r]) % self.modulus
        self._update(rand)
        return rand




















