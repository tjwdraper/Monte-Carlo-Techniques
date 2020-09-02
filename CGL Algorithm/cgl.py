# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 10:47:39 2020

This program test the CGL algoritm to approximate the integral of a simple test 
function. Throughout this file, I will use uniformly distributed random numbers

@author: Tom Draper (s4468201)
"""

import numpy as np
import matplotlib.pyplot as plt

#%% Define the CGL algorithm

class cgl():
    def __init__(self, point_set_size):
        self.vars = [0,0,0,0] # Initially [M,P,Q,R] = [0,0,0,0]
        self.count = 0 # keeps track of the number of updates
        self.pss = point_set_size
        self.estimators = np.zeros((point_set_size,3))
        
    def update_step(self, ew):
        # Update the internal parameters [M,P,Q,R] with 1 event_weight 
        # following the CGL algorithm
        [m,p,q,r] = self.vars
        N = self.count+1
        u = ew - m
        self.vars = [m+u/N,
                     (N-1)*(p+u**2/N)/N,
                     (N-1)*(q+(N-2)*u**3/N**2-3*p*u/N)/N,
                     (N-1)*(r+(p-(N-2)*u**2/N)**2/N-4*(q*u/N-p*u**2/N**2))/N]
        self.estimators[N-1] = [self.vars[0],
                                self.vars[1]/N,
                                self.vars[3]/N**3]
        self.count += 1
        
    def update(self, ews, show_progress = True):
        for ew in ews:
            self.update_step(ew)
        print('J1 = {} ± ({} ± {})'.format(self.estimators[self.pss-1,0],
                                                  self.estimators[self.pss-1,1]**0.5,
                                                  self.estimators[self.pss-1,2]**0.25))
        if show_progress:
            x = np.linspace(1, self.pss, self.pss)
            fig, (ax1, ax2) = plt.subplots(1,2, figsize = (15,8))
            # First subplot: show estimator of the integral
            ax1.set_title('Estimate of the integral', fontsize = 18)
            ax1.tick_params(labelsize = 14)
            ax1.set_xlabel('$\log(N)$', fontsize = 16)
            ax1.set_ylabel('Integral Estimator $\hat{E}_1$', fontsize = 16)
            ax1.semilogx(x,
                         self.estimators[:,0])
            # Show the relative errors of the estimates in a loglog plot.
            ax2.set_title('Estimate of the relative errors', fontsize = 18)
            ax2.tick_params(labelsize = 14)
            ax2.set_xlabel('$\log(N)$', fontsize = 16)
            ax2.set_ylabel('Relative Errors', fontsize = 16)
            ax2.loglog(x[1:],
                       self.estimators[1:,1]**0.5/self.estimators[1:,0],
                       label = '$\hat{E}_2^{1/2}/\hat{E}_1$')
            ax2.loglog(x[2:],
                       self.estimators[2:,2]**0.25/self.estimators[2:,1]**0.5,
                       label = '$\hat{E}_4^{1/4}/\hat{E}_2^{1/2}$')
            ax2.legend(loc = 'lower left', fontsize = 14)
            plt.show()
        

#%% Define a test function

# NOTE: I don't have to include the step function theta here. That has 
# automatically been taken care of by sampling random numbers in our desired 
# region Γ

def f(a,x):
    return (1+a)*x**a

#%% Set our parameters. The low and upper bound define our integration region Γ. 

low_bound = 0
up_bound = 1
a = -0.9
point_set_size = 10000

#%% Test the test function

# Initialize the parameters of the CGL algorithm.
alg = cgl(point_set_size)

# The probability density of a uniform distribution is just unity, so f alone
# generates to event_weights. This is of course different if the random numbers
# are sampled from another distribution.
event_weights = f(a,
                  np.random.uniform(low = low_bound, 
                                    high = up_bound, 
                                    size = point_set_size))

# Update the estimators with all event weights.
alg.update(event_weights)
            