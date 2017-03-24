# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 23:20:57 2017

@author: anie
"""

import numpy as np
import matplotlib.pyplot as plt

def g(x):
    return x**2*(2*x-1)

def T(cutoff, n=1000):
    if n <= cutoff:
        return g(n)
    if n%2==0:
        return 7*T(cutoff, n/2)+18*(n/2)**2
    else:
        #return T(cutoff, n+1)
        return 6 * T(cutoff, (n+1)/2) + T(cutoff, (n-1)/2) + 6*((n-1)/2)**2+8*(n**2-1)/4 + 4*((n+1)/2)**2
if __name__ == "__main__":
    bests = np.array([])
    ns = np.arange(10,15)
    for n in ns:
        t_best = np.array([])
        #t_worst = np.array([])
        for cutoff in np.arange(5,50):
            t_best = np.append(t_best,T(cutoff, 2**n+1))
            bests = np.append(bests, np.argmin(t_best))
            #t_worst = np.append(t_worst, T(cutoff, 2**n+1))
        plt.plot(np.log(t_best), '--')
        #plt.plot(np.log(t_worst), 'o')
    fig = plt.figure()
    plt.plot(bests)
