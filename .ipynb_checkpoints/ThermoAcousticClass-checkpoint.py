# importing libraries
import numpy as np
import matplotlib.pyplot as plt
import os
import sympy as sp

j = complex(0, 1)
e = np.e

class ThermoacousticModel():
    '''
    Thermo acoustic model class. 
    '''
    
    j = complex(0, 1)
    e = np.e
    
    def __init__(self, c, L, Ru, Rd, damping=True, temperature_gradient=True):
        self.c = c
        self.L = L
        self.Ru = Ru
        self.Rd = Rd
        self.damping = damping
        self.temperature_gradient = temperature_gradient
        
    
    def matrix_det(self, w, Ru, Rd):
        det = -sp.exp(j*w*self.L/self.c) + Rd*Ru*sp.exp(-j*w*self.L/self.c)
        return det
        
        
    def matrix_detzero(self, verbose=True):
        '''
        Find the value of frequency w (in rad/s) which makes the 
        determinent of the matrix M 0. 
        '''
        
        w = sp.Symbol('w')
        wl = sp.solveset(self.matrix_det(self, w, self.Ru, self.Rd), w, domain=sp.S.Complexes)
        

        # Infinite set into generator
        wl_iterable = iter(wl)
        # Taking the first 100 solutions, can take more if wanted. Goes in order 0, 1, -1, 2, -2, ...
        wl_all = [complex(next(wl_iterable)) for i in range(100)] 

        # taking positive values, want n = 1, 2, 3..., aka the odd indices in the wl_all list
        w_n = [wl_all[i] for i in range(len(wl_all)) if (i%2 == 1)]
        w_n = np.asarray(w_n)
        
        if verbose==True:
            sp.pprint(wl)
            
        return w_n
    
    
    def p_fluc(self, wl):
        self.x = np.linspace(0, self.L, 1000)[..., np.newaxis]
        return np.exp((j*wl*(-self.x/self.c))) + 1/self.Ru*e**(j*wl*self.x/self.c) # using BC and x = 0
        #return np.exp((j*wl*(-x/c))) + Rd*e**(-2*j*wl*L/c)*e**(j*wl*x/c) # using BC at x = L
    
    
    def plot_modes(self, n):
        w_n = self.matrix_detzero()
        p = self.p_fluc(w_n)
        
        labels = [str(f"n={i}") for i in range(1, n+1)]
        plt.plot(self.x, np.abs(p[:,:n]), label = labels)
        plt.ylabel('|p\'|')
        plt.xlabel('Length')
        plt.legend(loc='upper left')
        plt.show()
        