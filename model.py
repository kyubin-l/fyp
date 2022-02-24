import numpy as np
import matplotlib.pyplot as plt
import os
import sympy as sp
import h5py

from scipy.integrate import quad
from scipy.io import loadmat

j = complex(0, 1)
e = np.e

class TemperatureGradientModel:
    gamma = 1.4
    Rg = 287
    Rt = 0.001
    def __init__(self, L, x1, n, p1, M1, T1, T2, temp_dist = 'lin', Z1 = -1, fl_c = 1, viscous_effects=False):
        self.L = L
        self.x1 = x1
        self.p1 = p1
        self.x = np.linspace(x1, L, n)
        self.dx = self.x[1] - self.x[0]
        self.M1 = M1
        self.T1 = T1
        self.T2 = T2
        self.Z1 = Z1
        self.temp_dist = temp_dist
        self.fl_c = fl_c
        
        if self.temp_dist == 'lin':
            self.T = self.T1 + (self.T2-self.T1)/self.L*(self.x-self.x1)

        elif self.temp_dist == 'sin':
            self.T = (self.T1-self.T2)/2*np.sin(5*np.pi/4*(self.x-self.x1)/self.L + np.pi/4)+(self.T1+self.T2)/2

        self.c1 = np.sqrt(self.gamma * self.Rg * self.T[0])
        self.u1 = self.M1 * self.c1
        self.rho1 = self.p1/(self.Rg*self.T[0])
        
        self.A = self.rho1*self.u1
        self.B = self.p1 + self.rho1*self.u1**2
        
        self.w = self.fl_c*2*np.pi*self.c1/self.L
        self.viscous_effects = viscous_effects
        
    def calculate(self):
        self.calculate_mean_properties()
        self.calculate_flow_properties()
        
    
    def _grad(self, q):
        dqdx = np.zeros(q.shape)
        dqdx[1:-1] = (q[2:]-q[:-2])/(2*self.dx)
        dqdx[0] = (q[1]-q[0])/self.dx
        dqdx[-1] = (q[-1]-q[-2])/self.dx
        return dqdx
        
        
    def calculate_mean_properties(self):
        self.c = np.sqrt(self.gamma * self.Rg * self.T)
        self.u = (self.B/self.A - np.sqrt((self.B/self.A)**2 - 4*self.Rg * self.T))/2
        self.rho = self.A/self.u
        self.p = self.rho*self.Rg*self.T
        self.M = self.u/self.c
        self.alpha = 1/self.rho*self._grad(self.rho)
        
        def Dynamic_viscosity(T):
            if (100 <= T <= 1000):
                return 1.716e-5*(T/273.15)**(3/2)*383.55/(T+110.4)
            elif(1000 < T <= 3000):
                return 2.653e-8*T+1.573e-5
        
        self.v = np.array([Dynamic_viscosity(t) for t in self.T])/self.rho
        self.Pr = 0.71
        self.Sh = self.Rt * (self.w/self.v)**0.5
        
        if self.viscous_effects:
            self.k0 = self.w/(self.c*(1+self.M))*(1+(1-j)/np.sqrt(2)/self.Sh*(1+(self.gamma-1)/(self.Pr**(1/2)))-j/(self.Sh**2)*(1+(self.gamma-1)/(self.Pr**(1/2))-self.gamma/2*(self.gamma-1)/self.Pr))
        else:
            self.k0 = self.w/self.c
        return
    
    
    def plot_mean_properties(self):
        return
        
        
    def _numerical_integration(self):
        
        def c_T(T):
            return np.sqrt(self.gamma*self.Rg*T)

        def u_T(T):
            return (self.B/self.A - np.sqrt((self.B/self.A)**2-4*self.Rg*T))/2

        def T_x(x):
            if self.temp_dist == 'lin':
                return self.T1+(self.T2-self.T1)*(x-self.x1)/self.L
            else:
                return (self.T1-self.T2)/2*np.sin(5*np.pi/4*(x-self.x1)/self.L + np.pi/4)+(self.T1+self.T2)/2
        
        def integrand_plus(x):
            return 1/(c_T(T_x(x)) + u_T(T_x(x)))

        def integrand_minus(x):
            return 1/(c_T(T_x(x)) - u_T(T_x(x)))
        
        def integral_plus(x):
            return quad(integrand_plus, self.x1, x)[0]

        def integral_minus(x):
            return quad(integrand_minus, self.x1, x)[0]
        
        self.integrals_plus = np.array([integral_plus(value) for value in self.x])
        self.integrals_minus = np.array([integral_minus(value) for value in self.x])
        
    
    def _get_intermediate_vals(self, rho, M, rho1, M1, gamma, w, c, alpha, integrals_plus, integrals_minus):
        self.P1_plus = (rho/rho1)**(1/4)*(1+M1)/(1+M)*(np.exp(gamma*M1-gamma/4*M1**2-(gamma**2-1)/3*M1**3))/(np.exp(gamma*M-gamma/4*M**2-(gamma**2-1)/3*M**3))*np.exp(-j*w*integrals_plus)
        self.P1_minus= (rho/rho1)**(1/4)*(1-M1)/(1-M)*(np.exp(gamma*M+gamma/4*M**2-(gamma**2-1)/3*M**3))/(np.exp(gamma*M1+gamma/4*M1**2-(gamma**2-1)/3*M1**3))*np.exp(j*w*integrals_minus)
        self.K1 = (j*self.k0-(1+2*(1+gamma)*M+(3*gamma-7)*M**2)*alpha/4)/(j*self.k0-alpha*M)
        self.K2 = (j*self.k0+(1-2*(1+gamma)*M+(3*gamma-7)*M**2)*alpha/4)/(j*self.k0-alpha*M)
        
        
    def _get_flow_properties(self, C1_plus, C1_minus, P1_plus, P1_minus, K1, K2, rho, c):
        self.p_hat = C1_plus*P1_plus + C1_minus*P1_minus
        self.u_hat = K1*C1_plus*P1_plus/(rho*c) - K2*C1_minus*P1_minus/(rho*c)
        self.Fp = self.p_hat/(rho[0]*c[0]*self.u_hat[0])
        self.Fu = self.u_hat/self.u_hat[0]
        
    
    def calculate_flow_properties(self):
        self._numerical_integration()
        self._get_intermediate_vals(self.rho, self.M, self.rho1, self.M1, self.gamma, self.w, self.c, self.alpha, self.integrals_plus, self.integrals_minus)
        
        # Initial conditions don't matter, plots later on are initialised
        # self.p0 = 100
        # self.u0 = self.p0/(self.rho[0]*self.c[0]*self.Z1)
        
        self.u0 = 1
        self.p0 = self.rho[0] * self.c[0] * self.Z1 * self.u0
        
        self.C1_plus = self.p0*(self.K2[0] + 1/self.Z1)/(self.K1[0] + self.K2[0])
        self.C1_minus = self.p0*(self.K1[0] - 1/self.Z1)/(self.K1[0] + self.K2[0])
        
        self._get_flow_properties(self.C1_plus, self.C1_minus, self.P1_plus, self.P1_minus, self.K1, self.K2, self.rho, self.c)
        
        
    def plot_flow_properties(self):
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize = (14, 8))
        ax1.plot(self.x/self.L, abs(self.Fp))
        ax2.plot(self.x/self.L, abs(self.Fu))
        ax3.plot(self.x/self.L, np.angle(self.Fp)/np.pi)
        ax4.plot(self.x/self.L, np.angle(self.Fu)/np.pi)
        ax1.set_xlabel('x/l')
        ax2.set_xlabel('x/l')
        ax3.set_xlabel('x/l')
        ax4.set_xlabel('x/l')
        ax1.set_ylabel('|Fp|')
        ax2.set_ylabel('|Fu|')
        ax3.set_ylabel('<Fp/pi')
        ax4.set_ylabel('<Fu/pi')
        fig.show()
        
    
    def plot_mean_properties(self):
        # Need to plot temperature, velocity, mach number, speed of sound, density, pressure
        fig, ax = plt.subplots(ncols=2, nrows=3, figsize = (10, 10), constrained_layout=True, sharex=True)
        titles = ['Temperature', ]
        ax[0][0].plot(self.x, self.T)
        ax[0][1].plot(self.x, self.u)
        ax[1][0].plot(self.x, self.M)
        ax[1][1].plot(self.x, self.c)
        ax[2][0].plot(self.x, self.rho)
        ax[2][1].plot(self.x, self.p/100000)
        
        
        ax[0][0].set_title('Temperature (K)')
        ax[0][1].set_title('Velocity (m/s)')
        ax[1][0].set_title('Mach Number')
        ax[1][1].set_title('Speed of sound (m/s)')
        ax[2][0].set_title(r'Density ($\mathrm{kgm^{-3}}$)')
        ax[2][1].set_title(r'Pressure $10^{5} \mathrm{Pa}$')
        
        for axes in ax.flat:
            axes.grid(visible=True, which='major', linestyle='--')
        fig.supxlabel('x/L')