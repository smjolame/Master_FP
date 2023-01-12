import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import json
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import uarray
from uncertainties import unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,std_devs as stds)
from scipy.stats import sem


phi, U_min, U_max = np.genfromtxt("data/kontrast.txt",  unpack=True)

phi_rad = np.radians(phi)#rad

def K(U_min_,U_max_):
    return (U_max_- U_min_)/(U_max_ + U_min_)

def K_phi(phi,A):
    return A*np.abs(np.cos(phi)*np.sin(phi))

K = K(U_min, U_max)



params, cov = curve_fit(K_phi,phi_rad,K)
phi_rad_lin = np.linspace(0,np.pi,1000)

A_fit = ufloat(params,np.absolute(cov)**0.5)

plt.plot(phi_rad,K,'x',label="Messwerte")
#plt.xlabel()
#plt.ylabel()
plt.plot(phi_rad_lin,K_phi(phi_rad_lin,A_fit.n), label="Ausgleichsrechnung")
plt.legend()
plt.savefig("build/Kontrast.pdf")
plt.clf()

print(np.amax(K))
print(phi[np.where(K == np.amax(K))])

L=ufloat(0.1,0.0001) #m
T = 21.1+273.15 #K
p, M1 = np.genfromtxt("data/Brechungsindex_Luft.txt", delimiter=',', unpack=True)

lambda_0 = 632.99 * 10**(-9) #m 

def n_func2(M_):
    return (M_*lambda_0)/L+1
    

n_luft1 = n_func2(M1)
n_luft1_n = unp.nominal_values(n_luft1)

plt.plot(p,unp.nominal_values(n_luft1), 'x', label = 'Messwerte 1', c = 'b')
plt.savefig("build/Gas.pdf")