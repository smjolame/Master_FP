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

def abw(exact,approx):
    return (exact-approx)*100/exact  #Abweichnung

#Kontrast:

phi, U_min, U_max = np.genfromtxt("data/kontrastmessung_Dummy.txt", delimiter=',', unpack=True)

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


#Brechungsindex Glas
durchgang, M = np.genfromtxt("data/Brechungsindex_Glas_Dummy.txt", delimiter=',', unpack=True) 
lambda_0 = 632.99 * 10**(-9) #m 
theta_0 = np.radians(10) #rad
D = 1 * 10**(-3) #m
theta_rad = theta_0


def n_func(M_,theta):
    return 1/(1-(M_*lambda_0)/(2*D*theta*theta_0))

n = n_func(M,theta_rad)
print(n)



##Curvefit
#def BeispielFunktion(x,a,b):
#    return a*x+b 
#params, cov = curve_fit(BeispielFunktion, x-Werte, y-Werte,sigma=fehler_der_y_werte,p0=[schätzwert_1,#schätzwert_2])
#a = ufloat(params[0],np.absolute(cov[0][0])**0.5)
#b = ufloat(params[1],np.absolute(cov[1][1])**0.5)
#
#
##Json
#Ergebnisse = json.load(open('data/Ergebnisse.json','r'))
#if not 'Name' in Ergebnisse:
#    Ergebnisse['Name'] = {}
#Ergebnisse['Name']['Name des Wertes']=Wert
#
#json.dump(Ergebnisse,open('data/Ergebnisse.json','w'),indent=4)
