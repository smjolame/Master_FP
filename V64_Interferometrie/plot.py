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

print('A_fit:', A_fit)

plt.plot(phi_rad,K,'x',label="Messwerte")
#plt.xlabel()
#plt.ylabel()
plt.plot(phi_rad_lin,K_phi(phi_rad_lin,A_fit.n), label="Ausgleichsrechnung")
plt.legend(loc=1)
tick_pos= [0, np.pi/2 , np.pi]
labels = ['0', r'$\frac{\pi}{2}$', r'$\pi$']
plt.xticks(tick_pos,labels)
plt.xlabel(r'$\Theta \mathbin{/} \si{\radian}$')
plt.ylabel(r'$K$')
plt.ylim([0,1.4])
plt.grid()
plt.tight_layout()
plt.savefig("build/Kontrast.pdf")
plt.clf()


#Brechungsindex Glas
durchgang, M = np.genfromtxt("data/Brechungsindex_Glas_Dummy.txt", delimiter=',', unpack=True) 
lambda_0 = 632.99 * 10**(-9) #m 
theta_0 = np.radians(10) #rad
D = 1 * 10**(-3) #m
theta_rad = theta_0


def n_func(M_,theta):
    return 1/(1-(M_*lambda_0)/(2*D*theta*theta_0))

n_glas = n_func(M,theta_rad)
np.savetxt('build/n_Glas.txt',n_glas.T, delimiter = '  &  ')
#########Brechungsindex Gas:
####Luft
L=ufloat(0.1,0.0001) #m
T = 21.1+273.15 #K
p, M1,M2,M3,M4 = np.genfromtxt("data/Brechungsindex_Luft_Dummy.txt", delimiter=',', unpack=True)


def n_func2(M_):
    return (M_*lambda_0)/L+1


n_luft1 = n_func2(M1)
n_luft2 = n_func2(M2)
n_luft3 = n_func2(M3)
n_luft4 = n_func2(M4)


#Nominal werte
n_luft1_n = unp.nominal_values(n_luft1)
n_luft2_n = unp.nominal_values(n_luft2)
n_luft3_n = unp.nominal_values(n_luft3)
n_luft4_n = unp.nominal_values(n_luft4)

n_luft_Tabelle_n = np.array([n_luft1_n,n_luft2_n,n_luft3_n,n_luft4_n]).T

# n Daten abspeichern
np.savetxt('build/n_Luft.txt',n_luft_Tabelle_n, delimiter = '  &  ')


####Gas
R = const.R
def n_func_gas(p_,a,b):
    return a/(T*R)*p_+b

##### curvefits für Gas
params1, cov1 = curve_fit(n_func_gas, p, n_luft1_n)
a_gas1 = ufloat(params1[0],np.absolute(cov1[0][0])**0.5)
b_gas1 = ufloat(params1[1],np.absolute(cov1[1][1])**0.5)
print('a_gas_1:',a_gas1,'b_gas_1:',b_gas1)
params2, cov2 = curve_fit(n_func_gas, p, n_luft2_n)
a_gas2 = ufloat(params2[0],np.absolute(cov2[0][0])**0.5)
b_gas2 = ufloat(params2[1],np.absolute(cov2[1][1])**0.5)
print('a_gas_2:',a_gas2,'b_gas_2:',b_gas2)
params3, cov3 = curve_fit(n_func_gas, p, n_luft3_n)
a_gas3 = ufloat(params3[0],np.absolute(cov3[0][0])**0.5)
b_gas3 = ufloat(params3[1],np.absolute(cov3[1][1])**0.5)
print('a_gas_3:',a_gas3,'b_gas_3:',b_gas3)
params4, cov4 = curve_fit(n_func_gas, p, n_luft4_n)
a_gas4 = ufloat(params4[0],np.absolute(cov4[0][0])**0.5)
b_gas4 = ufloat(params4[1],np.absolute(cov4[1][1])**0.5)
print('a_gas_4:',a_gas4,'b_gas_4:',b_gas4)



######## Plots für Gas
p_lin=np.linspace(p[0], p[-1],1000)

plt.plot(p,unp.nominal_values(n_luft1), 'x', label = 'Messwerte 1', c = 'b')
plt.plot(p_lin,n_func_gas(p_lin,a_gas1.n,b_gas1.n), label = 'Fit 1', c = 'b')

plt.plot(p,unp.nominal_values(n_luft2), 'x', label = 'Messwerte 2',c = 'm')
plt.plot(p_lin,n_func_gas(p_lin,a_gas2.n,b_gas2.n), label = 'Fit 2', c = 'm')

plt.plot(p,unp.nominal_values(n_luft3), 'x', label = 'Messwerte 3',c = 'k')
plt.plot(p_lin,n_func_gas(p_lin,a_gas3.n,b_gas3.n), label = 'Fit 3',c = 'k')

plt.plot(p,unp.nominal_values(n_luft4), 'x', label = 'Messwerte 4',c = 'r')
plt.plot(p_lin,n_func_gas(p_lin,a_gas4.n,b_gas4.n), label = 'Fit 4',c = 'r')
plt.xlabel(r'$p \mathbin{/} \si{\milli\bar}$')
plt.ylabel(r'$n$')
plt.legend(fontsize='small')
plt.tight_layout()
plt.grid()
plt.savefig("build/Gas.pdf")

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
