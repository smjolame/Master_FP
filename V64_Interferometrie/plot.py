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

phi, U_min, U_max = np.genfromtxt("data/kontrast.txt", delimiter=',', unpack=True)

phi_rad = np.radians(phi)#rad

def K(U_min_,U_max_):
    return (U_max_- U_min_)/(U_max_ + U_min_)

def K_phi(phi,A):
    return A*np.abs(np.cos(phi)*np.sin(phi))

K = K(U_min, U_max)

np.savetxt('build/K.txt',K.T, fmt='%1.2f' , delimiter = '  &  ')


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
plt.xlabel(r'$\phi \mathbin{/} \si{\radian}$')
plt.ylabel(r'$K$')
plt.ylim([0,1.4])
plt.grid()
plt.tight_layout()
plt.savefig("build/Kontrast.pdf")
plt.clf()


#Brechungsindex Glas
durchgang, M = np.genfromtxt("data/Brechungsindex_Glas.txt", delimiter=',', unpack=True) 
lambda_0 = 632.99 * 10**(-9) #m 
theta_0 = np.radians(10) #rad
D = 1 * 10**(-3) #m
theta_rad = theta_0


def n_func(M_,theta):
    return 1/(1-(M_*lambda_0)/(2*D*theta*theta_0))

n_glas = n_func(M,theta_rad)
np.savetxt('build/n_Glas.txt',n_glas.T, fmt='%1.8f' , delimiter = '  &  ')


n_glas_mean = np.mean(n_glas)
n_glas_err = np.std(n_glas)
n_glas_exp=np.array([n_glas_mean,n_glas_err])
n_glas_u = ufloat(n_glas_mean,n_glas_err)
np.savetxt('build/n_Glas_mean.txt',n_glas_exp, fmt='%1.8f' , delimiter = '  &  ')


#########Brechungsindex Gas:
####Luft
L=ufloat(0.1,0.0001) #m
T = 20.7+273.15 #K
p, M1,M2,M3 = np.genfromtxt("data/Brechungsindex_Luft.txt", delimiter=',', unpack=True)


def n_func2(M_):
    return (M_*lambda_0)/L+1


n_luft1 = n_func2(M1)
n_luft2 = n_func2(M2)
n_luft3 = n_func2(M3)







#Nominal werte
n_luft1_n = unp.nominal_values(n_luft1)
n_luft2_n = unp.nominal_values(n_luft2)
n_luft3_n = unp.nominal_values(n_luft3)

n_luft1_mean = np.mean(unp.nominal_values(n_luft1))
n_luft2_mean = np.mean(unp.nominal_values(n_luft2))
n_luft3_mean = np.mean(unp.nominal_values(n_luft3))


n_luft1_std = np.std(unp.nominal_values(n_luft1))
n_luft2_std = np.std(unp.nominal_values(n_luft2))
n_luft3_std = np.std(unp.nominal_values(n_luft3))

n_luft = ufloat(np.mean([n_luft1_mean,n_luft2_mean,n_luft3_mean]),np.mean([n_luft1_std,n_luft2_std,n_luft3_std]))
print('n_luft',n_luft)
n_luft_mean = np.array([n_luft1_mean,n_luft1_std,n_luft2_mean,n_luft2_std,n_luft3_mean,n_luft3_std])


n_luft_Tabelle_n = np.array([n_luft1_n,n_luft2_n,n_luft3_n]).T

# n Daten abspeichern
np.savetxt('build/n_Luft.txt',n_luft_Tabelle_n,fmt='%1.8f', delimiter = '  &  ')
np.savetxt('build/n_Luft_mean.txt',n_luft_mean.T,fmt='%1.8f', delimiter = '  &  ')

####Gas
R = const.R
def n_func_gas(p_,a,b):
    return np.sqrt(1 + a*p_/(R*T)) + b


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



######## Plots für Gas
p_lin=np.linspace(p[0], p[-1],1000)



plt.plot(p,unp.nominal_values(n_luft1), 'x', label = 'Messwerte 1', c = 'b')
plt.plot(p_lin,n_func_gas(p_lin,a_gas1.n,b_gas1.n), label = 'Fit 1', c = 'b')

plt.plot(p,unp.nominal_values(n_luft2), 'x', label = 'Messwerte 2',c = 'r')
plt.plot(p_lin,n_func_gas(p_lin,a_gas2.n,b_gas2.n), label = 'Fit 2', c = 'r')

plt.plot(p,unp.nominal_values(n_luft3), 'x', label = 'Messwerte 3',c = 'k')
plt.plot(p_lin,n_func_gas(p_lin,a_gas3.n,b_gas3.n), label = 'Fit 3',c = 'k')

plt.xlabel(r'$p \mathbin{/} \si{\milli\bar}$')
plt.ylabel(r'$n$')
plt.legend(fontsize='small')
plt.tight_layout()
plt.grid()
plt.savefig("build/Gas.pdf")


#######abw

nL = 1.0003
nG = 1.45
print(nL,n_luft)
print(nG,n_glas_u)
print(abw(nL,n_luft))
print(abw(nG,n_glas_u))

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
