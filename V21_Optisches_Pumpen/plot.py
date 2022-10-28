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

from scipy.constants import mu_0 , h , elementary_charge , electron_mass, hbar , eV

mu_B = elementary_charge*hbar/(2*electron_mass)

def abw(exact,approx):
    return (exact-approx)*100/exact  #Abweichnung


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

#Spulen-Eigenschaften:
R_s = 0.1639 #m
N_s = 11
R_h = 0.1579 #m
N_h = 154
R_v = 0.11735 #m
N_v = 20

#Konvention: Rb87 = 1, Rb85 = 2 

#Daten imortieren:
f_kilo, b_s1, b_s2, b_h1, b_h2 = np.genfromtxt('data/B_felder.csv', delimiter=',', unpack=True)
U1, T1 = np.genfromtxt('data/Periodendauern1.csv', delimiter=',', unpack=True)
U2, T2 = np.genfromtxt('data/Periodendauern2.csv', delimiter=',', unpack=True)

#Frequenz in Hz:
f = f_kilo*10**3

bruchteile_period1 = np.array([1,3,4,5,6,9,10,10,11,10,13,15,16,18])
bruchteile_period2 = np.array([1,2,3,4,5,7,8,8,9,10,10,8,10,8])

I_s1 , I_s2 = b_s1*0.1 , b_s2*0.1
I_h1 , I_h2 = b_h1*0.3 , b_h2*0.3



T1 = T1/bruchteile_period1 
T2 = T2/bruchteile_period2 




# B-Felder berechnen:

def B_helmholtz(N,R,I):
    return 8/(np.sqrt(125))*mu_0*N*I/R


B_s1 = B_helmholtz(N_s,R_s,I_s1)
B_s2 = B_helmholtz(N_s,R_s,I_s2)
B_h1 = B_helmholtz(N_h,R_h,I_h1)
B_h2 = B_helmholtz(N_h,R_h,I_h2)

B_1 = (B_s1+B_h1)#T
B_2 = (B_s2+B_h2)#T

B_v = B_helmholtz(N_v,R_v,2.3*0.1)#T

#Test-Plots

plt.plot(f,B_1, 'x')
plt.plot(f,B_2, 'x')
plt.grid()
plt.show()
plt.clf()
plt.plot(U1,T1,'x')
plt.plot(U2,T2,'x')
plt.grid()
plt.show()
plt.clf()

def Fit_Gerade (x,a,b):
    return a*x+b

f_lin = np.linspace(0,1100000,1000)

params, cov = curve_fit(Fit_Gerade, f,B_1 ,p0=[1,1])
a1 = ufloat(params[0],np.absolute(cov[0][0])**0.5)
b1 = ufloat(params[1],np.absolute(cov[1][1])**0.5)

params, cov = curve_fit(Fit_Gerade, f,B_2 ,p0=[1,1])
a2 = ufloat(params[0],np.absolute(cov[0][0])**0.5)
b2 = ufloat(params[1],np.absolute(cov[1][1])**0.5)


plt.plot(f,B_1, 'x')
plt.plot(f_lin, Fit_Gerade(f_lin,a1.n,b1.n))
plt.plot(f,B_2, 'x')
plt.plot(f_lin, Fit_Gerade(f_lin,a2.n,b2.n))
plt.grid()
plt.show()
plt.clf()

#g-Faktoren
def g_faktor(a):
    return h/(mu_B*a)

g1 = g_faktor(a1)
g2 = g_faktor(a2)

print(g1,g2)

#Erdmagnetfeld
print(b1,b2)

#Kernspins:
def I(g_F):
    return 1/g_F-0.5

I1 = I(g1)
I2 = I(g2)

print(I1,I2)


#Isotopenverhältnis
Isotope_Ratio = 1/2


#Quadratischer Zeeman Effekt:
def delta_E(B,g_F,M_F,E_Hyperfein):
    return g_F*mu_B*B+(1-2*M_F)*(g_F*mu_B*B)**2/E_Hyperfein


# Hyperfeinaufspaltungen von Rb87 und Rb85
E_Hyper_1 = 4.53*10**(-24) #J
E_Hyper_2 = 2.01*10**(-24) #J
del_E1 = delta_E(B_1[-1],g1,2,E_Hyper_1)
del_E2 = delta_E(B_2[-1],g2,3,E_Hyper_2)

del_E1_eV = del_E1/eV
del_E2_eV = del_E2/eV

print(del_E1_eV, del_E2_eV)


######################
#Exp Sättigung

def expo(t,A,b,c):
    return A+b*np.exp(c*t)



t1_ ,data1_ = np.genfromtxt('data/Exp_Anst_1.csv', delimiter=',', unpack=True)
t2_ ,data2_ = np.genfromtxt('data/Exp_Anst_2.csv', delimiter=',', unpack=True)
data1 = np.cumsum(data1_)
data2 = np.cumsum(data2_)
t1 = t1_*25/20 #ms
t2 = t2_*25/20 #ms

params_exp1, cov_exp1 = curve_fit(expo,t1, data1, p0=[0,0,-1])
A_e1 = ufloat(params_exp1[0],np.absolute(cov_exp1[0][0])**0.5)
b_e1 = ufloat(params_exp1[1],np.absolute(cov_exp1[1][1])**0.5)
c_e1 = ufloat(params_exp1[2],np.absolute(cov_exp1[2][2])**0.5)

params_exp2, cov_exp2 = curve_fit(expo,t2, data2, p0=[0,0,-1])
A_e2 = ufloat(params_exp2[0],np.absolute(cov_exp2[0][0])**0.5)
b_e2 = ufloat(params_exp2[1],np.absolute(cov_exp2[1][1])**0.5)
c_e2 = ufloat(params_exp2[2],np.absolute(cov_exp2[2][2])**0.5)



t_lin1 = np.linspace(t1[0],t1[-1],1000)
plt.plot(t_lin1, expo(t_lin1, A_e1.n, b_e1.n, c_e1.n))
plt.plot(t1,data1,'x')
t_lin2 = np.linspace(t2[0],t2[-1],1000)
plt.plot(t_lin2, expo(t_lin2, A_e2.n, b_e2.n, c_e2.n))
plt.plot(t2,data2,'x')
plt.grid()
plt.show()
plt.clf()


#Periodendauern

#Hyperbolische Fitfunktion:
def hyper(T,a,b,c):
    return a+b/(T-c)

params_period1, cov_period1 = curve_fit(hyper,U1, T1, p0=[0,3,0])
a_p1 = ufloat(params_period1[0],np.absolute(cov_period1[0][0])**0.5)
b_p1 = ufloat(params_period1[1],np.absolute(cov_period1[1][1])**0.5)
c_p1 = ufloat(params_period1[2],np.absolute(cov_period1[2][2])**0.5)


params_period2, cov_period2 = curve_fit(hyper,U2, T2, p0=[0,3,0])
a_p2 = ufloat(params_period2[0],np.absolute(cov_period2[0][0])**0.5)
b_p2 = ufloat(params_period2[1],np.absolute(cov_period2[1][1])**0.5)
c_p2 = ufloat(params_period2[2],np.absolute(cov_period2[2][2])**0.5)


U_lin1 = np.linspace(U1[0],U1[-1],1000)
plt.plot(U1,T1,'x')
plt.plot(U_lin1,hyper(U_lin1,a_p1.n,b_p1.n,c_p1.n))
U_lin2 = np.linspace(U2[0],U2[-1],1000)
plt.plot(U2,T2,'x')
plt.plot(U_lin2,hyper(U_lin2,a_p2.n,b_p2.n,c_p2.n))
plt.grid()
plt.show()
print(a_p1, b_p1, c_p1)
print(a_p2, b_p2, c_p2)