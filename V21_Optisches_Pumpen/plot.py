from cgi import print_form
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
import matplotlib as mpl

mpl.use('pgf')
mpl.rcParams.update({
'font.family': 'serif',
'text.usetex': True,
'pgf.rcfonts': False,
'pgf.texsystem': 'lualatex',
'pgf.preamble': r'\usepackage{unicode-math}\usepackage{siunitx}',
})

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

#Daten importieren:
f_kilo, b_s1, b_s2, b_h1, b_h2 = np.genfromtxt('data/B_felder.csv', delimiter=',', unpack=True)
U1, T1 = np.genfromtxt('data/Periodendauern1.csv', delimiter=',', unpack=True)
U2, T2 = np.genfromtxt('data/Periodendauern2.csv', delimiter=',', unpack=True)

#Frequenz in Hz:
f = f_kilo*10**3

bruchteile_period1 = np.array([1,3,4,5,6,9,10,10,11,10,13,15,16,18])
bruchteile_period2 = np.array([1,2,3,4,5,7,8,8,9,10,10,8,10,8])

I_s1 , I_s2 = b_s1*0.1 , b_s2*0.1
I_h1 , I_h2 = b_h1*0.3 , b_h2*0.3

I_data = np.array([I_s1,I_s2,I_h1,I_h2]).T


np.savetxt('build/Data_Strom.txt', I_data,'%6.2f')

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
np.savetxt('build/Data_B_Felder.txt',np.array([B_1,B_2]).T*10**6,'%6.2f')

B_v = B_helmholtz(N_v,R_v,2.3*0.1)#T
print('B_v=', B_v*10**6,'muT')



#Bestimmung des g Faktors und des Erdmagnetfeldes
def Fit_Gerade (x,a,b):
    return a*x+b

f_lin = np.linspace(0,1100,1000) #kHz

params, cov = curve_fit(Fit_Gerade, f,B_1 ,p0=[1,1])
a1_ = ufloat(params[0],np.absolute(cov[0][0])**0.5)
b1_ = ufloat(params[1],np.absolute(cov[1][1])**0.5)

params, cov = curve_fit(Fit_Gerade, f,B_2 ,p0=[1,1])
a2_ = ufloat(params[0],np.absolute(cov[0][0])**0.5)
b2_ = ufloat(params[1],np.absolute(cov[1][1])**0.5)

a1=a1_*10**9 #muT/kHz
b1=b1_*10**6 #muT
a2=a2_*10**9 #muT/kHz
b2=b2_*10**6 #muT


print('a1=',a1,'muT/kHz')
print('b1=',b1,'muT')
print('a2=',a2,'muT/kHz')
print('b2=',b2,'muT')

plt.plot(f*10**(-3),B_1*10**(6), '+', label=r'Messwerte $^{87}$Rb')
plt.plot(f_lin, Fit_Gerade(f_lin,a1.n,b1.n), '--' ,label=r'Ausgleichsgerade $^{87}$Rb')
plt.plot(f*10**(-3),B_2*10**(6), 'x', label=r'Messwerte $^{85}$Rb')
plt.plot(f_lin, Fit_Gerade(f_lin,a2.n,b2.n), label=r'Ausgleichsgerade $^{85}$Rb')
plt.xlabel(r'$f \mathbin{/} \si{\kilo\hertz}$')
plt.ylabel(r'$B \mathbin{/} \si{\micro\tesla}$')
plt.legend()
plt.grid()
plt.savefig('build/B_Felder.pdf')
plt.clf()

#g-Faktoren
def g_faktor(a):
    return h/(mu_B*a)

g1 = g_faktor(a1_)
g2 = g_faktor(a2_)

print('g1=',g1)
print('g2=',g2)

#Erdmagnetfeld
print('b1=',b1,' muT')
print('b2=',b2,' muT')

#Kernspins:
def I(g_F):
    return 1/g_F-0.5

I1 = I(g1)
I2 = I(g2)

print('I1=',I1)
print('I2=',I2)


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

del_E1_eV = del_E1/eV*10**9
del_E2_eV = del_E2/eV*10**9
print('B_del1=',B_1[-1]*10**6,'muT')
print('B_del2=',B_1[-2]*10**6,'muT')
print('del_E1_eV=',del_E1_eV,'neV')
print('del_E2_eV=',del_E2_eV,'neV')


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

np.savetxt('build/expo.txt', np.array([t1,data1,data2]).T,'%3.2f')

params_exp1, cov_exp1 = curve_fit(expo,t1, data1, p0=[0,0,-1])
A_e1 = ufloat(params_exp1[0],np.absolute(cov_exp1[0][0])**0.5)
b_e1 = ufloat(params_exp1[1],np.absolute(cov_exp1[1][1])**0.5)
c_e1 = ufloat(params_exp1[2],np.absolute(cov_exp1[2][2])**0.5)

params_exp2, cov_exp2 = curve_fit(expo,t2, data2, p0=[0,0,-1])
A_e2 = ufloat(params_exp2[0],np.absolute(cov_exp2[0][0])**0.5)
b_e2 = ufloat(params_exp2[1],np.absolute(cov_exp2[1][1])**0.5)
c_e2 = ufloat(params_exp2[2],np.absolute(cov_exp2[2][2])**0.5)



t_lin1 = np.linspace(t1[0],t1[-1],1000)
plt.plot(t1,data1,'+', label='Messwerte $^{87}$Rb')
plt.plot(t_lin1, expo(t_lin1, A_e1.n, b_e1.n, c_e1.n),'--' ,label='Exponentieller Fit $^{87}$Rb')
t_lin2 = np.linspace(t2[0],t2[-1],1000)
plt.plot(t2,data2,'x', label='Messwerte $^{85}$Rb')
plt.plot(t_lin2, expo(t_lin2, A_e2.n, b_e2.n, c_e2.n), label='Exponentieller Fit $^{85}$Rb')
plt.grid()
plt.xlabel(r'$t \mathbin{/} \si{\milli\s}$')
plt.ylabel(r'$U \mathbin{/} \text{a.u.}$')
plt.legend()
plt.savefig('build/Expo.pdf')
plt.clf()


#Periodendauern

#Hyperbolische Fitfunktion:
def hyper(U,a,b,c):
    return a+b/(U-c)

params_period1, cov_period1 = curve_fit(hyper,U1, T1, p0=[0,3,0])
a_p1 = ufloat(params_period1[0],np.absolute(cov_period1[0][0])**0.5)
b_p1 = ufloat(params_period1[1],np.absolute(cov_period1[1][1])**0.5)
c_p1 = ufloat(params_period1[2],np.absolute(cov_period1[2][2])**0.5)


params_period2, cov_period2 = curve_fit(hyper,U2, T2, p0=[0,3,0])
a_p2 = ufloat(params_period2[0],np.absolute(cov_period2[0][0])**0.5)
b_p2 = ufloat(params_period2[1],np.absolute(cov_period2[1][1])**0.5)
c_p2 = ufloat(params_period2[2],np.absolute(cov_period2[2][2])**0.5)

np.savetxt('build/period.txt', np.array([U1,T1,T2]).T,'%3.2f')
print(a_p1,b_p1,c_p1,a_p2,b_p2,c_p2, b_p2/b_p1)

U_lin1 = np.linspace(U1[0],U1[-1],1000)
plt.plot(U1,T1,'+', label='Messwerte Periodendauer $^{87}$Rb')
plt.plot(U_lin1,hyper(U_lin1,a_p1.n,b_p1.n,c_p1.n),'--',label='Hyperbolischer Fit $^{87}$Rb')
U_lin2 = np.linspace(U2[0],U2[-1],1000)
plt.plot(U2,T2,'x', label='Messwerte Periodendauer $^{85}$Rb')
plt.plot(U_lin2,hyper(U_lin2,a_p2.n,b_p2.n,c_p2.n), label='Hyperbolischer Fit $^{85}$Rb')
plt.grid()
plt.xlabel(r'$\text{Amplitude} \mathbin{/} \si{\V}$')
plt.ylabel(r'$T \mathbin{/} \si{\milli\s}$')
plt.legend()
plt.savefig('build/Perioden.pdf')
print(a_p1, b_p1, c_p1)
print(a_p2, b_p2, c_p2)


####dikussion
#Erdmagnetfeld

B_erde_hori = 19.3449 #muT
B_erde_verti = 45.2309 #muT
abw_B_erde_hori1 = abw(B_erde_hori,b1)
abw_B_erde_hori2 = abw(B_erde_hori,b2)
abw_B_erde_verti = abw(B_erde_verti,B_v*10**6)

print(abw_B_erde_hori1,abw_B_erde_hori2)
print(abw_B_erde_verti)


#periodendauern
abw_period = abw(1.5, b_p2/b_p1)
print(abw_period)